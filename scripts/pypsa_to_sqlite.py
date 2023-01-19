import os
import sqlite3
import logging
from cmath import sqrt, sinh, tanh
from pprint import pprint

import entsoe.geo.utils
import numpy as np
import pandas as pd
import geopandas as gpd
from matplotlib import pyplot as plt
from pypsa import Network

from _helpers import configure_logging
logger = logging.getLogger(__name__)


def map_included_buses(n: Network, config: dict):
    # only include synchronous areas connected to slack_buses
    included_buses = pd.Series(config["build_network"]["slack_buses"])

    size = 0
    while size != included_buses.size:
        size = included_buses.size

        lines = n.lines[n.lines["bus0"].isin(included_buses.tolist()) | n.lines["bus1"].isin(included_buses.tolist())]
        transformers = n.transformers[n.transformers["bus0"].isin(included_buses.tolist()) | n.transformers["bus1"].isin(included_buses.tolist())]

        included_buses = pd.concat([lines["bus0"], lines["bus1"], transformers["bus0"], transformers["bus1"]]).unique()

    n.buses["in_synchronous_network"] = n.buses.index.isin(included_buses)


def set_bidding_zone(n: Network, config: dict, snapshot: pd.Timestamp):
    n.buses['bidding_zone'] = n.buses['country']  # default to country
    for country in config["bidding_zones_in_country"]:
        buses = n.buses.copy()
        buses = buses[buses['country'] == country]
        buses_geo = gpd.GeoDataFrame(buses, geometry=gpd.points_from_xy(buses.x, buses.y, crs='EPSG:4326').to_crs(3857))

        bidding_zones = entsoe.geo.utils.load_zones(config["bidding_zones_in_country"][country], snapshot.tz_convert(None)).to_crs(3857)
        buses_with_bidding_zone = gpd.sjoin_nearest(buses_geo, bidding_zones)
        n.buses.loc[buses_with_bidding_zone.index, "bidding_zone"] = buses_with_bidding_zone["index_right"]


def set_line_impedance_from_linetypes(n: Network, config: dict):
    Z_base = n.lines['v_nom'] ** 2 / config["M_base"]
    num_parallel = n.lines['num_parallel'].replace(0, 1)  # to avoid math problems

    # parameters
    l = n.lines["length"]  # [km]

    r_per_length = n.lines["type"].map(n.line_types.r_per_length) / num_parallel  # [Ohm/km]
    x_per_length = n.lines["type"].map(n.line_types.x_per_length) / num_parallel  # [Ohm/km]
    g_per_length = 0                                                              # [S/km]
    b_per_length = (                                                              # [S/km]
        2 * np.pi * n.lines["type"].map(n.line_types.f_nom)
        * n.lines["type"].map(n.line_types.c_per_length) * 1e-9
        * num_parallel
    )

    z_per_length = r_per_length + 1j*x_per_length  # [Ohm/km]
    y_per_length = g_per_length + 1j*b_per_length  # [S/km]

    z_medium = z_per_length * l
    y_medium = y_per_length * l

    # long line model
    gamma = (z_per_length * y_per_length)**(1/2)  # [1/m]

    f1 = (gamma * l).map(lambda gamma_l: sinh(gamma_l) / gamma_l)            # correction factor [p.u.]
    f2 = (gamma * l).map(lambda gamma_l: tanh(gamma_l / 2) / (gamma_l / 2))  # correction factor [p.u.]

    z_long = z_medium * f1
    y_long = y_medium * f2

    n.lines["r"] = np.real(z_long)
    n.lines["x"] = np.imag(z_long)

    n.lines["g"] = np.real(y_long)
    n.lines["b"] = np.imag(y_long)

    n.lines["r_pu"] = n.lines["r"] / Z_base
    n.lines["x_pu"] = n.lines["x"] / Z_base
    n.lines["g_pu"] = n.lines["g"] * Z_base
    n.lines["b_pu"] = n.lines["b"] * Z_base


def set_dynamic_attributes(n: Network, components: dict[str, str], dynamic_attributes: dict[str, list[str]], snapshot: pd.Timestamp):
    for component in components:
        if component in dynamic_attributes:
            for attribute in dynamic_attributes[component]:
                component_dynamic = getattr(n, component + '_t')
                values = getattr(component_dynamic, attribute)
                getattr(n, component)[attribute].update(values.loc[snapshot])


def set_generator_p_set(n: Network, config: dict):
    p_max_generators = n.generators["p_nom"] * n.generators["p_max_pu"]
    p_max_generators.loc[n.generators["carrier"] == 'nuclear'] = n.generators["p_nom"]  # do not limit nuclear
    p_max_generators[n.generators["p_nom"] < config["generator_p_min"]] = 0

    n.generators["p_set"] = p_max_generators

    p_max_storage = n.storage_units["p_nom"] * n.storage_units["p_max_pu"]
    p_max_storage[n.storage_units["p_nom"] < config["generator_p_min"]] = 0

    n.storage_units["p_set"] = p_max_storage


def create_case_database(n: Network, components: dict[str, str], db_raw: str, db: str, config: dict):
    open(db, 'w').close()      # create empty sqlite
    open(db_raw, 'w').close()  # create empty sqlite
    c_raw = sqlite3.connect(db_raw)
    c = sqlite3.connect(db)

    for component in components:
        df = getattr(n, component)
        df.to_sql(component, con=c_raw, dtype={components[component]: 'TEXT PRIMARY KEY'})
        df[config["database"][component]].to_sql(component, con=c, dtype={components[component]: 'TEXT PRIMARY KEY'})
    c_raw.close()
    c.close()


# !!!! Check these manual corrections when changing config or data gets updated.
# !!!! IDs might not stay the same for instance and corrections are not necessarily required anymore
def apply_parameter_corrections(n: Network, config: dict):
    # map 380 kV to 400 kV
    n.buses["v_nom"].replace(380, 400, inplace=True)
    n.lines["v_nom"].replace(380, 400, inplace=True)

    # SydVästlänken link is not commissioned before ~2022
    n.links.loc['14806', "under_construction"] = True
    n.links.loc['14806', "p_nom"] = True

    # move 3/4 of load from one bus to another to fix loadflow problems
    factor = 1/2
    total_load = n.loads.loc['7010', 'p_set']

    n.loads.loc['7010', 'p_set'] = total_load * (1-factor)

    moved_load = n.loads.loc['7010'].copy()
    moved_load['bus'] = '7017'
    moved_load.p_set = total_load * factor
    n.loads.loc['7017'] = moved_load


def set_generation_type(n: Network, config: dict):
    n.generators["type"] = n.generators["carrier"].map(
        lambda carrier: config["generation"]["pypsa_map"].get(carrier, 'other')
    )

    n.storage_units["type"] = n.storage_units["carrier"].map(
        lambda carrier: config["generation"]["pypsa_map"].get(carrier, 'other')
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("to_sqlite", network="elec")
    configure_logging(snakemake)

    components = {
        'lines': 'Line',
        'transformers': 'Transformer',
        'links': 'Link',
        'generators': 'Generator',
        'loads': 'Load',
        'shunt_impedances': 'ShuntImpedance',
        'storage_units': 'StorageUnit',
        'buses': 'Bus',
    }
    dynamic_attributes = {
        'generators': ['p_max_pu'],
        'loads': ['p_set'],
        'storage_units': ['inflow'],
    }

    n = Network(snakemake.input.network)

    snapshots = pd.read_csv(snakemake.input.snapshots, index_col=0)
    snapshot = snapshots.loc[snakemake.wildcards.case, 'snapshot']
    snapshot = pd.Timestamp(snapshot, tz='UTC')

    logger.info('Setting dynamic attributes for snapshot')
    set_dynamic_attributes(n, components, dynamic_attributes, snapshot)

    logger.info('Applying manual parameter corrections')
    apply_parameter_corrections(n, snakemake.config)

    logger.info('Adding in_synchronous_area data')
    map_included_buses(n, snakemake.config)

    logger.info('Adding bidding_zone info (for DK)')
    set_bidding_zone(n, snakemake.config, snapshot)
    logger.info('Setting line impedance')
    set_line_impedance_from_linetypes(n, snakemake.config)
    logger.info('Setting generation p_set')
    set_generator_p_set(n, snakemake.config)
    logger.info('Setting generation type')
    set_generation_type(n, snakemake.config)

    # creates the 'raw' database, still saved to keep
    logger.info('Saving to database')
    create_case_database(n, components, snakemake.output.database_raw, snakemake.output.database, snakemake.config)
