import math
import sqlite3
import logging
from cmath import sinh, tanh
from operator import xor

import entsoe.geo.utils
import numpy as np
import pandas as pd
import geopandas as gpd
from pypsa import Network

from _helpers import configure_logging, bus_included

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


def create_case_database(n: Network, components: dict[str, str], external_links: pd.DataFrame, db_raw: str, db: str, config: dict):
    open(db, 'w').close()      # create empty sqlite
    open(db_raw, 'w').close()  # create empty sqlite
    c_raw = sqlite3.connect(db_raw)
    c = sqlite3.connect(db)

    for component in components:
        df = getattr(n, component)
        df.to_sql(component, con=c_raw, dtype={components[component]: 'TEXT PRIMARY KEY'})
        df[config["database"][component]].to_sql(component, con=c, dtype={components[component]: 'TEXT PRIMARY KEY'})

    external_links.to_sql('external_links', c)
    external_links.to_sql('external_links', c_raw)

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


def get_external_links(cross_border_flow: pd.DataFrame, n: Network, config: dict) -> dict[tuple, list[tuple]]:
    external_links = {}
    for link_id, link in n.links.iterrows():
        # external connection added as p-load. Multiple links get aggregated.
        # xor: one of buses is in nordics, one of buses is outside nordics
        bus0 = n.buses.loc[link.bus0]
        bus1 = n.buses.loc[link.bus1]

        bus0_included = bus_included(bus0, config)
        bus1_included = bus_included(bus1, config)

        if xor(bus0_included, bus1_included):
            # one bus is excluded, other is included -> external link
            inside_bus = bus0 if bus0_included else bus1
            outside_bus = bus1 if bus0_included else bus0

            link_key = (inside_bus.bidding_zone, outside_bus.bidding_zone)

            if inside_bus.bidding_zone in cross_border_flow.index and outside_bus.bidding_zone in cross_border_flow.columns:
                if link_key not in external_links:
                    external_links[link_key] = []

                external_links[link_key].append(
                    (inside_bus.name, outside_bus.name)
                )

    return external_links


def insert_external_links(external_links: dict[tuple, list[tuple]], external_flow: pd.DataFrame):
    df = pd.DataFrame(columns=["bus", "bus_outside", "link", "p_set", "q_set"])
    for link in external_links:
        num_connections = len(external_links[link])
        if num_connections > 0:
            flow = external_flow.loc[link[0], link[1]]
            flow = 0 if math.isnan(flow) else flow
            flow_per_connection = flow / num_connections  # evenly distribute over multiple links
            for (bus0, bus1) in external_links[link]:
                d = pd.DataFrame([dict(bus=bus0, bus_outside=bus1, link='-'.join(link), p_set=flow_per_connection, q_set=0)])
                df = pd.concat([df, d], ignore_index=True)

    df.index.names = ['ExternalLink']
    return df


def scale_generation_unit(generation_unit: pd.Series, scale: pd.DataFrame):
    if generation_unit.bidding_zone not in scale_factors:
        return generation_unit.p_set

    return generation_unit.p_set * scale.loc[generation_unit.type, generation_unit.bidding_zone]


def scale_load_per_bidding_zone(n: Network, load_per_bidding_zone: pd.Series):
    n.loads["p_set_country"] = n.loads["p_set"]  # save for comparison

    n.loads["country"] = n.loads["bus"].map(lambda bus: n.buses.loc[bus, "country"])
    n.loads["bidding_zone"] = n.loads["bus"].map(lambda bus: n.buses.loc[bus, "bidding_zone"])

    n.loads["country_factor"] = n.loads.groupby("country", group_keys=False)['p_set'].apply(lambda x: x / x.sum())
    n.loads["bidding_zone_factor"] = n.loads.groupby("bidding_zone", group_keys=False)['country_factor'].apply(lambda x: x / x.sum())

    n.loads["p_set"] = n.loads.groupby("bidding_zone", group_keys=False)['bidding_zone_factor'].apply(
        lambda x: x * (load_per_bidding_zone[x.name] if (x.name in load_per_bidding_zone) else 0)
    )

def remove_load_generation_outside_sync_area(n: Network):
    n.loads.loc[n.loads["bus"].map(lambda bus: n.buses.loc[bus, "in_synchronous_network"]) == 0, 'p_set'] = 0
    n.generators.loc[n.generators["bus"].map(lambda bus: n.buses.loc[bus, "in_synchronous_network"]) == 0, 'p_set'] = 0
    n.storage_units.loc[n.storage_units["bus"].map(lambda bus: n.buses.loc[bus, "in_synchronous_network"]) == 0, 'p_set'] = 0



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

    snapshot = pd.Timestamp(snakemake.config["snapshot"], tz='UTC')

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

    #
    # find and add external links
    #
    logger.info('Finding and adding interconnection links')

    cross_border_flow = pd.read_csv(snakemake.input.cross_border_flows, index_col=0)

    external_links = get_external_links(cross_border_flow, n, snakemake.config)
    external_links = insert_external_links(external_links, cross_border_flow)

    #
    # add link p_set inside nordics
    # !!!! This is manual, if more HVDC links between countries arise, manual fix required !!!!
    #

    n.links["bidding_zone0"] = n.links["bus0"].map(lambda bus: n.buses.loc[bus, "bidding_zone"])
    n.links["bidding_zone1"] = n.links["bus1"].map(lambda bus: n.buses.loc[bus, "bidding_zone"])

    se_fi_links = n.links[(n.links["bidding_zone0"] == 'SE_3') & (n.links["bidding_zone1"] == "FI")]
    se_fi_p_nom = pd.to_numeric(se_fi_links["p_nom"], errors='coerce')
    se_fi_capacity = se_fi_p_nom.sum()
    se_fi_flow = cross_border_flow.loc["FI", "SE_3"]
    n.links.loc[se_fi_links.index, "p_set"] = se_fi_p_nom / se_fi_capacity * se_fi_flow

    #
    # scale generation per bidding_zone and per generation type with entsoe data
    #
    logger.info('Scaling generation')

    n.generators["country"] = n.generators["bus"].map(lambda bus: n.buses.loc[bus, 'country'])
    n.storage_units["country"] = n.storage_units["bus"].map(lambda bus: n.buses.loc[bus, 'country'])

    n.generators["bidding_zone"] = n.generators["bus"].map(lambda bus: n.buses.loc[bus, 'bidding_zone'])
    n.storage_units["bidding_zone"] = n.storage_units["bus"].map(lambda bus: n.buses.loc[bus, 'bidding_zone'])

    all_generation_units = pd.concat([n.generators, n.storage_units])

    # set bidding_zone and in_synchronous_network from bus it's attached to
    all_generation_units["in_synchronous_network"] = all_generation_units["bus"].map(
        lambda bus: n.buses.loc[bus, 'in_synchronous_network'])

    # filter
    generation_units = all_generation_units[all_generation_units["in_synchronous_network"] == 1]
    generation_units = generation_units[generation_units["p_set"] >= snakemake.config["generator_p_min"]]

    # actual generation in sweden is only available for the whole country
    sweden = snakemake.config['bidding_zones_in_country']['SE']
    generation_units.loc[generation_units['bidding_zone'].isin(sweden), 'bidding_zone'] = 'SE'

    # find scale factor per bidding_zone and generation type
    entsoe_generation = pd.read_csv(snakemake.input.actual_generation).set_index('type')
    model_generation = generation_units.groupby(['bidding_zone', 'type'])['p_set'].sum().unstack(level=0)

    scale_factors = (entsoe_generation / model_generation).fillna(0)
    for bidding_zone_in_sweden in sweden:
        scale_factors[bidding_zone_in_sweden] = scale_factors['SE']

    # scale p_set of all generation, including hydro which is in storage_units
    n.generators["p_set"] = n.generators.apply(lambda unit: scale_generation_unit(unit, scale_factors), axis=1)
    n.storage_units["p_set"] = n.storage_units.apply(lambda unit: scale_generation_unit(unit, scale_factors), axis=1)

    #
    # Set p_set to 0 outside sync area
    #
    remove_load_generation_outside_sync_area(n)

    #
    # Scale load per bidding zone
    #
    logger.info('Scaling load per bidding zone')
    entsoe_load = pd.read_csv(snakemake.input.load, index_col=0)['load']
    scale_load_per_bidding_zone(n, entsoe_load)

    #
    # Imbalance per country
    # The imbalance per country is: load - generation + external_flows since load includes losses
    #
    # total_flows = cross_border_flow.sum(axis=1).fillna(0).to_frame('flow')
    # total_flows['country'] = total_flows.index.str[:2]
    # total_flows_country = total_flows.groupby('country')['flow'].sum()
    #
    # generation = n.generators.groupby('country')['p_set'].sum().fillna(0)
    # storage = n.storage_units.groupby('country')['p_set'].sum().fillna(0)
    # all_generation = (generation+storage).fillna(generation)
    # load = n.loads.groupby('country')['p_set'].sum().fillna(0)
    #
    # imbalance = (load + total_flows_country - all_generation).dropna()
    # scale_factors = ((load - imbalance) / load).fillna(1)
    #
    # # scale loads
    # n.loads['p_set'] = n.loads.groupby('country', group_keys=False)['p_set'].apply(lambda x, factor: x * factor[x.name], scale_factors)

    #
    # Load data includes losses, but it shouldn't
    #
    logger.info('Scaling load to compensate for losses')
    scale_factor = (n.loads["p_set"].sum() - snakemake.config["losses"]) / n.loads["p_set"].sum()
    n.loads["p_set"] = scale_factor * n.loads["p_set"]

    # creates the 'raw' database, still saved to keep
    logger.info('Saving to database')
    create_case_database(n, components, external_links, snakemake.output.database_raw, snakemake.output.database, snakemake.config)
