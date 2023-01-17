import logging
import math
import shutil
import sqlite3
from operator import xor
from sqlite3 import Connection

import pandas as pd

from _helpers import configure_logging, bus_included

logger = logging.getLogger(__name__)


def get_external_links(external_flow: pd.Series, buses: pd.DataFrame, links: pd.DataFrame, config: dict) -> dict[str, tuple]:
    external_links = external_flow.index.values.tolist()
    external_links = {k: [] for k in external_links}

    logger.info('Finding interconnection links')
    for bus_id, bus in buses.iterrows():
        links_from_bus = links[links['bus0'] == bus_id]
        for link_id, link in links_from_bus.iterrows():
            # external connection added as p-load. Multiple links get aggregated.
            # xor: one of buses is in nordics, one of buses is outside nordics
            bus0 = buses.loc[link.bus0]
            bus1 = buses.loc[link.bus1]

            bus0_included = bus_included(bus0, config)
            bus1_included = bus_included(bus1, config)
            if xor(bus0_included, bus1_included):
                inside_bus = bus0 if bus0_included else bus1
                outside_bus = bus1 if bus0_included else bus0

                link_name = inside_bus.country + '-' + outside_bus.country

                if link_name in external_links:
                    external_links[inside_bus.country + '-' + outside_bus.country].append(
                        (inside_bus.name, outside_bus.name)
                    )
    return external_links


def insert_external_links(c: Connection, external_links: dict[str, tuple], external_flow: pd.Series):
    df = pd.DataFrame(columns=["bus", "bus_outside", "link", "p_set", "q_set"])
    for link in external_links:
        flow = external_flow.loc[link]
        flow = 0 if math.isnan(flow) else flow
        num_connections = len(external_links[link])
        if num_connections > 0:
            flow_per_connection = flow / num_connections  # evenly distribute over multiple links
            for (bus0, bus1) in external_links[link]:
                d = pd.DataFrame([dict(bus=bus0, bus_outside=bus1, link=link, p_set=flow_per_connection, q_set=0)])
                df = pd.concat([df, d], ignore_index=True)

    df.index.names = ['ExternalLink']
    df.to_sql('external_links', c)


def scale_generation_unit(generation_unit: pd.Series, scale: pd.DataFrame):
    scale_factor = scale.loc[generation_unit.type, generation_unit.country]
    scale_factor = 1 if math.isnan(scale_factor) else scale_factor  # dont scale NaN
    return generation_unit.p_set * scale_factor


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake('external_links_to_sqlite', case='high')
    configure_logging(snakemake)

    # copy database and connect
    shutil.copy2(snakemake.input.database, snakemake.output.database)
    c = sqlite3.connect(snakemake.output.database)

    snapshots = pd.read_csv(snakemake.input.snapshots, index_col=0)
    snapshot = snapshots.loc[snakemake.wildcards.case, 'snapshot']

    #
    # find and add external links
    #

    external_flows = pd.read_csv(snakemake.input.cross_border_flows, index_col=0)
    external_flows.index = pd.to_datetime(external_flows.index)
    external_flow = external_flows.loc[snapshot]

    buses = pd.read_sql('SELECT * FROM buses', c).set_index('Bus')
    links = pd.read_sql('SELECT * FROM links', c).set_index('Link')

    external_links = get_external_links(external_flow, buses, links, snakemake.config)
    insert_external_links(c, external_links, external_flow)

    #
    # add link p_set inside nordics
    # !!! This is manual, if more HVDC links between countries arise, manual fix required !!!
    #

    links["country0"] = links["bus0"].map(lambda bus: buses.loc[bus, "country"])
    links["country1"] = links["bus1"].map(lambda bus: buses.loc[bus, "country"])

    se_fi_links = links[(links["country0"] == 'SE') & (links["country1"] == "FI")]
    se_fi_p_nom = pd.to_numeric(se_fi_links["p_nom"], errors='coerce')
    se_fi_capacity = se_fi_p_nom.sum()
    se_fi_flow = external_flow["FI-SE_3"]
    links.loc[se_fi_links.index, "p_set"] = se_fi_p_nom / se_fi_capacity * se_fi_flow
    links.to_sql('links', con=c, if_exists='replace')

    #
    # scale generation per country and per generation type with entsoe data
    #

    # get generation from database
    generators = pd.read_sql('SELECT * FROM generators', c).set_index('Generator')
    storage_units = pd.read_sql('SELECT * FROM storage_units', c).set_index('StorageUnit')

    generators["country"] = generators["bus"].map(lambda bus: buses.loc[bus, 'country'])
    storage_units["country"] = storage_units["bus"].map(lambda bus: buses.loc[bus, 'country'])

    all_generation_units = pd.concat([generators, storage_units])

    # set country and in_synchronous_network from bus it's attached to
    all_generation_units["in_synchronous_network"] = all_generation_units["bus"].map(lambda bus: buses.loc[bus, 'in_synchronous_network'])

    # filter
    generation_units = all_generation_units[all_generation_units["in_synchronous_network"] == 1]
    generation_units = generation_units[generation_units["p_set"] >= snakemake.config["generator_p_min"]]

    # find scale factor per country and generation type
    entsoe_generation = pd.read_csv(snakemake.input.actual_generation).set_index('type')
    model_generation = generation_units.groupby(['country', 'type'])['p_set'].sum().unstack(level=0)

    scale_factors = entsoe_generation / model_generation

    # scale p_set of all generation, including hydro which is in storage_units
    generators["p_set"] = generators.apply(lambda unit: scale_generation_unit(unit, scale_factors), axis=1)
    storage_units["p_set"] = storage_units.apply(lambda unit: scale_generation_unit(unit, scale_factors), axis=1)

    #
    # save to database
    #
    generators[snakemake.config["database"]["generators"]].to_sql('generators', con=c, dtype={'Generator': 'TEXT PRIMARY KEY'}, if_exists='replace')
    storage_units[snakemake.config["database"]["storage_units"]].to_sql('storage_units', con=c, dtype={'StorageUnit': 'TEXT PRIMARY KEY'}, if_exists='replace')

    c.close()
