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


#Matches the geotag of each bus with the coordinates of a bidding zone
def set_bidding_zone(n: Network, config: dict, snapshot: pd.Timestamp):
    n.buses['bidding_zone'] = n.buses['country']  # default to country
    for country in config["bidding_zones_in_country"]:
        buses = n.buses.copy()
        buses = buses[buses['country'] == country]
        buses_geo = gpd.GeoDataFrame(buses, geometry=gpd.points_from_xy(buses.x, buses.y, crs='EPSG:4326').to_crs(3857))

        bidding_zones = entsoe.geo.utils.load_zones(config["bidding_zones_in_country"][country], snapshot.tz_convert(None)).to_crs(3857)
        buses_with_bidding_zone = gpd.sjoin_nearest(buses_geo, bidding_zones)
        n.buses.loc[buses_with_bidding_zone.index, "bidding_zone"] = buses_with_bidding_zone["index_right"]

#merges existing database (nordics_db) with the database that includes a region for each bus (example_db)

def merge_buses(nordics_db, example_db, merged_db):
    # Create a SQLite database connection for nordics.sqlite
    nordics_conn = sqlite3.connect(nordics_db)

    # Read the "buses" table from nordics.sqlite into a pandas DataFrame
    nordics_buses_df = pd.read_sql_query("SELECT * from buses", nordics_conn)

    # Close the nordics.sqlite connection
    nordics_conn.close()

    # Create a SQLite database connection for example.sqlite
    example_conn = sqlite3.connect(example_db)

    # Read the "Buses" table from example.sqlite into a pandas DataFrame
    example_buses_df = pd.read_sql_query("SELECT * from Buses", example_conn)

    # Close the example.sqlite connection
    example_conn.close()

    example_buses_df['Bus'] = example_buses_df['Bus'].astype(str)
    example_buses_df = example_buses_df.drop(columns=['x', 'y'], axis=1)

    # Merge the two DataFrames based on the "buses" column
    merged_buses_df = pd.merge(nordics_buses_df, example_buses_df, on="Bus", how="left")

    # Create a SQLite database connection for the merged database
    merged_conn = sqlite3.connect(merged_db)

    # Save the merged DataFrame to a new table in the merged database
    merged_buses_df.to_sql('buses', merged_conn, if_exists='replace', index=False)

    # Close the merged database connection
    merged_conn.close()


def set_line_impedance_from_linetypes(n: Network, config: dict):
    for i, line in n.lines.iterrows():
        ll = int(line.iloc[4])
        v_nom = float(line["v_nom"])
        num_parallel = n.lines.iloc[:, 3].astype(int)

        if (num_parallel == 0).any():
            num_parallel[num_parallel == 0] = 1

        if (ll > 250) and v_nom == 400.0:
            num_parallel = n.lines.loc[i, "num_parallel"]
            if num_parallel == 0:
                num_parallel = 1

            Z_base = n.lines.loc[i, 'v_nom'] ** 2 / config["M_base"]
            r_per_length = 0.0326 / num_parallel
            x_per_length = 0.27 / num_parallel
            n.lines.loc[i, "s_nom"] = 2245 * 0.7
            n.lines.loc[i, "s_nom2"] = 2245

            b_per_length = (  # [S/km]
                    2 * np.pi * 50
                    * 13.22 * 1e-9
                    * num_parallel
            )
            g_per_length = 0

            z_per_length = r_per_length + 1j * x_per_length  # [Ohm/km]
            y_per_length = g_per_length + 1j * b_per_length  # [S/km]

            z_medium = z_per_length * ll
            y_medium = y_per_length * ll

            gamma = (z_per_length * y_per_length)**(1/2)  # [1/m]

            f1 = sinh(gamma * ll) / (gamma * ll)
            f2 = tanh((gamma * ll) / 2) / ((gamma * ll) / 2)

            z_long = z_medium * f1
            y_long = y_medium * f2

            r = np.real(z_long)
            x = np.imag(z_long)

            g = np.real(y_long)
            b = np.imag(y_long)

            n.lines.loc[i, "r_pu"] = r / Z_base
            n.lines.loc[i, "x_pu"] = x / Z_base
            n.lines.loc[i, "g_pu"] = g * Z_base
            n.lines.loc[i, "b_pu"] = b * Z_base


        elif v_nom == 400.00:
            num_parallel = n.lines.loc[i, "num_parallel"]
            if num_parallel == 0:
                num_parallel = 1
            Z_base = n.lines.loc[i, 'v_nom'] ** 2 / config["M_base"]
            r_per_length = 0.0326 / num_parallel
            x_per_length = 0.27 / num_parallel
            n.lines.loc[i, "s_nom"] = 2245 * 0.7
            n.lines.loc[i, "s_nom2"] = 2245

            b_per_length = (  # [S/km]
                    2 * np.pi * 50
                    * 13.22 * 1e-9
                    * num_parallel
            )
            g_per_length = 0  # [S/km]
            z_per_length = r_per_length + 1j * x_per_length  # [Ohm/km]
            y_per_length = g_per_length + 1j * b_per_length  # [S/km]

            z_medium = z_per_length * ll
            y_medium = y_per_length * ll

            r = np.real(z_medium)
            x = np.imag(z_medium)

            g = np.real(y_medium)
            b = np.imag(y_medium)

            n.lines.loc[i, "r_pu"] = r / Z_base
            n.lines.loc[i, "x_pu"] = x / Z_base
            n.lines.loc[i, "g_pu"] = g * Z_base
            n.lines.loc[i, "b_pu"] = b * Z_base
        elif v_nom == 301.00:
            num_parallel = n.lines.loc[i, "num_parallel"]
            if num_parallel == 0:
                num_parallel = 1
            Z_base = n.lines.loc[i, 'v_nom'] ** 2 / config["M_base"]
            r_per_length = 0.0477 / num_parallel
            x_per_length = 0.42 / num_parallel
            n.lines.loc[i, "s_nom"] = 450 * 2 * 0.7
            n.lines.loc[i, "s_nom2"] = 900

            b_per_length = (  # [S/km]
                    2 * np.pi * 50
                    * 8.62 * 1e-9
                    * num_parallel
            )
            g_per_length = 0  # [S/km]
            z_per_length = r_per_length + 1j * x_per_length  # [Ohm/km]
            y_per_length = g_per_length + 1j * b_per_length  # [S/km]

            z_medium = z_per_length * ll
            y_medium = y_per_length * ll

            r = np.real(z_medium)
            x = np.imag(z_medium)

            g = np.real(y_medium)
            b = np.imag(y_medium)

            n.lines.loc[i, "r_pu"] = r / Z_base
            n.lines.loc[i, "x_pu"] = x / Z_base
            n.lines.loc[i, "g_pu"] = g * Z_base
            n.lines.loc[i, "b_pu"] = b * Z_base

        elif ll == 1:
            n.lines.loc[i, "r_pu"] = 0.0001
            n.lines.loc[i, "x_pu"] = 0.0001
            n.lines.loc[i, "g_pu"] = 0.0001
            n.lines.loc[i, "b_pu"] = 0.0001
            n.lines.loc[i, "s_nom"] = 3000
            n.lines.loc[i, "s_nom2"] = 3000

        else:
            num_parallel = n.lines.loc[i, "num_parallel"]
            if num_parallel == 0:
                num_parallel = 1
            Z_base = n.lines.loc[i, 'v_nom'] ** 2 / config["M_base"]
            r_per_length = 0.0384 / num_parallel
            x_per_length = 0.30 / num_parallel
            n.lines.loc[i, "s_nom"] = 746 * 0.7
            n.lines.loc[i, "s_nom2"] = 746

            b_per_length = (  # [S/km]
                    2 * np.pi * 50
                    * 11.84 * 1e-9
                    * num_parallel
            )
            g_per_length = 0  # [S/km]
            z_per_length = r_per_length + 1j * x_per_length  # [Ohm/km]
            y_per_length = g_per_length + 1j * b_per_length  # [S/km]

            z_medium = z_per_length * ll
            y_medium = y_per_length * ll

            r = np.real(z_medium)
            x = np.imag(z_medium)

            g = np.real(y_medium)
            b = np.imag(y_medium)

            n.lines.loc[i, "r_pu"] = r / Z_base
            n.lines.loc[i, "x_pu"] = x / Z_base
            n.lines.loc[i, "g_pu"] = g * Z_base
            n.lines.loc[i, "b_pu"] = b * Z_base
    i = 0


def set_dynamic_attributes(n: Network, components: dict[str, str], dynamic_attributes: dict[str, list[str]], snapshot: pd.Timestamp):
    for component in components:
        if component in dynamic_attributes:
            for attribute in dynamic_attributes[component]:
                component_dynamic = getattr(n, component + '_t')
                values = getattr(component_dynamic, attribute)
                getattr(n, component)[attribute].update(values.loc[snapshot])


#do not change here, change in def add_transformer in build_psse_network.py
def set_transformer_x_pu(n: Network):
    n.transformers['x_pu'] = n.transformers['x']
    n.transformers['x'] = None


#sets p_set to p_nom before scaling p_set in def scale_generation_per_bidding_zone
def set_generator_p_set(n: Network, config: dict):
    p_max_generators = n.generators["p_nom"] * n.generators["p_max_pu"]
    p_max_generators.loc[n.generators["carrier"] == 'nuclear'] = n.generators["p_nom"]  # do not limit nuclear
    p_max_generators[n.generators["p_nom"] < config["generator_p_min"]] = 0

    n.generators["p_set"] = p_max_generators
    n.generators["p_set_before_scaling"] = p_max_generators

    p_max_storage = n.storage_units["p_nom"] * n.storage_units["p_max_pu"]
    p_max_storage[n.storage_units["p_nom"] < config["generator_p_min"]] = 0

    n.storage_units["p_set"] = p_max_storage
    n.storage_units["p_set_before_scaling"] = p_max_storage


#creates and fills the database
#the raw database includes all information avaliable, the regular only includes what is used to build the model in psse
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


# changes some inaccuracies from the pypsa model. The load distribution is in the north of norway/finland
# !!!! Check these manual corrections when changing config or data gets updated.
# !!!! IDs might not stay the same for instance and corrections are not necessarily required anymore
def apply_parameter_corrections(n: Network, config: dict):
    # map 380 kV to 400 kV
    n.buses["v_nom"].replace(380, 400, inplace=True)
    n.lines["v_nom"].replace(380, 400, inplace=True)

    n.buses["v_nom"].replace(300, 400, inplace=True)
    n.lines["v_nom"].replace(300, 400, inplace=True)

    # move part of load from one bus to another to fix loadflow problems
    factor = 7/10
    total_load = n.loads.loc['7010', 'p_set']

    #norge-finland
    n.loads.loc['7010', 'p_set'] = total_load * (1-factor)

    moved_load = n.loads.loc['7010'].copy()
    moved_load['bus'] = '7017'
    moved_load.p_set = total_load * factor
    n.loads.loc['7017'] = moved_load


#redistributes load from low voltage side of transformer to the high voltage side
def redistribute_load(n: Network, config: dict):
    factor = 3/4
    df = pd.read_excel('Transformer/transformer_list.xlsx')
    word_to_remove = "HYDRO"
    df = df[~df['From Bus  Name'].str.contains(word_to_remove)]
    df = df.reset_index(drop=True)
    df = df.drop(['From Bus  Name', 'To Bus  Name'], axis=1)

    for index, row in df.iterrows():
        from_string = str(row['From Bus  Number'])
        to_string = str(row['To Bus  Number'])
        total_load = n.loads.loc[from_string, 'p_set']
        n.loads.loc[from_string, 'p_set'] = total_load * (1-factor)
        moved_load = n.loads.loc[from_string].copy()
        moved_load['bus'] = to_string
        moved_load.p_set = total_load * factor
        n.loads.loc[to_string] = moved_load


# corrects generation type
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


#to control the NO_2-GB flow, change the value in the NO_2xNO_2 cell in the cross border flow sheet
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


#Here the hvdc links within the nordics are set.
def set_link_p_set(n: Network, cross_border_flow: pd.DataFrame):
    n.links["bidding_zone0"] = n.links["bus0"].map(lambda bus: n.buses.loc[bus, "bidding_zone"])
    n.links["bidding_zone1"] = n.links["bus1"].map(lambda bus: n.buses.loc[bus, "bidding_zone"])

    se_fi_links = n.links[(n.links["bidding_zone0"] == 'SE_3') & (n.links["bidding_zone1"] == "FI")]
    se_fi_p_nom = pd.to_numeric(se_fi_links["p_nom"], errors='coerce')
    se_fi_capacity = se_fi_p_nom.sum()
    se_fi_flow = cross_border_flow.loc["FI", "SE_3"]
    n.links.loc[se_fi_links.index, "p_set"] = se_fi_p_nom / se_fi_capacity * se_fi_flow

    # SydVästlänken link is not commissioned before ~2022
    n.links.loc['14806', "under_construction"] = False
    n.links.loc['14806', "p_nom"] = 1200
    n.links.loc['14806', 'p_set'] = -200


def scale_generation_per_bidding_zone(n: Network, bz_generation: pd.DataFrame, generators_case):
    n.generators = n.generators.drop('6937 onwind', axis=0)
    n.generators = n.generators.drop('7008 onwind', axis=0)
    n.generators = n.generators.drop('7010 onwind', axis=0)
    n.generators = n.generators.drop('7019 onwind', axis=0)
    n.generators = n.generators.drop('6655 onwind', axis=0)
    n.generators = n.generators.drop('6660 onwind', axis=0)
    n.generators = n.generators.drop('6672 onwind', axis=0)
    n.generators = n.generators.drop('6716 onwind', axis=0)
    n.generators = n.generators.drop('6744 onwind', axis=0)
    n.generators = n.generators.drop('6746 onwind', axis=0)
    n.generators = n.generators.drop('6747 onwind', axis=0)
    n.generators = n.generators.drop('6748 onwind', axis=0)
    #n.generators = n.generators.drop('Loviisa-1', axis=0)
    #n.generators = n.generators.drop('Loviisa-2 2', axis=0)

    all_generation_units = pd.concat([n.generators, n.storage_units])

    # set bidding_zone and in_synchronous_network from bus it's attached to
    all_generation_units["in_synchronous_network"] = all_generation_units["bus"].map(
        lambda bus: n.buses.loc[bus, 'in_synchronous_network'])

    # add manual generator included in the scaling
    all_generation_units = pd.concat([all_generation_units, generators_case])
    n.generators = pd.concat([n.generators, generators_case])
    n.generators.index.names = ['Generator']

    # filter
    generation_units = all_generation_units[all_generation_units["in_synchronous_network"] == 1]
    generation_units = generation_units[generation_units["p_set"] > 0]

    # find scale factor per bidding_zone and generation type
    model_generation = generation_units.groupby(['bidding_zone', 'type'])['p_set'].sum().unstack(level=0)
    model_generation = model_generation.fillna(0)

    bz_generation = bz_generation.fillna(0)
    scale_factors = bz_generation.set_index('type').div(model_generation)

    n.generators["p_set"] = n.generators.apply(lambda unit: scale_generation_unit(unit, scale_factors), axis=1)
    n.storage_units["p_set"] = n.storage_units.apply(lambda unit: scale_generation_unit(unit, scale_factors), axis=1)


def scale_generation_unit(generation_unit: pd.Series, scale: pd.DataFrame):
    if generation_unit.bidding_zone not in scale:
        return generation_unit.p_set

    return generation_unit.p_set * scale.loc[generation_unit.type, generation_unit.bidding_zone]


# this uses the pypsa way of distributing the loads
def scale_load_per_bidding_zone(n: Network, load_per_bidding_zone: pd.Series):
    n.loads["p_set_country"] = n.loads["p_set"]  # save for comparison

    n.loads["country"] = n.loads["bus"].map(lambda bus: n.buses.loc[bus, "country"])
    n.loads["bidding_zone"] = n.loads["bus"].map(lambda bus: n.buses.loc[bus, "bidding_zone"])

    n.loads["country_factor"] = n.loads.groupby("country", group_keys=False)['p_set'].apply(lambda x: x / x.sum())
    n.loads["bidding_zone_factor"] = n.loads.groupby("bidding_zone", group_keys=False)['country_factor'].apply(lambda x: x / x.sum())

    n.loads["p_set"] = n.loads.groupby("bidding_zone", group_keys=False)['bidding_zone_factor'].apply(
        lambda x: x * (load_per_bidding_zone[x.name] if (x.name in load_per_bidding_zone) else 0)
    )

    dk_2 = n.loads.loc[n.loads['bidding_zone'] == 'DK_2']
    total_dk2 = dk_2["p_set_country"].sum()

    # update only DK2 because it's different
    n.loads.loc[n.loads["bidding_zone"] == 'DK_2', "p_set"] = n.loads["p_set_country"] / total_dk2 * load_per_bidding_zone['DK_2']

    #scalar bidding zones
    n.loads.loc[n.loads["bidding_zone"] == "NO_1", "p_set"] = n.loads["p_set"]

    n.loads.loc[n.loads["bidding_zone"] == "NO_5", "p_set"] = n.loads["p_set"]


def remove_load_generation_outside_sync_area(n: Network):
    n.loads.loc[n.loads["bus"].map(lambda bus: n.buses.loc[bus, "in_synchronous_network"]) == 0, 'p_set'] = 0
    try:
        n.generators.loc[n.generators["bus"].map(lambda bus: n.buses.loc[bus, "in_synchronous_network"]) == 0, 'p_set'] = 0
    except KeyError:
        pass
    n.storage_units.loc[n.storage_units["bus"].map(lambda bus: n.buses.loc[bus, "in_synchronous_network"]) == 0, 'p_set'] = 0


# removes losses from loads per bidding zone defined in the config.yaml file
def subtract_losses_from_loads(n: Network, config: dict):
    losses = pd.Series(config["losses"])
    load_per_bidding_zone = n.loads.groupby('bidding_zone')["p_set"].sum()

    scale_factor = (load_per_bidding_zone - losses) / load_per_bidding_zone

    n.loads["scale_factor_losses"] = n.loads["bidding_zone"].map(
        lambda bidding_zone: scale_factor.get(bidding_zone, 1.0))
    n.loads["p_set_with_losses"] = n.loads["p_set"]
    n.loads["p_set"] = n.loads["scale_factor_losses"] * n.loads["p_set_with_losses"]


# adds loads to some of the new buses in the model. More can be added if necessary
def add_new_loads(n):

    #häradsbo-alvesta
    factor = 3/10

    total_load_1 = n.loads.loc['6349', 'p_set']
    n.loads.loc['6349', 'p_set'] = total_load_1 * (1-factor)

    moved_load_1 = n.loads.loc['6349'].copy()
    moved_load_1['bus'] = '1053'
    moved_load_1.p_set = total_load_1 * factor
    n.loads.loc['1053'] = moved_load_1

    #breared-söderåsen
    factor = 5/10
    total_load = n.loads.loc['6335', 'p_set']
    n.loads.loc['6335', 'p_set'] = total_load * (1-factor)

    moved_load = n.loads.loc['6335'].copy()
    moved_load['bus'] = '1052'
    moved_load.p_set = total_load * factor
    n.loads.loc['1052'] = moved_load

    #sauda-forre-lyse
    factor = 5/10
    total_load = n.loads.loc['6694', 'p_set']
    n.loads.loc['6694', 'p_set'] = total_load * (1-factor * (3/2))

    moved_load = n.loads.loc['6694'].copy()
    moved_load['bus'] = '6682'
    moved_load.p_set = total_load * factor / 2
    n.loads.loc['6682'] = moved_load

    moved_load['bus'] = '6681'
    moved_load.p_set = total_load * factor / 2
    n.loads.loc['6681'] = moved_load

    moved_load['bus'] = '6698'
    moved_load.p_set = total_load * factor / 2
    n.loads.loc['6698'] = moved_load

    #finland
    factor = 9/10

    total_load = n.loads.loc['6903', 'p_set']
    n.loads.loc['6903', 'p_set'] = total_load * (1-factor)

    moved_load = n.loads.loc['6903'].copy()
    moved_load['bus'] = '6920'
    moved_load.p_set = total_load * factor
    n.loads.loc['6920'] = moved_load


def add_bus_names(n, bus_names):
    bus_names['Bus'] = bus_names['Bus'].astype(n.buses.index.dtype)
    n.buses = n.buses.merge(bus_names, on='Bus', how='left')

    bus_names['Bus'] = bus_names['Bus'].astype(n.buses.index.dtype)
    bus_names['Bus'] = bus_names['Bus'].astype(str)

    for index, row in n.buses.iterrows():
        for i, bus_row in bus_names.iterrows():
            n.buses['Bus'] = n.buses['Bus'].astype(int)
            if row['Bus'] == bus_row['Bus']:
                n.buses.loc[index, 'symbol'] = bus_row['symbol_name']
                break

    n.buses=n.buses.set_index('Bus')[[ 'v_nom', 'symbol', 'under_construction', 'tags', 'x', 'y',
       'carrier', 'country', 'lv', 'onshore', 'has_connections',
       'substation_lv', 'generator_substation_lv', 'load_substation_lv',
       'substation_off', 'base_only', 'type', 'unit', 'v_mag_pu_set',
       'v_mag_pu_min', 'v_mag_pu_max', 'control', 'sub_network',
       'in_synchronous_network', 'bidding_zone', 'symbol_name']]


def add_manual_buses(n, manual_buses):
    n.buses = pd.concat([n.buses, manual_buses])
    n.buses.index.names = ['Bus']


#Adding a manual load
def add_manual_loads(n, manual_loads):
    n.loads = pd.concat([n.loads, manual_loads])
    n.loads.index.names = ['Load']


def add_manual_shunts(n, manual_shunts):
    n.shunt_impedances = pd.concat([n.shunt_impedances, manual_shunts])
    n.shunt_impedances.index.names = ['ShuntImpedance']
    n.shunt_impedances['bus']= n.shunt_impedances['bus'].astype(str)


def add_manual_generators(n, manual_generators):
    n.generators = pd.concat([n.generators, manual_generators])
    n.generators.index.names = ['Generator']


def add_manual_lines(n, manual_lines, case_lines):
    n.lines = pd.concat([n.lines, manual_lines])
    n.lines.index.names = ['Line']
    n.lines = pd.concat([n.lines, case_lines])
    n.lines.index.names = ['Line']


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
    #py-psa snapshot above

    logger.info('Setting x_pu on transformers')
    set_transformer_x_pu(n)

    logger.info('Applying manual parameter corrections')
    apply_parameter_corrections(n, snakemake.config)

    logger.info('Redistributing load from 220kV to 400kV')
    redistribute_load(n, snakemake.config)

    logger.info('Adding in_synchronous_area data')
    map_included_buses(n, snakemake.config)

    logger.info('Adding bidding_zone info')
    set_bidding_zone(n, snakemake.config, snapshot)

    # set country and biddingzone columns
    n.generators["country"] = n.generators["bus"].map(lambda bus: n.buses.loc[bus, 'country'])
    n.storage_units["country"] = n.storage_units["bus"].map(lambda bus: n.buses.loc[bus, 'country'])

    n.generators["bidding_zone"] = n.generators["bus"].map(lambda bus: n.buses.loc[bus, 'bidding_zone'])
    n.storage_units["bidding_zone"] = n.storage_units["bus"].map(lambda bus: n.buses.loc[bus, 'bidding_zone'])

    logger.info('Adding manual buses')
    manual_buses = pd.read_excel(snakemake.input.manual_buses, index_col=0)
    add_manual_buses(n, manual_buses)

    logger.info('Adding manual lines')
    manual_lines = pd.read_excel(snakemake.input.manual_lines, index_col=0)
    case_lines = pd.read_excel(snakemake.input.case_lines, index_col=0)
    add_manual_lines(n, manual_lines, case_lines)

    logger.info('Setting line impedance')
    set_line_impedance_from_linetypes(n, snakemake.config)

    logger.info('Setting generation type')
    set_generation_type(n, snakemake.config)

    logger.info('Setting generation p_set')
    set_generator_p_set(n, snakemake.config)

    logger.info('Finding and adding interconnection links')
    cross_border_flow = pd.read_csv(snakemake.input.cross_border_flows, index_col=0)
    external_links = get_external_links(cross_border_flow, n, snakemake.config)
    external_links = insert_external_links(external_links, cross_border_flow)

    # !!!! This is manual, if more HVDC links between countries or biddingzones arise, manual fix required !!!!
    logger.info('Setting p_set for inter-area HVDC-links')
    set_link_p_set(n, cross_border_flow)

    logger.info('Scaling generation per bidding zone')
    bz_generation = pd.read_excel(snakemake.input.actual_generation)
    generators_case = pd.read_excel(snakemake.input.generators_case, index_col=0)
    scale_generation_per_bidding_zone(n, bz_generation, generators_case)

    logger.info('Removing generation and load outside Nordic area')
    remove_load_generation_outside_sync_area(n)

    logger.info('Scaling load per bidding zone')
    entsoe_load = pd.read_csv(snakemake.input.load, index_col=0)['load']
    scale_load_per_bidding_zone(n, entsoe_load)

    logger.info('Scaling load to compensate for losses')
    subtract_losses_from_loads(n, snakemake.config)

    logger.info('Moving loads')
    add_new_loads(n)

    logger.info('Adding manual loads')
    manual_loads = pd.read_excel(snakemake.input.manual_loads, index_col=0)
    add_manual_loads(n, manual_loads)

    logger.info('Adding manual generators')
    manual_generators = pd.read_excel(snakemake.input.manual_generators, index_col=0)
    add_manual_generators(n, manual_generators)

    logger.info('Adding manual shunts')
    manual_shunts = pd.read_excel(snakemake.input.manual_shunts, index_col=0)
    add_manual_shunts(n, manual_shunts)

    # if no bus names are needed, this can be commented out, it takes a long time to run
    logger.info('Adding bus names')
    bus_names = pd.read_excel(snakemake.input.bus_names)
    add_bus_names(n, bus_names)

    # creates the database with useful information and a 'raw' database with more information
    logger.info('Saving to database')
    create_case_database(n, components, external_links, snakemake.output.database_raw, snakemake.output.database, snakemake.config)

    # here you need to make sure that the database/INSERT_OWN/nordics is matching your case name
    logger.info('Adding region info')
    merge_buses("database/ren_hh/nordics.sqlite", 'example.sqlite', "database/ren_hh/nordics.sqlite")