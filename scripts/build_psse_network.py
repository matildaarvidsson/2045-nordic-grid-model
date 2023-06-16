"""
Creates a PSS/E model from a sqlite database.
"""
import logging
import math
import sqlite3

import pandas as pd

from _helpers import configure_logging, bus_included

import psse35
import psspy
from psspy import _i, _f, _s

logger = logging.getLogger(__name__)


def add_area(id: int, name: str):
    psspy.area_data(id, _i, [_f, _f], name)


def add_region(id: int, LnNamn: str):
    psspy.zone_data(id, LnNamn)


def add_bus(bus: pd.Series, country_id: int, LnNamn_id: int, has_generation_units: bool, config: dict):
    if bus.carrier == 'DC' or bus.v_nom is None:
        return

    if bus.name in config["build_network"]["slack_buses"]:
        code = 3  # swing bus
    elif has_generation_units:
        code = 2  # generator bus
    else:
        code = 1  # non-gen bus

    try:
        psspy.bus_data_4(
            int(bus.name),
            0,
            [code, country_id, LnNamn_id, _i],
            [bus.v_nom, _f, _f, _f, _f, _f, _f],
            bus.symbol[0:12],
        )
    except Exception as e:
        logger.info(bus)
        raise e

    psspy.bus_chng_4(6875, 0, [2, _i, _i, _i], [301.0, _f, _f, _f, _f, _f, _f], _s)
    psspy.bus_chng_4(6667, 0, [2, _i, _i, _i], [301.0, _f, _f, _f, _f, _f, _f], _s)
    psspy.bus_chng_4(6671, 0, [2, _i, _i, _i], [301.0, _f, _f, _f, _f, _f, _f], _s)
    psspy.bus_chng_4(6751, 0, [2, _i, _i, _i], [301.0, _f, _f, _f, _f, _f, _f], _s)
    psspy.bus_chng_4(6699, 0, [2, _i, _i, _i], [301.0, _f, _f, _f, _f, _f, _f], _s)
    psspy.bus_chng_4(6706, 0, [2, _i, _i, _i], [301.0, _f, _f, _f, _f, _f, _f], _s)
    psspy.bus_chng_4(6707, 0, [2, _i, _i, _i], [301.0, _f, _f, _f, _f, _f, _f], _s)


def add_machine(machine: pd.Series, bus: pd.Series, id: int, in_service: bool, config: dict) -> bool:
    # TODO r_source, x_source, r_tran, x_tran
    # TODO machines at DC-buses
    if bus.carrier == 'DC':
        return False

    p = machine.p_set

    if p <= 0:
        return False

    pf = 0.94  # Check value, maybe change per generation type ?
    angle = math.acos(pf)
    m_base = machine.p_nom / pf  # do not take into account setpoint, only actual capacity [MW]

    p_min = 0
    p_max = machine.p_nom

    #q_max = m_base * math.sin(angle)
    q_max = 9999
    q_min = -q_max

    try:
        psspy.plant_data_4(
            int(machine.bus),
            0,
            [_i, _i],
            [config["build_network"]["scheduled_voltage"], _f]
        )

        psspy.machine_data_4(
            int(machine.bus),
            str(id),
            [int(in_service), _i, _i, _i, _i, _i, _i],
            [p, 0, q_max, q_min, p_max, p_min, m_base, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f],
            _s,
        )
    except Exception as e:
        logger.info(machine)
        raise e

    return True



def add_load(load: pd.Series, bus: pd.Series, id: int, in_service: bool):
    # TODO: validate data
    try:
        add_load_psse(int(load.bus), str(id), load.p_set, load.q_set, in_service=in_service, scalable=True)
    except Exception as e:
        logger.info(load)
        raise e


def add_external_link(external_link: pd.Series, bus: pd.Series, id: int, in_service: bool, area: int):
    try:
        add_load_psse(int(external_link.bus), 'E'+str(id), external_link.p_set, external_link.q_set, external_link.link, in_service=in_service, scalable=False, area=area)
    except Exception as e:
        logger.info(external_link)
        raise e


def add_load_psse(bus_id: int, id: str, p_set: float, q_set: float, description: str = "", in_service: bool = True, scalable: bool = True, area: int = _i):
    psspy.load_data_6(
        bus_id,
        id,
        [int(in_service), area, _i, _i, int(scalable), _i, _i],
        [p_set, q_set, _f, _f, _f, _f, _f, _f],
        description
    )


def add_shunt_impedance(shunt_impedance: pd.Series, id: int):

    int_shunt_bus = int(shunt_impedance.bus)
    int_bl1_steps = int(shunt_impedance.bl1_steps)
    int_bl2_steps = int(shunt_impedance.bl2_steps)

    psspy.switched_shunt_data_5(
        int_shunt_bus, str(id),
        [int_bl1_steps, int_bl2_steps, _i, _i, _i, _i, _i, _i, 1, _i, _i, _i, _i, _i, _i, _i, _i, _i, _i, _i, _i],
        [float(-20.0), float(20.0), _f, _f, _f, _f, _f, _f, shunt_impedance.v_high, shunt_impedance.v_low, _f, _f], _s
    )


def add_branch(line: pd.Series, id: int):
    # TODO long transmission line model, how does that work in PSS/E?
    # TODO do the PU calculation here instead of in pypsa-eur
    try:
        psspy.branch_data_3(
            int(line.bus0),
            int(line.bus1),
            str(id),
            [_i, _i, _i, _i, _i, _i],
            [line.r_pu, line.x_pu, line.b_pu, _f, _f, _f, _f, line.length, _f, _f, _f, _f],
            [line.s_nom, line.s_nom2, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f],
            "",
        )
    except Exception as e:
        logger.info(line)
        raise e


def remove_branch(remove_lines, remove_buses, remove_transformers):

    for index, row in remove_lines.iterrows():
        psspy.purgbrn(int(row['bus0']), int(row['bus1']), str(row['id']))

    #for index, row in inactive_lines.iterrows():
     #   psspy.branch_chng_3(
      #      int(row['bus0']),
      #      int(row['bus1']),
     #       str(row['id']),
      #      [0, _i, _i, _i, _i, _i],
    #        [_f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f],
     #       [_f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f], _s)

    for index, row in remove_transformers.iterrows():
        psspy.purgbrn(int(row['bus0']), int(row['bus1']), str(row['id']))

    for index, row in remove_buses.iterrows():
        psspy.bsysinit(1)
        psspy.bsyso(1, int(row['bus']))
        psspy.extr(1, 0, [0, 0])


def add_transformer(transformer: pd.Series, id: int):
    try:
        psspy.two_winding_data_6(
            int(transformer.bus0),
            int(transformer.bus1),
            str(id),
            [_i, _i, _i, _i, _i, _i, _i, _i, _i, _i, _i, _i, _i, _i, 2, _i],
            [transformer.r_pu, 0.12, 400.0, _f, _f, _f, _f, _f, _f, _f, _f, transformer.g_pu, transformer.b_pu, _f, _f, _f, _f,
             _f, _f, _f],
            [transformer.s_nom, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f],
            _s,
            _s,
        )

    except Exception as e:
        logger.info(transformer)
        raise e

    psspy.two_winding_data_6(6643, 9411, r"""1""", [1, 6643, 1, 0, 0, 0, 33, 0, 6643, 0, 0, 1, 0, 1, 2, 1],
                             [0.0, 0.12, 400.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.1, 0.9, 1.1,
                              0.9, 0.0, 0.0, 0.0], [2000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0, 0.0, 0.0], "", "")
    psspy.two_winding_data_6(6635, 9410, r"""1""", [1, 6635, 1, 0, 0, 0, 33, 0, 6635, 0, 0, 1, 0, 1, 2, 1],
                             [0.0, 0.12, 400.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.1, 0.9, 1.1,
                              0.9, 0.0, 0.0, 0.0], [2000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0, 0.0, 0.0], "", "")
    psspy.two_winding_data_6(6636, 9410, r"""1""", [1, 6636, 1, 0, 0, 0, 33, 0, 6636, 0, 0, 1, 0, 1, 2, 1],
                             [0.0, 0.12, 400.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.1, 0.9, 1.1,
                              0.9, 0.0, 0.0, 0.0], [2000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0, 0.0, 0.0], "", "")
    #HEMSIL
    psspy.two_winding_data_6(6714, 9015, r"""1""", [1, 6714, 1, 0, 0, 0, 33, 0, 6714, 0, 0, 1, 0, 1, 2, 1],
                             [0.0, 0.12, 400.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.1, 0.9, 1.1,
                              0.9, 0.0, 0.0, 0.0], [2000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0, 0.0, 0.0], "", "")
    #SAUDA
    psspy.two_winding_data_6(6682, 9019, r"""1""", [1, 6682, 1, 0, 0, 0, 33, 0, 6682, 0, 0, 1, 0, 1, 2, 1],
                             [0.0, 0.12, 400.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.1, 0.9, 1.1,
                              0.9, 0.0, 0.0, 0.0], [2000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0, 0.0, 0.0], "", "")
    #FÅBERG
    psspy.two_winding_data_6(6718, 9016, r"""1""", [1, 6718, 1, 0, 0, 0, 33, 0, 6718, 0, 0, 1, 0, 1, 2, 1],
                             [0.0, 0.12, 400.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.1, 0.9, 1.1,
                              0.9, 0.0, 0.0, 0.0], [2000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0, 0.0, 0.0], "", "")
    #OSLO
    psspy.two_winding_data_6(6735, 9017, r"""1""", [1, 6735, 1, 0, 0, 0, 33, 0, 6735, 0, 0, 1, 0, 1, 2, 1],
                             [0.0, 0.12, 400.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.1, 0.9, 1.1,
                              0.9, 0.0, 0.0, 0.0], [2000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0, 0.0, 0.0], "", "")
    #TEGNEBY
    psspy.two_winding_data_6(6739, 9018, r"""1""", [1, 6739, 1, 0, 0, 0, 33, 0, 6739, 0, 0, 1, 0, 1, 2, 1],
                             [0.0, 0.12, 400.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.1, 0.9, 1.1,
                              0.9, 0.0, 0.0, 0.0], [2000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0, 0.0, 0.0], "", "")
    #FJOTLAND
    psspy.two_winding_data_6(6708, 9021, r"""1""", [1, 6708, 1, 0, 0, 0, 33, 0, 6708, 0, 0, 1, 0, 1, 2, 1],
                             [0.0, 0.12, 400.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.1, 0.9, 1.1,
                              0.9, 0.0, 0.0, 0.0], [2000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0, 0.0, 0.0], "", "")
    #LYSE
    psspy.two_winding_data_6(6697, 9020, r"""1""", [1, 6697, 1, 0, 0, 0, 33, 0, 6697, 0, 0, 1, 0, 1, 2, 1],
                             [0.0, 0.12, 400.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.1, 0.9, 1.1,
                              0.9, 0.0, 0.0, 0.0], [2000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0, 0.0, 0.0], "", "")
    #Seitenoikea
    psspy.two_winding_data_6(9101, 6903, r"""1""", [1, 9101, 1, 0, 0, 0, 33, 0, 9101, 0, 0, 1, 0, 1, 2, 1],
                             [0.0, 0.12, 400.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.1, 0.9, 1.1,
                              0.9, 0.0, 0.0, 0.0], [2000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0, 0.0, 0.0], "", "")

    #Nuoja
    psspy.two_winding_data_6(9102, 6923, r"""1""", [1, 9102, 1, 0, 0, 0, 33, 0, 9102, 0, 0, 1, 0, 1, 2, 1],
                             [0.0, 0.12, 400.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.1, 0.9, 1.1,
                              0.9, 0.0, 0.0, 0.0], [2000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0, 0.0, 0.0], "", "")



def add_link(link: pd.Series, bus0: pd.Series, bus1: pd.Series, config: dict):
    # Two possible connections:
    # AC <-> AC: VSC
    # AC <-> DC: Multi-DC --> change topology to combine links so this does not happen
    if bus0.carrier == 'AC' and bus1.carrier == 'AC':
        enabled = ~link.under_construction and link.p_set != 0
        try:
            psspy.vsc_dc_line_data(
                str(link.name),
                [enabled, _i, _i, _i, _i],
                [_f, _f, _f, _f, _f]
            )
            psspy.vsc_dc_converter_data_3(
                str(link.name),  # TODO
                1,
                [int(link.bus0), 2, _i, _i, _i],
                [float(link.p_set), config["build_network"]["scheduled_voltage"], _f, _f, _f, float(link.p_nom), _f, _f, _f, _f, _f]
            )
            psspy.vsc_dc_converter_data_3(
                str(link.name),  # TODO
                2,
                [int(link.bus1), 1, _i, _i, _i],
                [400, config["build_network"]["scheduled_voltage"], _f, _f, _f, float(link.p_nom), _f, _f, _f, _f, _f]
            )
        except Exception as e:
            logger.debug(link)
            raise e


def disconnect_bus(bus: pd.Series):
    try:
        psspy.dscn(int(bus.name))
    except Exception as e:
        logger.info(bus)
        raise e


def solve_newton_raphson() -> bool:
    status_map = {
        0: 'Met convergence tolerance',
        1: 'Iteration limit exceeded',
        2: 'Blown up',
    }

    # increase iteration limit
    psspy.solution_parameters_5([_i, 500, _i, _i, 10, _i, _i, 20, 0],
                                [_f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f, _f,
                                 _f])

    psspy.fnsl([0,0,0,1,1,1,99,0])  # with flat start
    n = psspy.iterat()
    status_val = psspy.solved()
    status = status_map.get(status_val, 'Other - see API documentation of psspy.solved()')
    logger.info(f"{n} iterations used to solve: {status}")
    return status_val == 0


def first_free_index(a: list, b: list) -> int:
    for i in range(1, len(a) + len(b) + 2):
        if i not in a and i not in b:
            return i


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake('build_psse_network', case='high')
    configure_logging(snakemake)

    # when PSS/E encounters an error, don't ignore it
    psspy.throwPsseExceptions = True

    # Initialize PSSE
    psspy.psseinit()

    if snakemake.config["build_network"]["supress_psse_output"]:
        # this prevents PSS/E from spamming a lot of output
        psspy.report_output(6, '', [])
        psspy.progress_output(6, '', [])
        psspy.alert_output(6, '', [])
        psspy.prompt_output(6, '', [])

    # create empty .sav
    psspy.newcase_2([0, 1], snakemake.config["M_base"], 50.0, "Nordics grid", "By DNV (Niek Brekelmans)")

    # Load all data from the sqlite database
    c = sqlite3.connect(snakemake.input.database)
    buses = pd.read_sql("SELECT * FROM buses", con=c).set_index('Bus')

    generators = pd.read_sql("SELECT * FROM generators WHERE p_set > 0", con=c).set_index('Generator')
    storage_units = pd.read_sql("SELECT * FROM storage_units WHERE p_set > 0", con=c).set_index('StorageUnit')
    loads = pd.read_sql("SELECT * FROM loads WHERE p_set > 0", con=c).set_index('Load')
    shunt_impedances = pd.read_sql("SELECT * FROM shunt_impedances", con=c).set_index('ShuntImpedance')
    external_links = pd.read_sql("SELECT * FROM external_links", con=c).set_index('ExternalLink')

    lines = pd.read_sql("SELECT * FROM lines", con=c).set_index('Line')
    transformers = pd.read_sql("SELECT * FROM transformers", con=c).set_index('Transformer')
    links = pd.read_sql("SELECT * FROM links", con=c).set_index('Link')

    # get a sorted list of bidding zones for the area id:  in nordics first, sort the rest alphabetically
    all_bidding_zones = pd.Series(buses["bidding_zone"].unique())
    bidding_zone_in_config = all_bidding_zones.map(lambda country: country in snakemake.config["bidding_zones"])
    all_bidding_zones = pd.concat([
        all_bidding_zones.rename('bidding_zone'),
        bidding_zone_in_config.rename('bidding_zone_in_config')
    ], axis=1)\
        .sort_values(["bidding_zone_in_config", "bidding_zone"], ascending=[False, True])\
        .reset_index()['bidding_zone']

    # get a sorted list of regions for the area id:  in nordics first, sort the rest alphabetically
    all_regions = pd.Series(buses["LnNamn"].unique())
    region_in_config = all_regions.map(lambda country: country in snakemake.config["regions"])
    all_regions = pd.concat([
        all_regions.rename('LnNamn'),
        region_in_config.rename('region_in_config')
    ], axis=1) \
        .sort_values(["region_in_config", "LnNamn"], ascending=[False, True]) \
        .reset_index()['LnNamn']

    # find all buses with external link (only for visualization, 'outside' buses will be disabled)
    buses_with_external_links = set()
    external_links.apply(lambda external_link: buses_with_external_links.add(external_link.bus), axis=1)
    external_links.apply(lambda external_link: buses_with_external_links.add(external_link.bus_outside), axis=1)

    # keep track of buses that are already added in PSS/E
    buses_in_psse = set()

    # add all buses with 1-terminal connections
    logger.info('Adding buses with generators, storage units, loads, shunt impedances and external links')
    for bus_id, bus in buses.iterrows():
        status = bus_included(bus, snakemake.config)
        if not status and bus.name not in buses_with_external_links:
            # don't draw buses in countries which are not in the nordics, except to show external links
            continue

        generators_bus = generators[generators["bus"] == bus.name]
        storage_units_bus = storage_units[storage_units["bus"] == bus.name]
        loads_bus = loads[loads["bus"] == bus.name]
        shunt_impedances_bus = shunt_impedances[shunt_impedances["bus"] == bus.name]
        external_links_bus = external_links[external_links["bus"] == bus.name]

        has_generators = len(generators[generators["bus"] == bus.name].index) > 0 or len(storage_units[storage_units["bus"] == bus.name].index) > 0

        try:
            region_index = all_regions[all_regions == bus.LnNamn].index[0] + 1
        except IndexError:
            # region name not found, set region index to 0 or handle the error
            region_index = 0
        bidding_zone_index = all_bidding_zones[all_bidding_zones == bus.bidding_zone].index[0] + 1
        add_bus(bus, bidding_zone_index, region_index, has_generators, snakemake.config)
        buses_in_psse.add(bus.name)

        # add_bus(bus, all_bidding_zones[all_bidding_zones == bus.bidding_zone].index[0] + 1, has_generators, snakemake.config)
        # buses_in_psse.add(bus.name)

        index = 1
        for _, generator in generators_bus.iterrows():
            if add_machine(generator, bus, index, status, snakemake.config):
                index += 1

        for _, storage_unit in storage_units_bus.iterrows():
            if add_machine(storage_unit, bus, index, status, snakemake.config):
                index += 1

        index = 1
        for _, load in loads_bus.iterrows():
            add_load(load, bus, index, status)
            index += 1

        index = 1
        for _, shunt_impedance in shunt_impedances_bus.iterrows():
            add_shunt_impedance(shunt_impedance, index)
            index += 1

        index = 1
        for _, external_link in external_links_bus.iterrows():
            bus_outside = buses.loc[external_link.bus_outside]
            area = all_bidding_zones[all_bidding_zones == bus_outside.bidding_zone].index[0] + 1
            add_external_link(external_link, bus, index, status, area)
            index += 1

    logger.info('Adding lines')
    # indexes are meant to find a 'free' id (for both buses) when connecting 2 buses
    indexes = {k: [] for k in buses.index.unique()}
    for line_id, line in lines.iterrows():
        if line.bus0 not in indexes.keys() or line.bus1 not in indexes.keys() \
                or line.bus0 not in buses_in_psse or line.bus1 not in buses_in_psse:
            continue

        index = first_free_index(indexes[line.bus0], indexes[line.bus1])
        add_branch(line, index)
        indexes[line.bus0].append(index)
        indexes[line.bus1].append(index)

    logger.info('Adding transformers')
    # indexes are meant to find a 'free' id (for both buses) when connecting 2 buses
    indexes = {k: [] for k in buses.index.unique()}
    for transformer_id, transformer in transformers.iterrows():
        if transformer.bus0 not in indexes.keys() or transformer.bus1 not in indexes.keys() \
                or transformer.bus0 not in buses_in_psse or transformer.bus1 not in buses_in_psse:
            continue

        index = first_free_index(indexes[transformer.bus0], indexes[transformer.bus1])
        add_transformer(transformer, index)
        indexes[transformer.bus0].append(index)
        indexes[transformer.bus1].append(index)

    logger.info('Removing lines and buses')
    remove_lines = pd.read_excel(snakemake.input.remove_lines)
    #inactive_lines = pd.read_excel(snakemake.input.inactive_lines)
    remove_buses = pd.read_excel(snakemake.input.remove_buses)
    remove_transfomers = pd.read_excel(snakemake.input.remove_transformers)
    remove_branch(remove_lines, remove_buses, remove_transfomers)


    logger.info('Adding links')
    # indexes are meant to find a 'free' id (for both buses) when connecting 2 buses
    indexes = {k: [] for k in buses.index.unique()}
    for link_id, link in links.iterrows():
        if link.bus0 not in indexes.keys() or link.bus1 not in indexes.keys() \
                or link.bus0 not in buses_in_psse or link.bus1 not in buses_in_psse:
            continue

        bus0 = buses.loc[link.bus0]
        bus1 = buses.loc[link.bus1]
        add_link(link, bus0, bus1, snakemake.config)

    # add countries as area; TODO (?): trade zones within countries -> manual work
    logger.info('Adding areas')
    for i, country in enumerate(all_bidding_zones, start=1):
        add_area(i, country)

    logger.info('Adding regions')
    for i, country in enumerate(all_regions, start=1):
        add_region(i, country)

    logger.info('Disabling excluded buses')
    for _, bus in buses.iterrows():
        if not bus_included(bus, snakemake.config) and bus.name in buses_in_psse:
            disconnect_bus(bus)

    # run load flow
    logger.info('Solving using Newton-Raphson (flat start)')
    solve_newton_raphson()

    # save case and diagram to file
    logger.info('Saving ' + snakemake.output.sav)
    error = psspy.save(snakemake.output.sav)

    # close PSS/E
    logger.info('Closing PSS/E')
    psspy.deltmpfiles()

    psspy.pssehalt_2()
