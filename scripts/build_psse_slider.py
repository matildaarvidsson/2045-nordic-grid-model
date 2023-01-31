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


def draw_bus(bus: pd.Series):
    psspy.updatebuslocdiagfile()
    try:
        psspy.growdiagram_2(
            0,
            1,
            ["BU   {0}".format(bus.name)],
            _f,
            _f,
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        )
    except Exception as e:
        logger.info(bus)
        raise e


def create_loc_file(all_buses: pd.DataFrame, filename: str):
    lines = [
        "GEOPHYSICAL LATLON\n",
        "/BusID, Y location, X location\n",
    ]

    all_buses.apply(lambda bus: lines.append("{id}, {lat}, {lon}, ,\n".format(id=bus.name, lat=bus.y, lon=bus.x)), axis=1)

    logger.info('Saving ' + filename)
    f = open(filename, "w")
    f.writelines(lines)
    f.close()

    psspy.openbuslocfile(filename)
    psspy.updatebuslocdiagfile()


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

    # Open case
    psspy.case(snakemake.input.sav)

    # create empty .sld
    psspy.newdiagfile()

    # from recorded action; sets the colors of voltage levels to distinguish between 400 kV, 300 kV and 220 kV
    psspy.setdiagresvrcs_2(2, 1, [1, 1, 1, 1, 1, 1, 1, 1], [-1, -1, -1, -1, -1, -1, -1, -1],
                           [0.0, 13.8, 18.0, 20.0, 220.0, 300.0, 380.0], [0, 255, 0, 0, 255, 64, 0, 0],
                           [0, 0, 0, 255, 0, 128, 0, 0], [0, 0, 255, 0, 128, 128, 219, 255], 1, 1, 63736, [64, 0, 64],
                           [0, 0, 0], [255, 0, 255], [0, 255, 0], [255, 0, 0], _i, _i, _i, _i, [_i, 0, 0])

    # Load all data from the sqlite database
    c = sqlite3.connect(snakemake.input.database)
    buses = pd.read_sql("SELECT * FROM buses", con=c).set_index('Bus')

    # a .loc file contains the coordinates of all buses to draw the SLD diagram geographically correct
    create_loc_file(buses, snakemake.output.loc)

    # draw the sld
    logger.info('Drawing buses (can take some time)')
    buses.apply(draw_bus, axis=1)

    logger.info('Saving ' + snakemake.output.sld)
    psspy.savediagfile(snakemake.output.sld)

    # close PSS/E
    logger.info('Closing PSS/E')
    psspy.deltmpfiles()

    psspy.pssehalt_2()
