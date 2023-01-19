import logging
import sqlite3

import pandas as pd

from _helpers import configure_logging

import psse35
import psspy
from psspy import _i, _f, _s
import pssgrpg

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("validate_psse", case="high")
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

    _, numbers = psspy.aareaint(string='NUMBER')
    _, names = psspy.aareachar(string='AREANAME')

    numbers = numbers[0]
    names = names[0]

    areas = dict(map(lambda i, j: (i, j), numbers, names))
    areas = {number: name.strip() for number, name in areas.items() if name.strip() in snakemake.config["bidding_zones"]}

    inter_area_flow = pd.DataFrame()
    for area_from in areas:
        for area_to in areas:
            p, _ = (0, 0)
            if area_from != area_to:
                try:
                    p, _ = pssgrpg.area_interchange_ij_v(area_from, area_to)
                except AttributeError:
                    pass

            inter_area_flow.loc[areas[area_from], areas[area_to]] = round(p, 0)

    print(inter_area_flow)
