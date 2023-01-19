import logging
import sqlite3

import numpy as np
import pandas as pd

from _helpers import configure_logging

import psse35
import psspy
from psspy import _i, _f, _s
import pssgrpg

logger = logging.getLogger(__name__)


def psse_get_areas() -> dict[int, str]:
    _, numbers = psspy.aareaint(string='NUMBER')
    _, names = psspy.aareachar(string='AREANAME')

    numbers = numbers[0]
    names = names[0]

    areas = dict(map(lambda i, j: (i, j), numbers, names))
    return {number: name.strip() for number, name in areas.items() if name.strip() in snakemake.config["bidding_zones"]}


def psse_get_inter_area_flow(areas: dict[int, str]) -> pd.DataFrame:
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

    inter_area_flow.index.names = ['country_from']
    return inter_area_flow


def compare_cross_border_flows(entsoe_flows: pd.DataFrame, psse_flows: pd.DataFrame):
    with pd.ExcelWriter(snakemake.output.cross_border_flow, engine='xlsxwriter') as writer:
        for bidding_zone in snakemake.config["bidding_zones"]:
            df = pd.DataFrame()
            df['model'] = psse_flows.loc[bidding_zone]
            df['entsoe'] = entsoe_flows.loc[bidding_zone]

            df.fillna(0, inplace=True)

            df.loc['total'] = df.sum(axis=0)

            df['difference'] = df['model'] - df['entsoe']
            df['percentage'] = df['model'] / df['entsoe']

            df.to_excel(writer, sheet_name=bidding_zone)

        # format sheets
        number_format = writer.book.add_format({'num_format': '0'})
        percent_format = writer.book.add_format({'num_format': '0%'})
        for bidding_zone in writer.sheets:
            worksheet = writer.sheets[bidding_zone]
            worksheet.set_column(1, 3, None, number_format)
            worksheet.set_column(4, 4, None, percent_format)


def fill_reverse(entsoe_flows: pd.DataFrame):
    for rowIndex, row in entsoe_flows.iterrows():
        for columnIndex, value in row.items():
            if ~np.isnan(value):
                entsoe_flows.loc[columnIndex, rowIndex] = 0 if value == 0 else -value  # set reverse flow
    return entsoe_flows


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

    #
    # Compare cross-border flows
    #
    areas = psse_get_areas()

    inter_area_flow = psse_get_inter_area_flow(areas)
    entsoe_cross_border_flows = pd.read_csv(snakemake.input.cross_border_flows, index_col=0)
    entsoe_cross_border_flows = fill_reverse(entsoe_cross_border_flows)

    compare_cross_border_flows(
        entsoe_cross_border_flows,
        inter_area_flow,
    )
