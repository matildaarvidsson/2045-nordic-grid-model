import logging
import sqlite3
from sqlite3 import Connection

import pandas as pd

from _helpers import configure_logging, bus_included, get_empty_generation_df

logger = logging.getLogger(__name__)


def get_generation_totals(column: str, group_by: str, c: Connection) -> pd.DataFrame:
    buses = pd.read_sql("SELECT * FROM buses", con=c).set_index('Bus')
    generators = pd.read_sql("SELECT * FROM generators", con=c).set_index('Generator')
    storage_units = pd.read_sql("SELECT * FROM storage_units", con=c).set_index('StorageUnit')
    generation_units = pd.concat([generators, storage_units])

    generation_units[group_by] = generation_units["bus"].map(lambda bus: buses.loc[bus, group_by])
    generation_units["in_synchronous_network"] = generation_units["bus"].map(
        lambda bus: buses.loc[bus, "in_synchronous_network"]
    )

    # filter out not in nordics sync network
    generation_units = generation_units[generation_units["in_synchronous_network"] == True]

    return generation_units.groupby([group_by, 'type'])[column].sum().unstack(level=0)


def validate_generation_capacity(c: Connection, country: bool) -> None:
    groupby = 'country' if country else 'bidding_zone'
    output = snakemake.output.generation_capacity_country if country else snakemake.output.generation_capacity_bidding_zone
    areas = snakemake.config["countries"] if country else snakemake.config["bidding_zones"]

    empty_df = get_empty_generation_df(snakemake.config["generation"])

    totals = empty_df.copy()
    totals["model"] = 0
    totals["entsoe"] = 0

    entsoe_capacity = pd.read_csv(snakemake.input.entsoe_capacities).set_index('type')
    model_capacity = get_generation_totals('p_nom', groupby, c)

    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        for area in areas:
            df = empty_df.copy()
            df['model'] = model_capacity[area]
            df['entsoe'] = entsoe_capacity[area]
            df.fillna(0, inplace=True)

            df.loc['total'] = df.sum(axis=0)

            df['difference'] = df['model'] - df['entsoe']
            df['percentage'] = df['model'] / df['entsoe']

            totals["model"] = totals["model"] + df["model"]
            totals["entsoe"] = totals["entsoe"] + df["entsoe"]

            df.to_excel(writer, sheet_name=area)

        totals.loc['total'] = totals.sum(axis=0)
        totals['difference'] = totals['model'] - totals['entsoe']
        totals['percentage'] = totals['model'] / totals['entsoe']
        totals.to_excel(writer, sheet_name='Total')

        # format sheets
        number_format = writer.book.add_format({'num_format': '0'})
        percent_format = writer.book.add_format({'num_format': '0%'})
        for area in writer.sheets:
            worksheet = writer.sheets[area]
            worksheet.set_column(1, 3, None, number_format)
            worksheet.set_column(4, 4, None, percent_format)


def validate_balance():
    cross_border_flow = pd.read_csv(snakemake.input.entsoe_cross_border_flow, index_col=0)
    generation = pd.read_csv(snakemake.input.entsoe_generation, index_col=0)
    load = pd.read_csv(snakemake.input.entsoe_load, index_col=0)

    balance = pd.DataFrame()

    balance['generation'] = generation.sum() * -1
    balance['load'] = load['load']
    balance['cross_border_flow'] = cross_border_flow.sum(axis=1)

    balance['balance'] = balance.sum(axis=1)

    balance.to_excel(snakemake.output.balance)


def validate_actual_generation(c: Connection, country: bool) -> None:
    groupby = 'country' if country else 'bidding_zone'
    output = snakemake.output.actual_generation_country if country else snakemake.output.actual_generation_bidding_zone
    areas = snakemake.config["countries"] if country else snakemake.config["bidding_zones"]

    empty_df = get_empty_generation_df(snakemake.config["generation"])

    totals = empty_df.copy()
    totals["model"] = 0
    totals["entsoe"] = 0

    entsoe_generation = pd.read_csv(snakemake.input.entsoe_generation).set_index('type')
    model_generation = get_generation_totals('p_set', groupby, c)

    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        for area in areas:
            df = empty_df.copy()
            df['model'] = model_generation[area]
            df['entsoe'] = entsoe_generation[area]
            df.fillna(0, inplace=True)

            df.loc['total'] = df.sum(axis=0)

            df['difference'] = df['model'] - df['entsoe']
            df['percentage'] = df['model'] / df['entsoe']

            df.to_excel(writer, sheet_name=area)

            totals["model"] = totals["model"] + df["model"]
            totals["entsoe"] = totals["entsoe"] + df["entsoe"]

        totals.loc['total'] = totals.sum(axis=0)
        totals['difference'] = totals['model'] - totals['entsoe']
        totals['percentage'] = totals['model'] / totals['entsoe']
        totals.to_excel(writer, sheet_name='Total')

        # format sheets
        number_format = writer.book.add_format({'num_format': '0'})
        percent_format = writer.book.add_format({'num_format': '0%'})
        for area in writer.sheets:
            worksheet = writer.sheets[area]
            worksheet.set_column(1, 3, None, number_format)
            worksheet.set_column(4, 4, None, percent_format)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("validate_database", case="high")
    configure_logging(snakemake)

    c = sqlite3.connect(snakemake.input.database)

    validate_generation_capacity(c, True)
    validate_generation_capacity(c, False)

    validate_actual_generation(c, True)
    validate_actual_generation(c, False)

    validate_balance()
