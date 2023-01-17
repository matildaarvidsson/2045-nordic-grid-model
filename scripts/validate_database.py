import logging
import sqlite3
from sqlite3 import Connection

import pandas as pd

from _helpers import configure_logging, bus_included, get_empty_generation_df

logger = logging.getLogger(__name__)


def get_generation_totals(column: str, c: Connection) -> pd.DataFrame:
    buses = pd.read_sql("SELECT * FROM buses", con=c).set_index('Bus')
    generators = pd.read_sql("SELECT * FROM generators", con=c).set_index('Generator')
    storage_units = pd.read_sql("SELECT * FROM storage_units", con=c).set_index('StorageUnit')
    generation_units = pd.concat([generators, storage_units])

    generation_units["country"] = generation_units["bus"].map(lambda bus: buses.loc[bus, "country"])
    generation_units["in_synchronous_network"] = generation_units["bus"].map(
        lambda bus: buses.loc[bus, "in_synchronous_network"]
    )

    # filter out not in nordics sync network
    generation_units = generation_units[generation_units["in_synchronous_network"] == True]

    return generation_units.groupby(['country', 'type'])[column].sum().unstack(level=0)


def validate_generation_capacity(c: Connection) -> None:
    empty_df = get_empty_generation_df(snakemake.config["generation"])

    totals = empty_df.copy()
    totals["model"] = 0
    totals["entsoe"] = 0

    entsoe_capacity = pd.read_csv(snakemake.input.entsoe_capacities).set_index('type')
    model_capacity = get_generation_totals('p_nom', c)

    with pd.ExcelWriter(snakemake.output.installed_generation_capacity, engine='xlsxwriter') as writer:
        for country in snakemake.config["countries"]:
            df = empty_df.copy()
            df['model'] = model_capacity[country]
            df['entsoe'] = entsoe_capacity[country]
            df.fillna(0, inplace=True)

            df.loc['total'] = df.sum(axis=0)

            df['difference'] = df['model'] - df['entsoe']
            df['percentage'] = df['model'] / df['entsoe']

            totals["model"] = totals["model"] + df["model"]
            totals["entsoe"] = totals["entsoe"] + df["entsoe"]

            df.to_excel(writer, sheet_name=country)

        totals.loc['total'] = totals.sum(axis=0)
        totals['difference'] = totals['model'] - totals['entsoe']
        totals['percentage'] = totals['model'] / totals['entsoe']
        totals.to_excel(writer, sheet_name='Total')

        # format sheets
        number_format = writer.book.add_format({'num_format': '0'})
        percent_format = writer.book.add_format({'num_format': '0%'})
        for country in writer.sheets:
            worksheet = writer.sheets[country]
            worksheet.set_column(1, 3, None, number_format)
            worksheet.set_column(4, 4, None, percent_format)


def validate_actual_generation(c: Connection) -> None:
    empty_df = get_empty_generation_df(snakemake.config["generation"])

    totals = empty_df.copy()
    totals["model"] = 0
    totals["entsoe"] = 0

    entsoe_generation = pd.read_csv(snakemake.input.entsoe_generation).set_index('type')
    model_capacity = get_generation_totals('p_set', c)

    with pd.ExcelWriter(snakemake.output.actual_generation, engine='xlsxwriter') as writer:
        for country in snakemake.config["countries"]:
            df = empty_df.copy()
            df['model'] = model_capacity[country]
            df['entsoe'] = entsoe_generation[country]
            df.fillna(0, inplace=True)

            df.loc['total'] = df.sum(axis=0)

            df['difference'] = df['model'] - df['entsoe']
            df['percentage'] = df['model'] / df['entsoe']

            df.to_excel(writer, sheet_name=country)

            totals["model"] = totals["model"] + df["model"]
            totals["entsoe"] = totals["entsoe"] + df["entsoe"]

        totals.loc['total'] = totals.sum(axis=0)
        totals['difference'] = totals['model'] - totals['entsoe']
        totals['percentage'] = totals['model'] / totals['entsoe']
        totals.to_excel(writer, sheet_name='Total')

        # format sheets
        number_format = writer.book.add_format({'num_format': '0'})
        percent_format = writer.book.add_format({'num_format': '0%'})
        for country in writer.sheets:
            worksheet = writer.sheets[country]
            worksheet.set_column(1, 3, None, number_format)
            worksheet.set_column(4, 4, None, percent_format)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("validate_database", case="high")
    configure_logging(snakemake)

    c = sqlite3.connect(snakemake.input.database)
    validate_generation_capacity(c)
    validate_actual_generation(c)
