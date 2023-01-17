import logging

import pandas as pd
from entsoe import EntsoePandasClient

from _helpers import configure_logging, get_start_and_end_of_year, get_empty_generation_df

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake('retrieve_actual_generation', case="high")
    configure_logging(snakemake)

    client = EntsoePandasClient(api_key=snakemake.config["entsoe"]["security_token"])

    generation = get_empty_generation_df(snakemake.config["generation"])

    snapshots = pd.read_csv(snakemake.input.snapshots, index_col=0)
    snapshot = snapshots.loc[snakemake.wildcards.case, 'snapshot']
    snapshot = pd.Timestamp(snapshot, tz="UTC")

    for country in snakemake.config["countries"]:
        country_entsoe = snakemake.config["entsoe"]["country_map"].get(country, country)  # for mapping DK to DK2
        df = client.query_generation(
            country_entsoe,
            start=snapshot,
            end=snapshot + pd.DateOffset(hours=1),
        )

        df = df.loc[snapshot].reset_index()
        df = df.set_axis(['production_type', 'generation'], axis=1)

        df["production_type_mapped"] = df["production_type"].map(
            lambda production_type: snakemake.config["generation"]["entsoe_map"].get(production_type, 'other')
        )
        generation[country] = df.groupby('production_type_mapped')['generation'].sum()

    generation.sort_index(inplace=True)
    generation.to_csv(snakemake.output.actual_generation)




