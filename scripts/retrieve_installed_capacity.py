import logging

from entsoe import EntsoePandasClient

from _helpers import configure_logging, get_start_and_end_of_year, get_empty_generation_df, get_entsoe_client

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake('retrieve_installed_capacity')
    configure_logging(snakemake)

    client = get_entsoe_client(snakemake.config)
    capacity = get_empty_generation_df(snakemake.config["generation"])
    (start, end) = get_start_and_end_of_year(snakemake.config["year"])  # to get installed capacity at END of the year
    for bidding_zone in snakemake.config["bidding_zones"]:
        df = client.query_installed_generation_capacity(
            bidding_zone,
            start=start,
            end=end,
        ).iloc[0].reset_index()

        df = df.set_axis(['production_type', 'capacity'], axis=1)

        df["production_type_mapped"] = df["production_type"].map(
            lambda production_type: snakemake.config["generation"]["entsoe_map"].get(production_type, 'other')
        )
        capacity[bidding_zone] = df.groupby('production_type_mapped')['capacity'].sum()

    capacity.sort_index().to_csv(snakemake.output.installed_capacity)
