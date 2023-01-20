import logging

import pandas as pd
from entsoe import EntsoePandasClient

from _helpers import configure_logging, get_start_and_end_of_year, get_empty_generation_df, get_entsoe_client, \
    all_areas_entsoe

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake('retrieve_actual_generation', case="high")
    configure_logging(snakemake)

    client = get_entsoe_client(snakemake.config)

    load = pd.DataFrame(dtype='float64')
    load.index.name = 'bidding_zone'

    snapshots = pd.read_csv(snakemake.input.snapshots, index_col=0)
    snapshot = snapshots.loc[snakemake.wildcards.case, 'snapshot']
    snapshot = pd.Timestamp(snapshot, tz="UTC")

    for area in all_areas_entsoe(snakemake.config):
        df = client.query_load(
            area,
            start=snapshot,
            end=snapshot + pd.DateOffset(hours=1),
        )

        load.loc[area, 'load'] = df.loc[snapshot, 'Actual Load']

    load.sort_index(inplace=True)
    load.to_csv(snakemake.output.load)
