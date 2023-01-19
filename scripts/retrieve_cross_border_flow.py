import logging
from operator import xor
from typing import List

import entsoe.mappings
import numpy as np
import pandas as pd
from entsoe.exceptions import NoMatchingDataError, InvalidBusinessParameterError

from _helpers import configure_logging, get_start_and_end_of_year, all_combinations, get_entsoe_client

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake('retrieve_cross_border_flow')
    configure_logging(snakemake)

    client = get_entsoe_client(snakemake.config)

    snapshots = pd.read_csv(snakemake.input.snapshots, index_col=0)
    snapshot = snapshots.loc[snakemake.wildcards.case, 'snapshot']

    start = pd.Timestamp(snapshot, tz='UTC')
    end = start + pd.DateOffset(hours=1)

    flows = pd.DataFrame()
    flows.index.name = 'bidding_zone_from'

    # entsoe package knows neighbouring bidding zones
    included_cross_border_flows = { zone: entsoe.mappings.NEIGHBOURS[zone] for zone in snakemake.config["bidding_zones"] }

    for bidding_zone_from in included_cross_border_flows:
        for bidding_zone_to in included_cross_border_flows[bidding_zone_from]:
            index = bidding_zone_from + "-" + bidding_zone_to
            logger.info(f"Retrieving cross border flow {index}")
            try:
                flow_from_nordics = client.query_crossborder_flows(
                    bidding_zone_from,
                    bidding_zone_to,
                    start=start,
                    end=end
                ).loc[start]

                flow_to_nordics = client.query_crossborder_flows(
                    bidding_zone_to,
                    bidding_zone_from,
                    start=start,
                    end=end
                ).loc[start]

                flow = flow_from_nordics - flow_to_nordics
            except InvalidBusinessParameterError:  # link really does not exist
                flow = np.nan
                pass
            except NoMatchingDataError:  # link has no data for this timestamp
                flow = 0
                pass

            flows.loc[bidding_zone_from, bidding_zone_to] = flow
    flows.to_csv(snakemake.output.cross_border_flows)
