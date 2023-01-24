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

    snapshot = snakemake.config["snapshot"]
    start = pd.Timestamp(snapshot, tz='UTC')
    end = start + pd.DateOffset(hours=1)

    flows = pd.DataFrame()
    flows.index.name = 'bidding_zone_from'

    # entsoe package knows neighbouring bidding zones
    included_cross_border_flows = {
        zone: entsoe.mappings.NEIGHBOURS[zone] for zone in snakemake.config["bidding_zones"]
    }
    included_cross_border_flows['SE_4'].append('SE_3')  # missing in the package list

    for bidding_zone_from in included_cross_border_flows:
        for bidding_zone_to in included_cross_border_flows[bidding_zone_from]:
            if bidding_zone_from == bidding_zone_to:
                continue
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
                logger.warning(f"InvalidBusinessParameterError in {index}")
                flow = np.nan
                pass
            except NoMatchingDataError:  # link has no data for this timestamp
                logger.warning(f"NoMatchingDataError in {index}")
                flow = 0
                pass

            flows.loc[bidding_zone_from, bidding_zone_to] = flow

    # DE_AT_LU for date < october 2018
    # DE_LU    for date > october 2018
    # so manual fix required
    flows['DE'] = flows['DE_AT_LU'] + flows['DE_LU']
    flows.drop(columns=['DE_AT_LU', 'DE_LU'], inplace=True)
    flows.to_csv(snakemake.output.cross_border_flows)
