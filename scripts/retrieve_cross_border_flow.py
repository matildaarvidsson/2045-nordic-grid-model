import logging
from operator import xor
from typing import List

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

    # year = snakemake.config["year"]
    # (start, end) = get_start_and_end_of_year(year)
    # date_range = pd.date_range(start=start, end=end, freq="H", tz=None, name='datetime')

    snapshots = pd.read_csv(snakemake.input.snapshots, index_col=0)
    snapshot = snapshots.loc[snakemake.wildcards.case, 'snapshot']

    start = pd.Timestamp(snapshot, tz='UTC')
    end = start + pd.DateOffset(hours=1)

    flows = pd.DataFrame()
    flows.index.name = 'country_from'
    included_cross_border_flows = snakemake.config["cross_border_flows"]
    for country_from in included_cross_border_flows:
        for country_to in included_cross_border_flows[country_from]:
            index = country_from + "-" + country_to
            logger.info(f"Retrieving cross border flow {index}")

            # only possible to query country <=> country or zone <=> zone
            if xor('_' in country_from, '_' in country_to):
                bidding_zone_combinations = all_combinations(
                    snakemake.config["bidding_zones_in_country"].get(country_from, [country_from]),
                    snakemake.config["bidding_zones_in_country"].get(country_to, [country_to]),
                )
            else:
                bidding_zone_combinations = [(country_from, country_to)]

            flow = 0
            exists = False
            for bidding_zone_from, bidding_zone_to in bidding_zone_combinations:
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

                    exists = True

                    flow += flow_from_nordics - flow_to_nordics
                except InvalidBusinessParameterError:  # link really does not exist
                    pass
                except NoMatchingDataError:  # link has no data for this timestamp
                    exists = True
                    pass

            flows.loc[country_from, country_to] = flow if exists else np.nan
    flows.to_csv(snakemake.output.cross_border_flows)
