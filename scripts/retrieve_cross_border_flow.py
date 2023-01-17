import logging
from typing import List

import numpy as np
import pandas as pd
from entsoe.exceptions import NoMatchingDataError

from _helpers import configure_logging, get_start_and_end_of_year
from entsoe import EntsoePandasClient

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake('retrieve_cross_border_flow')
    configure_logging(snakemake)

    client = EntsoePandasClient(api_key=snakemake.config["entsoe"]["security_token"])

    year = snakemake.config["year"]
    (start, end) = get_start_and_end_of_year(year)
    date_range = pd.date_range(start=start, end=end, freq="H", tz=None)

    flows = pd.DataFrame(index=date_range)
    flows.index.rename('datetime')
    external_links = snakemake.config["external_links"]
    for country_in_nordics in external_links:
        for country_outside_nordics in external_links[country_in_nordics]:
            index = country_in_nordics + "-" + country_outside_nordics
            logger.info(f"Retrieving cross border flow {index} ({year})")
            try:
                flow_from_nordics = client.query_crossborder_flows(country_in_nordics, country_outside_nordics, start=start, end=end)
                flow_to_nordics = client.query_crossborder_flows(country_outside_nordics, country_in_nordics, start=start, end=end)

                flows[index] = flow_from_nordics - flow_to_nordics
            except NoMatchingDataError:
                flows[index] = np.nan

    flows.to_csv(snakemake.output.cross_border_flows)
