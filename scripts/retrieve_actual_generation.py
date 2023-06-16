import logging

import pandas as pd
from entsoe import EntsoePandasClient

from _helpers import configure_logging, get_start_and_end_of_year, get_empty_generation_df, get_entsoe_client, \
    all_areas_entsoe

logger = logging.getLogger(__name__)


def get_data_from_entsoe(config: dict):
    client = get_entsoe_client(config)

    generation = get_empty_generation_df(config["generation"])

    snapshot = pd.Timestamp(config["snapshot"], tz='UTC')

    # gather data from ENTSO-E TP for each area
    for area in all_areas_entsoe(config):
        df = client.query_generation(
            area,
            start=snapshot,
            end=snapshot + pd.DateOffset(hours=1),
        )

        df = df.loc[snapshot].reset_index()
        df = df.set_axis(['production_type', 'generation'], axis=1)

        df["production_type_mapped"] = df["production_type"].map(
            lambda production_type: config["generation"]["entsoe_map"].get(production_type, 'other')
        )
        generation[area] = df.groupby('production_type_mapped')['generation'].sum()
    return generation


# def apply_manual_correction_for_sweden(generation: pd.DataFrame, sweden: pd.DataFrame):
#     # rescale with total entsoe data so SE is equal to the sum of all swedish bidding zones
#     sweden_totals = sweden.sum(axis=1)
#     sweden_totals_entsoe = generation['SE']
#
#     scale_factors_series = (sweden_totals_entsoe / sweden_totals).fillna(0)
#
#     #scale_factors_series = (sweden_totals/ sweden_totals).fillna(0)
#
#     # same scale factors for each zone
#     scale_factors = pd.DataFrame()
#     scale_factors['SE_1'] = scale_factors_series
#     scale_factors['SE_2'] = scale_factors_series
#     scale_factors['SE_3'] = scale_factors_series
#     scale_factors['SE_4'] = scale_factors_series
#
#     # scale
#     sweden = sweden * scale_factors
#
#     # overwrite missing entsoe data with scaled SvK data
#     generation['SE_1'] = sweden['SE_1']
#     generation['SE_2'] = sweden['SE_2']
#     generation['SE_3'] = sweden['SE_3']
#     generation['SE_4'] = sweden['SE_4']
#
#     # there are no 'other' generators in sweden. Move that power to hydro (largest type), so the power is not 'lost'
#     generation.loc['hydro', 'SE_2'] += generation.loc['other', 'SE_2']
#     generation.loc['other', 'SE_2'] = 0
#
#     return generation


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake('retrieve_actual_generation', case="high")
    configure_logging(snakemake)

    generation = get_data_from_entsoe(snakemake.config)

    # read and apply manual SvK data, from 'Mimer' portal:
    # https://mimer.svk.se/ProductionConsumption/ProductionIndex
    #sweden = pd.read_excel(snakemake.input.generation_sweden, index_col=0)
    #generation = apply_manual_correction_for_sweden(generation, sweden)

    # save the ENTSO-E generation data (with the manual corrections applied)
    generation.to_csv(snakemake.output.actual_generation)
