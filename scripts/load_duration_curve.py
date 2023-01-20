import logging

import pandas as pd
from matplotlib import pyplot as plt
from pypsa import Network

from _helpers import configure_logging

logger = logging.getLogger(__name__)


def plot_load_duration_curve(load: pd.Series, load_duration_curve_file: str):
    p = load.plot(x="percentage", y="total_load")

    p.set_title("Load-Duration Curve", fontsize=20)
    p.set_xlabel("Time (%)", fontsize=15)
    p.set_ylabel("Load (MW)", fontsize=15)

    plt.axvline(x=10, color='red')
    plt.axvline(x=90, color='red')
    plt.savefig(load_duration_curve_file)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake('pypsa_snapshots')
    configure_logging(snakemake)

    n = Network(snakemake.input.network)

    # load duration curve
    load = n.loads_t.p_set
    load['total_load'] = load.sum(axis=1)
    load = load.sort_values(by=['total_load'], ascending=False)
    load['interval'] = 1
    load['duration'] = load['interval'].cumsum()
    load['percentage'] = load['duration'] * 100 / len(load)

    snapshot_high = load["percentage"].ge(10).idxmax()
    snapshot_low = load["percentage"].ge(90).idxmax()

    logger.info('Snapshot high: '+snapshot_high.strftime("%Y-%m-%d, %H:%M"))
    logger.info('Snapshot low: '+snapshot_low.strftime("%Y-%m-%d, %H:%M"))

    plot_load_duration_curve(load, snakemake.output.load_duration_curve)


