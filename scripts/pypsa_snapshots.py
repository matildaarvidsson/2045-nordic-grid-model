import logging

import pandas as pd
from matplotlib import pyplot as plt
from pypsa import Network

from _helpers import configure_logging

logger = logging.getLogger(__name__)


def plot_load_duration_curve(load: pd.Series, load_options: pd.DataFrame, load_duration_curve_file: str):
    p = load.plot(x="percentage", y="total_load")

    p.set_title("Load-Duration Curve", fontsize=20)
    p.set_xlabel("Time (%)", fontsize=15)
    p.set_ylabel("Load (MW)", fontsize=15)

    load_options.apply(lambda row: plt.axvline(x=row["percentage"], color='red'), axis=1)
    plt.savefig(load_duration_curve_file)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake('pypsa_snapshots')
    configure_logging(snakemake)

    n = Network(snakemake.input.network)

    load_options = pd.DataFrame(snakemake.config["snapshots"])
    load_options = load_options.transpose()
    load_options["percentage"] = pd.to_numeric(load_options["percentage"])

    # load duration curve
    load = n.loads_t.p_set
    load['total_load'] = load.sum(axis=1)
    load = load.sort_values(by=['total_load'], ascending=False)
    load['interval'] = 1
    load['duration'] = load['interval'].cumsum()
    load['percentage'] = load['duration'] * 100 / len(load)

    load_options["snapshot"] = load_options["snapshot"].fillna(load_options["percentage"].map(lambda x: load["percentage"].ge(x).idxmax()))
    load_options["snapshot"].to_csv(snakemake.output.snapshots)

    plot_load_duration_curve(load, load_options, snakemake.output.load_duration_curve)


