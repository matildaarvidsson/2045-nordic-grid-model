import itertools
from pathlib import Path
from typing import Tuple

import pandas as pd
import logging

from entsoe import EntsoePandasClient

logger = logging.getLogger(__name__)


def configure_logging(snakemake, skip_handlers=False):
    """
    Configure the basic behaviour for the logging module.

    Note: Must only be called once from the __main__ section of a script.

    The setup includes printing log messages to STDERR and to a log file defined
    by either (in priority order): snakemake.log.python, snakemake.log[0] or "logs/{rulename}.log".
    Additional keywords from logging.basicConfig are accepted via the snakemake configuration
    file under snakemake.config.logging.

    Parameters
    ----------
    snakemake : snakemake object
        Your snakemake object containing a snakemake.config and snakemake.log.
    skip_handlers : True | False (default)
        Do (not) skip the default handlers created for redirecting output to STDERR and file.
    """

    kwargs = snakemake.config.get("logging", dict()).copy()
    kwargs.setdefault("level", "INFO")

    if skip_handlers is False:
        fallback_path = Path(__file__).parent.joinpath(
            "..", "logs", f"{snakemake.rule}.log"
        )
        logfile = snakemake.log.get(
            "python", snakemake.log[0] if snakemake.log else fallback_path
        )
        kwargs.update(
            {
                "handlers": [
                    # Prefer the 'python' log, otherwise take the first log for each
                    # Snakemake rule
                    logging.FileHandler(logfile),
                    logging.StreamHandler(),
                ]
            }
        )
    logging.basicConfig(**kwargs)


def mock_snakemake(rulename, **wildcards):
    """
    This function is expected to be executed from the 'scripts'-directory of '
    the snakemake project. It returns a snakemake.script.Snakemake object,
    based on the Snakefile.

    If a rule has wildcards, you have to specify them in **wildcards.

    Parameters
    ----------
    rulename: str
        name of the rule for which the snakemake object should be generated
    **wildcards:
        keyword arguments fixing the wildcards. Only necessary if wildcards are
        needed.
    """
    import os

    import snakemake as sm
    from packaging.version import Version, parse
    from addict import Dict
    from snakemake.script import Snakemake

    script_dir = Path(__file__).parent.resolve()
    assert (
        Path.cwd().resolve() == script_dir
    ), f"mock_snakemake has to be run from the repository scripts directory {script_dir}"
    os.chdir(script_dir.parent)
    for p in sm.SNAKEFILE_CHOICES:
        if os.path.exists(p):
            snakefile = p
            break
    kwargs = dict(rerun_triggers=[]) if parse(sm.__version__) > Version("7.7.0") else {}
    workflow = sm.Workflow(snakefile, overwrite_configfiles=[], **kwargs)
    workflow.include(snakefile)
    workflow.global_resources = {}
    rule = workflow.get_rule(rulename)
    dag = sm.dag.DAG(workflow, rules=[rule])
    wc = Dict(wildcards)
    job = sm.jobs.Job(rule, dag, wc)

    def make_accessable(*ios):
        for io in ios:
            for i in range(len(io)):
                io[i] = os.path.abspath(io[i])

    make_accessable(job.input, job.output, job.log)
    snakemake = Snakemake(
        job.input,
        job.output,
        job.params,
        job.wildcards,
        job.threads,
        job.resources,
        job.log,
        job.dag.workflow.config,
        job.rule.name,
        None,
    )
    # create log and output dir if not existent
    for path in list(snakemake.log) + list(snakemake.output):
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    os.chdir(script_dir)
    return snakemake


def throw_on_psspy_error(error_code: int, detail=None) -> None:
    """
       This function makes sure the script stops on PSS/E errors and helps debugging the error.

       Parameters
       ----------
       error_code: int
           error code returned by PSS/E, which can be looked up in the PSS/E API Manual
       detail
           extra information to be printed to help debugging
       """
    if error_code:
        if detail:
            logger.debug(detail)
        raise Exception("PSSPy API error: {0} (see API documentation for details)".format(error_code))


def get_start_and_end_of_year(year: int) -> Tuple[pd.Timestamp, pd.Timestamp]:
    end = pd.Timestamp(year=year, month=12, day=31, hour=23, tz="UTC")  # end of the year
    start = end - pd.DateOffset(years=1)  # end of last year (to have a period of exactly one year including 1st of Jan)

    return start, end


# get all mapped values in 2 config (carrier in pypsa and entsoe)
def get_empty_generation_df(generation_config: dict) -> pd.DataFrame:
    set1 = set(val for val in generation_config["pypsa_map"].values())
    set2 = set(val for val in generation_config["entsoe_map"].values())
    return pd.DataFrame(
        set1.union(set2),
        columns=['type']
    ).set_index('type').sort_index()


def bus_included(bus: pd.Series, config: dict) -> bool:
    if bus is None:
        return False

    if bus.bidding_zone not in config["bidding_zones"]:
        return False

    if bus.carrier != 'AC':
        return False

    if config["build_network"]["exclude_under_construction"] and bus.under_construction:
        return False

    if not bus.in_synchronous_network:
        return False

    return True


def all_combinations(list1: list, list2: list) -> list:
    l = [list(zip(each_permutation, list2)) for each_permutation in itertools.permutations(list1, len(list2))]
    return [item for sublist in l for item in sublist]  # flatten


def all_areas_entsoe(config: dict):
    bidding_zones = config["bidding_zones"]
    countries = config["countries"]

    return sorted(set(bidding_zones) | set(countries))


def get_entsoe_client(config: dict) -> EntsoePandasClient:
    import requests_cache

    requests_cache.install_cache('entsoe_api')
    return EntsoePandasClient(api_key=config["entsoe"]["security_token"])
