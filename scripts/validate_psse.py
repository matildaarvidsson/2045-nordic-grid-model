import logging
from _helpers import configure_logging

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("validate_psse", case="high")
    configure_logging(snakemake)