import argparse

from ripta.config import read_config
from ripta.io import process_vcfs
from ripta.results import compile_results


def main():
    """RIPTA command-line interface.

    Args:
        None

    Returns:
        None
    """
    parser = argparse.ArgumentParser(
        description="RIPTA: Re-evaluating Introgression site Patterns Through Ancestral alleles.",
    )
    # Path to the YAML configuration file.
    parser.add_argument(
        "-y",
        "--yaml_path",
        required=True,
        type=str,
        action="store",
        help="Path to the YAML configuration file.",
    )
    args = parser.parse_args()

    # Load and validate the configuration file.
    config = read_config(args.yaml_path)
    # Process all VCF files and return per-block allele count sums.
    blocks = process_vcfs(config)
    # Aggregate block results, compute statistics, and write the CSV output.
    compile_results(config, blocks)
