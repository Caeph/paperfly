import argparse
import pandas as pd
import edlib

parser = argparse.ArgumentParser()
parser.add_argument("--assembled_profiles",
                    default=None,
                    type=str,
                    help="Path to the input file: output of exact match mappings."
                    )
parser.add_argument("--low_abundance_kmers",
                    default=None,
                    type=str,
                    help="Path to the input file: low abundance kmers"
                    )
parser.add_argument("--variants_from_ec",
                    default=None,
                    type=str,
                    help="Path to the input file: tsv from variants found during error corrections."
                    )


def main(args):
    # todo
    # align variants with count excluded (IF POSSIBLE -- if fragments are too short, throw out)
    # indicate variants for assembled sequences
    # align low abundant kmers (even if overlaps)
    ...


if __name__ == '__main__':
    args = parser.parse_args([] if "__file__" not in globals() else None)
    main(args)