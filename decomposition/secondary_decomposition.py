import preprocessing as prep
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--corrected_graph",
                    default=None,
                    type=str,
                    help="Path to input: fasta of a corrected graph."
                    )
parser.add_argument("--corrected_components",
                    default=None,
                    type=str,
                    help="Path to output: directory name of corrected components."
                    )
parser.add_argument("--draw",
                    dest='draw',
                    action='store_true',
                    help="Draw components with more than two k-mers. "
                         "Takes time, but can be useful for visual analysis of the component. Optional. "
                         "Has a time threshold, so extra large component may not be drawn.",
                    )


def main(args):
    corrected_file = args.corrected_graph
    output_dir = args.corrected_components

    if (corrected_file is None) or (output_dir is None):
        raise Exception("Insufficient input -- either input or output is missing.")

    prep.divide_to_components(corrected_file, comp_dir=output_dir)
    if args.draw:
        prep.draw_components(output_dir, 1, f"{output_dir}_pictures")


if __name__ == '__main__':
    args = parser.parse_args([] if "__file__" not in globals() else None)
    main(args)
