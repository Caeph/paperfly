import argparse
import gzip
import math
import os
import datetime
import re
import shutil
import subprocess
import pandas as pd
from pyfastaq.tasks import to_fasta
from pyfastaq.sequences import file_reader as fasta_reader
import preprocessing as prep
# import training_weighted as hmm
import pomegranate as pm
import numpy as np
from itertools import product
import time
from statistics import strongly_connected_components_description

parser = argparse.ArgumentParser()
parser.add_argument("--input_fastq",
                    default=None,
                    type=str,
                    help="Input FASTQ file. Can be gzipped, in that case, \".gz\" suffix is expected. "
                         "If both fasta and fastq file are set, the fasta file will be used as input.")
parser.add_argument("--input_fasta",
                    default=None,
                    type=str,
                    help="Input FASTA file. If both fasta and fastq file are set, "
                         "the fasta file will be used as input. The fasta file cannot be gzipped.")
parser.add_argument("--working_dir",
                    default=None,
                    type=str,
                    help="Path to the working directory. If directory exists, error will be thrown."
                    )
parser.add_argument("--k",
                    default=31,
                    type=int,
                    help="Size of the k-mer created by BCALM. Defaults at 31."
                    )
parser.add_argument("--minimal_abundance",
                    default=None,
                    type=int,
                    help="Minimal abundance of a k-mer. Should be lower for higher k-mers. "
                         "Note that lower abundance numbers lead to "
                         "higher data complexity and longer runtime."
                         "If not specified, it is set as 75th percentile in the (non-unique) k-mer counts."
                         "Percentile can be adjusted by the <--minimal_abundance_percentile> parameter (def. 75)"
                    )
parser.add_argument("--minimal_abundance_percentile",
                    default=75.0,
                    type=float,
                    help="Defined percentile of sufficiently abundant kmers."  
                         "This calculation is time-consuming, "
                         "but comes in handy if you don't know much about the input fasta/fastq size. "
                         "BCALM also makes a number of temporary files during the calculation. "
                         "These will be removed afterwards."
                    )
parser.add_argument("--draw",
                    dest='draw',
                    action='store_true',
                    help="Draw components with more than two k-mers. "
                         "Takes time, but can be useful for visual analysis of the component. Optional. "
                         "Has a time threshold, so extra large component may not be drawn.",
                    )
parser.add_argument("--unwrap",
                    dest='unwrap',
                    action='store_true',
                    help="Unwrap merged k-mers.",
                    )


def calculate_minimal_abundance(args, fastapath, percentile=75, cleanup=False):
    print(f"calculating minimal abundance with percentile {percentile}")
    unitigs_abu = prep.run_bcalm(fastapath, args.k, 2)

    items = []
    for entry in fasta_reader(unitigs_abu):
        id_, ln, kc, km, *_ = entry.id.split(' ')
        num_vars = []
        for var in [ln, kc, km]:
            num_vars.append(float(var.split(':')[-1]))
        items.append([id_, *num_vars])

    df = pd.DataFrame(items)
    df.columns = ['id', 'length', 'count_complete', 'count_average']
    min_abu = math.ceil(np.percentile(df.count_average.values, percentile))

    if cleanup:
        os.remove(unitigs_abu)

    return int(min_abu)


def get_fasta_from_fastq(fastqfile, args):
    zipped = False
    if fastqfile[-3:] == ".gz":  # unzip
        print("unzipping file...")
        zipped = True
        unzipped_filename = os.path.join(args.working_dir, os.path.split(fastqfile)[1].split(".")[0]) + ".fastq"
        with gzip.open(args.input_fastq, 'rb') as f_in:
            with open(unzipped_filename, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        unzipped_filename = args.input_fastq

    print("extracting FASTQ file to FASTA...")
    path_to_fasta = os.path.join(args.working_dir, os.path.split(fastqfile)[1].split(".")[0]) + ".fasta"
    to_fasta(unzipped_filename, path_to_fasta)

    if zipped:
        print("removing the superfluous FASTQ (only the original gzipped version is kept)...")
        os.remove(unzipped_filename)

    return path_to_fasta


def run_prep(args, fastapath):
    base_path = os.path.join(args.working_dir, os.path.split(fastapath)[1].split(".")[0])

    if not args.canonical:
        unitigs = prep.run_bcalm(fastapath, args.k, args.minimal_abundance * 2)
    else:
        unitigs = prep.run_bcalm(fastapath, args.k, args.minimal_abundance)
    os.rename(unitigs, base_path + ".unitigs.fa")
    unitigs = base_path + ".unitigs.fa"

    if not args.canonical:
        moved = prep.replace_zero_index(unitigs)
        expanded = prep.expand(moved)
        jf_csv = prep.run_jellyfish_bcalm_kmers(fastapath,
                                                expanded,
                                                args.k,
                                                directory=args.working_dir)
        linked_unitigs = prep.link_jellyfish_counts(expanded, jf_csv, args.k)
        filtered = prep.filter_abundance(linked_unitigs, args.minimal_abundance, args.k)
        if args.unwrap:
            prep.split_unitigs(filtered, base_path+".split.fa", args.k)
            filtered = base_path+".split.fa"
    else:
        filtered = unitigs

    compdir = prep.divide_to_components(filtered, comp_dir=os.path.join(args.working_dir, "components"),
                                        canonical=args.canonical)
    if args.draw:
        prep.draw_components(compdir, args.minimal_abundance, canonical=args.canonical)

    return compdir


def print_singles(singles_path, output_path):
    with open(output_path, mode='w') as writer:
        for entry in fasta_reader(singles_path):
            count = min(int(x) for x in entry.id.split(' ')[4].split(','))
            writer.write(f"{entry.seq};{count}\n")


def run_assemblation(args, components_dir, report=5,
                     path_to_sampler="~/Projects/MatrixMotif/debruijn_pipeline/Sampler/Sampler/bin/Release/Sampler.exe"):
    # TODO deal with path in a reasonable way...
    assemblies_dir = os.path.join(args.working_dir, "assemblies")
    process = subprocess.run([f"{path_to_sampler} --input_path {components_dir} "
                              f"--samplingDepth {args.sampling_depth} "
                              f"--k {args.k} "
                              f"--output_path {assemblies_dir} "
                              f"--samplerType sequence "
                              f"--report_step {report}"
                              ],
                             shell=True)
    prep.evaluate_process(process, continue_if_bad=False)

    print("printing singles...")
    print_singles(os.path.join(components_dir, "singles.comp.fasta"),
                  os.path.join(assemblies_dir, "singles.comp.fasta.csv"))

    return assemblies_dir


def pool_assemblies(args, assemblies_dir):
    lst = []
    print("pooling...")
    for filename in os.listdir(assemblies_dir):
        df = pd.read_csv(os.path.join(assemblies_dir, filename), index_col=None, header=None, sep=';')
        lst.append(df)
    all_assembled_df = pd.concat(lst, axis=0, ignore_index=True)
    all_assembled_df.columns = ['sequence', 'count']

    all_assembled_file = os.path.join(args.working_dir, "all_assembled.csv")

    all_assembled_df.to_csv(all_assembled_file, index=False)
    print(f"pooling done, stored in {all_assembled_file}")
    return all_assembled_file


def main(args):
    # check if input exists
    if (args.input_fastq is None) and (args.input_fasta is None):
        print("No input file given, exiting...")
        return

    # create working_directory
    if args.working_dir is None:
        if args.input_fasta is not None:
            input_short = f"fasta={os.path.split(args.input_fasta)[1]}"
        else:  # got fastq
            input_short = f"fastq={os.path.split(args.input_fastq)[1]}"
        args.working_dir = "{}-{}-{}".format(
            "motif_finding",
            datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S"),
            "_".join(
                [*("{}={}".format(re.sub("(.)[^_]*_?", r"\1", key), value) for key, value in sorted(vars(args).items())
                   if
                   key not in ["working_dir", "input_fasta", "input_fastq"]), input_short]
            )
        )
        print(f"The working dir was set as {args.working_dir}")

    args.canonical = False  # assemblation is difficult
    os.makedirs(args.working_dir, exist_ok=False)
    #
    # unzip fastq, convert to fasta
    if args.input_fasta is not None:
        fastapath = args.input_fasta
    else:  # got fastq
        fastapath = get_fasta_from_fastq(args.input_fastq, args)

    # minimal abundance calculation
    if args.minimal_abundance is None:
        args.minimal_abundance = calculate_minimal_abundance(args, fastapath, percentile=args.minimal_abundance_percentile)
        print(f"Minimal abundance was set as {args.minimal_abundance}")

    # # expand canonical, divide to (weakly connected) components, optionally draw
    components_path = run_prep(args, fastapath)

    strongly_connected_components_description(components_path)

    # run on each component the C# sampler. Runs in parallel.
    # assemblies_dir = run_assemblation(args, components_path)

    # assemblies_csv = pool_assemblies(args, assemblies_dir)
    # print(f"All pseudo-assembled sequences are stored in {assemblies_dir}.")

    # learn the model on the sampled data
    # run the learnt model on the original fasta
    # cluster identified models
    # train_model_multinomial(args, assemblies_csv, fastapath)


if __name__ == '__main__':
    args = parser.parse_args([] if "__file__" not in globals() else None)
    main(args)
