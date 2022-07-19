import argparse
import gzip
import math
import os
import datetime
import re
import shutil
import pandas as pd
import pyfastaq.sequences
from pyfastaq.tasks import to_fasta
from pyfastaq.sequences import file_reader as fasta_reader
import preprocessing as prep
import numpy as np
from statistics import strongly_connected_components_description

parser = argparse.ArgumentParser()
parser.add_argument("--input_description", default=None, type=str,
                    help="Path to a description of the experiment in the TSV format (columns FASTQ and CONTROL, "
                         "tab separator). In both columns, filenames (without path) should be indicated.")
parser.add_argument("--input_directory", default=None, type=str,
                    help="Path to the directory that contains ALL the inputs -- both treatment and control files.")
parser.add_argument("--working_dir", default=None, type=str, help="Path to the working directory. If directory "
                                                                  "exists, error will be raised.")
parser.add_argument("--k", default=31, type=int, help="Size of the k-mer created by BCALM. Defaults at 31.")
parser.add_argument("--minimal_abundance", default=None, type=int,
                    help="Minimal abundance of a k-mer. Should be lower for higher k-mers. "
                         "Note that lower abundance numbers lead to "
                         "higher data complexity and longer runtime."
                         "If not specified, it is set as 90th percentile in the (non-unique) k-mer counts."
                         "Percentile can be adjusted by the <--minimal_abundance_percentile> parameter (def. 75)")
parser.add_argument("--minimal_abundance_percentile", default=90.0, type=float,
                    help="Defined percentile of sufficiently abundant kmers."
                         "This calculation is time-consuming, "
                         "but comes in handy if you don't know much about the input fasta/fastq size. "
                         "BCALM also makes a number of temporary files during the calculation. "
                         "These will be removed afterwards.")
parser.add_argument("--draw", dest='draw', action='store_true',
                    help="Draw components with more than two k-mers. "
                         "Takes time, but can be useful for visual analysis of the component. Optional. "
                         "Has a time threshold, so extra large component may not be drawn.", )
parser.add_argument("--no_store_low", dest='low_store', action='store_false',
                    help="Timesaving option. Skip saving low abundant k-mers.", )
parser.add_argument("--no_controls", dest='available_controls', action='store_false', # todo check if works
                    help="Option to switch off inclusion of control experiment. "
                         "Links with empty column CONTROL in input description CSV")


def get_fasta_from_fastq(fastqfile, args):
    zipped = False
    if fastqfile[-3:] == ".gz":  # unzip
        print(f"unzipping file {fastqfile}...")
        zipped = True
        unzipped_filename = os.path.join(args.working_dir, os.path.split(fastqfile)[1].split(".")[0]) + ".fastq"
        with gzip.open(fastqfile, 'rb') as f_in:
            with open(unzipped_filename, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        unzipped_filename = args.input_fastq

    print(f"extracting FASTQ file {unzipped_filename} to FASTA...")
    path_to_fasta = os.path.join(args.working_dir, os.path.split(fastqfile)[1].split(".")[0]) + ".fasta"
    to_fasta(unzipped_filename, path_to_fasta)

    if zipped:
        print(f"removing the superfluous FASTQ {unzipped_filename} (only the original gzipped version is kept)...")
        os.remove(unzipped_filename)

    return path_to_fasta


def get_path_to_fasta(ident, args):
    path = f"{args.input_directory}/{ident}"
    if (path[-9:] == ".fastq.gz") or (path[-6:] == ".fastq"):
        path_to_fasta = get_fasta_from_fastq(path, args)
    elif path[-6:] == ".fasta":
        path_to_fasta = path
    else:
        raise Exception(f"Unrecognized file extension: {path}")

    return path_to_fasta


def pool_sequences(pooled_fasta, fasta_path_dict):
    with open(pooled_fasta, mode="wb") as writer:
        for file in fasta_path_dict:
            f = fasta_path_dict[file]
            with open(f, 'rb') as fd:
                shutil.copyfileobj(fd, writer)


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


def run_prep(args, fastapath, controlpath, scaling_coef=0.99):
    base_path = os.path.join(args.working_dir, os.path.split(fastapath)[1].split(".")[0])
    unitigs = prep.run_bcalm(fastapath, args.k, args.minimal_abundance * 2)
    shutil.move(unitigs, base_path + ".unitigs.fa")
    unitigs = base_path + ".unitigs.fa"

    moved = prep.replace_zero_index(unitigs)
    expanded = prep.expand(moved)
    jf_csv = prep.run_jellyfish_bcalm_kmers(fastapath,
                                            expanded,
                                            args.k,
                                            directory=args.working_dir,
                                            cleanup=False
                                            )

    if controlpath is not None:
        control_counts_filename = prep.get_control_counts(controlpath, args.k, directory=args.working_dir)

        control_counts = pd.read_csv(control_counts_filename, sep=' ', header=None)
        control_counts.columns = ['kmer', 'count']

        current_counts = pd.read_csv(jf_csv, sep=' ', header=None)
        current_counts.columns = ['kmer', 'count']

        total_current = current_counts['count'].sum()
        total_control = control_counts['count'].sum()

        # scale control to match the magnitudes in current counts
        scale = total_current / total_control
        scale *= scaling_coef
        control_counts['count'] = control_counts['count'] * scale
        control_counts['count'] = control_counts['count'].round()

        merged = pd.merge(current_counts, control_counts, how='left', on='kmer', suffixes=['', '_control'])
        merged = merged.fillna(0)

        merged["updated_count"] = np.fmax(merged['count'] - merged['count_control'], 0).astype(int)
        merged = merged.drop(columns=['count', 'count_control'])
        merged = merged[merged["updated_count"] > 0]
        merged.to_csv(jf_csv, header=False, index=False, sep=' ')

    if args.low_store:
        low_abund_file = prep.store_low_abundace_kmers(jf_csv,
                                                       args.minimal_abundance,
                                                       args.k,
                                                       directory=args.working_dir
                                                       )
        print(f'lowly abundant k-mers stored in {low_abund_file}')

    linked_unitigs = prep.link_jellyfish_counts(expanded, jf_csv, args.k)
    filtered = prep.filter_abundance(linked_unitigs, args.minimal_abundance, args.k)

    filtered = prep.check_edges(filtered, args.k)

    compdir = prep.divide_to_components(filtered, comp_dir=os.path.join(args.working_dir, "components"),
                                        canonical=False)
    if args.draw:
        prep.draw_components(compdir, args.minimal_abundance, canonical=args.canonical)

    return compdir


def main(args):
    for item in [args.input_description, args.input_directory]:
        if item is None:
            print("No input given, exiting.")
            exit(1)

    # create working_directory
    if args.working_dir is None:
        input_short = f"fasta={os.path.split(args.input_directory)[-1]}"
        args.working_dir = "{}-{}-{}".format(
            "motif_finding",
            datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S"),
            "_".join(
                [*("{}={}".format(re.sub("(.)[^_]*_?", r"\1", key), value) for key, value in sorted(vars(args).items())
                   if
                   key not in ["working_dir", "input_fasta", "input_fastq", "control_filename"]), input_short]
            )
        )
        print(f"The working dir was set as {args.working_dir}")
    os.makedirs(args.working_dir, exist_ok=False)

    description = pd.read_csv(args.input_description, header=0, sep="\t")
    #
    # # to fasta
    fastq2fasta = {}
    for fastq in description["fastq"]:
        fastq2fasta[fastq] = get_path_to_fasta(fastq, args)
    control2fasta = {}
    control_groups = description.groupby(by="control")
    for name, gr in control_groups:
        control2fasta[name] = get_path_to_fasta(name, args)

    pooled_fasta = f"{args.working_dir}/pooled_sequences.fasta"
    pooled_fasta = os.path.abspath(pooled_fasta)
    print(f"pooling treatment sequences --> {pooled_fasta}")
    pool_sequences(pooled_fasta, fastq2fasta)

    pooled_control = f"{args.working_dir}/pooled_control_sequences.fasta"
    pooled_control = os.path.abspath(pooled_control)
    print(f"pooling control sequences --> {pooled_control}")
    pool_sequences(pooled_control, control2fasta)

    # # calculate minimal abundance
    if args.minimal_abundance is None:
        overall_minimal_abundance = calculate_minimal_abundance(args,
                                                                pooled_fasta,
                                                                percentile=args.minimal_abundance_percentile,
                                                                cleanup=True)
        args.minimal_abundance = overall_minimal_abundance
    else:
        args.minimal_abundance = args.minimal_abundance
        # prune based on counts

    components_path = run_prep(args, pooled_fasta, pooled_control)
    strongly_connected_components_description(components_path)

    # prepare counts on individual treatment and control files
    if_filename = f"{args.working_dir}/present_kmers.fa"
    with open(if_filename, 'w') as iffile:
        if args.low_store:
            for entry in pyfastaq.sequences.file_reader(f"{args.working_dir}/low_abund_kmers.fasta"):
                iffile.write(f">{entry.id}\n{entry.seq}\n")
        for entry in pyfastaq.sequences.file_reader(f"{args.working_dir}/mers.unitigs.{args.k}.fa"):
            iffile.write(f">{entry.id}\n{entry.seq}\n")

    for item in fastq2fasta:
        path = fastq2fasta[item]
        prep.get_all_kmer_counts(path,
                                 args.k,
                                 args.working_dir,
                                 f"{args.working_dir}/{item}.csv",
                                 if_filename
                                 )
    for item in control2fasta:
        path = control2fasta[item]
        prep.get_all_kmer_counts(path,
                                 args.k,
                                 args.working_dir,
                                 f"{args.working_dir}/{item}.csv",
                                 if_filename
                                 )

    # do cleanup
    seen = os.listdir(args.working_dir)

    to_remove = [x for x in seen if (x[-5:] == "fasta") and (x != "low_abund_kmers.fasta")]
    to_remove.extend([x for x in seen if x[-3:] == "tmp"])
    to_remove.extend([x for x in seen if x[-2:] == "fa"])
    to_remove.extend([x for x in seen if x[-2:] == "jf"])
    for item in to_remove:
        path = os.path.join(args.working_dir, item)
        os.remove(path)


if __name__ == '__main__':
    args = parser.parse_args([] if "__file__" not in globals() else None)
    main(args)
