from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
import os
import pandas as pd
import argparse
import numpy as np
import pyfastaq

parser = argparse.ArgumentParser()
parser.add_argument("--k", default=None, type=int,
                    help="K-mer size.")
parser.add_argument("--profiles_filename", default="aligned", type=str,
                    help="Name of the input profiles file.")
parser.add_argument("--min_abundance_mapping", default=5, type=int, help="Minimal abundance of a kmer to be mapped to "
                                                                 "partially assembled sequences.")
parser.add_argument("--low_abund_filename", default="low_abund_kmers.fasta", type=str,
                    help="Name of the input file with low abundance k-mers.")
parser.add_argument("--working_dir", default=None, type=str,
                    help="Name of the input working directory.")
parser.add_argument("--perc_identity", default=80, type=int,
                    help="Low k-mer identity percentage threshold.")
parser.add_argument("--output_filename", default="aligned_corrected.csv", type=str,
                    help="Name of the input profiles file.")


def run_blast(working_dir, low_abund_file, profiles_file, perc_identity, k):
    db = os.path.join(working_dir, 'aligned.faa')  # BLAST database
    with open(db, mode='w') as writer:
        for i, line in enumerate(open(os.path.join(working_dir, profiles_file))):
            sq, counts = line.strip('\n').split(';')
            print(f">{i}", file=writer)
            print(sq, file=writer)

    query = os.path.join(working_dir, low_abund_file)  # query sequence(s)
    blastout = os.path.join(working_dir, 'low_abund_mapped.tab')  # BLAST output

    cmd_db = NcbimakeblastdbCommandline(dbtype="nucl", input_file=db)
    cmd_db()
    print("Blast database created")

    cmd_blastx = NcbiblastnCommandline(query=query,
                                       out=blastout,
                                       outfmt=6, db=db, window_size=0, word_size=k // 3, strand="plus", num_threads=8,
                                       evalue=0.000001, perc_identity=perc_identity)
    cmd_blastx()
    print("Blast analysis finished")

    return blastout


def process_blast_output(blastout, perc_identity, k, working_dir, profiles_file, outfile):
    try:
        results = pd.read_csv(blastout, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        df_profiles = pd.read_csv(os.path.join(working_dir, profiles_file), sep=";", header=None)
        df_profiles.columns = ["sq", "counts"]
        df_profiles.index = df_profiles["sq"]
        df_profiles = df_profiles[["counts"]]
        df_profiles.to_csv(os.path.join(working_dir, outfile), header=None, sep=';')
        return
    results.columns = ["query", "target", "percentage", "length", "mismatch",
                       "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"
                       ]
    results["same"] = results["percentage"] * results["length"] / 100

    results = results[results["same"] >= (perc_identity / 100 * k)]

    gr = results.groupby("query").max("same").reset_index()[["query", "same"]]
    gr = pd.merge(gr, results, on=["query", "same"], how="inner"
                  ).groupby("query").max("length").reset_index()[["query", "same", "length"]]

    fin = pd.merge(gr, results, on=["query", "same", "length"], how="inner").drop_duplicates("query")
    fin = fin[["query", "target", "qstart", "qend", "sstart", "send", "same"]]

    fin["abs_sstart"] = np.fmax(0, fin["sstart"] - fin["qstart"])
    fin["used_length"] = k - np.abs(np.fmin(fin["sstart"] - fin["qstart"], 0))

    fin["count_to_add"] = fin["query"].str.split("_").str[-1].astype(int)

    df_profiles = pd.read_csv(os.path.join(working_dir, profiles_file), sep=";", header=None)
    df_profiles.columns = ["sq", "counts"]
    df_profiles = df_profiles.reset_index()
    merged = pd.merge(fin, df_profiles, left_on="target", right_on="index", how="right")

    unaffected = merged[merged["used_length"].isna()].copy()
    unaffected["adjusted_counts"] = unaffected["counts"]
    unaffected.index = unaffected["sq"]

    merged = merged[~merged["used_length"].isna()]
    for col in ["used_length", "count_to_add", "abs_sstart"]:
        merged[col] = merged[col].astype(int)

    merged = merged.groupby("sq").agg(
        {
            "used_length": lambda x: list(x),
            "counts": lambda x: list(x)[0],
            "count_to_add": lambda x: list(x),
            "abs_sstart": lambda x: list(x),
        }
    )

    def add(counts, abs_sstart_col, used_length_col, count_to_add_col):
        new_counts = [int(x) for x in counts.split(',')]
        for abs_sstart, used_length, count_to_add in zip(abs_sstart_col, used_length_col, count_to_add_col):
            for i in range(abs_sstart, min(abs_sstart + used_length, len(new_counts))):
                new_counts[i] += count_to_add
        return ",".join([str(x) for x in new_counts])

    merged["adjusted_counts"] = merged.apply(lambda row: add(row["counts"],
                                                             row["abs_sstart"],
                                                             row["used_length"],
                                                             row["count_to_add"]
                                                             ), axis=1
                                             )
    merged = pd.concat([merged[["adjusted_counts"]], unaffected[["adjusted_counts"]]])
    merged.to_csv(os.path.join(working_dir, outfile), header=None, sep=';')


def main(args):
    working_dir = args.working_dir
    profiles_file = args.profiles_filename
    low_abund_file = args.low_abund_filename
    outfile = args.output_filename

    perc_identity = args.perc_identity
    k = args.k

    low_abund_file_filtered = low_abund_file + ".fil"
    with open(os.path.join(working_dir, low_abund_file_filtered), mode='w') as writer:
        for entry in pyfastaq.sequences.file_reader(os.path.join(working_dir, low_abund_file)):
            count = int(entry.id.split("_")[-1])
            if count >= args.min_abundance_mapping:
                writer.write(f">{entry.id}\n{entry.seq}\n")

    blastout = run_blast(working_dir, low_abund_file_filtered, profiles_file, perc_identity, k)
    if blastout is None:
        return
    process_blast_output(blastout, perc_identity, k, working_dir, profiles_file, outfile)

    seen = os.listdir(args.working_dir)

    to_remove = [x for x in seen if (x[-3:] in ["faa", "nhr", "nin", "nsq", "fil"])]
    # to_remove.extend([x for x in seen if x[-3:] == "tmp"])
    for item in to_remove:
        path = os.path.join(working_dir, item)
        os.remove(path)


if __name__ == '__main__':
    args = parser.parse_args([] if "__file__" not in globals() else None)
    main(args)
