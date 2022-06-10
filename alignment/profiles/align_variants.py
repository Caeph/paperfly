import argparse
import pandas as pd
import edlib
import numpy as np

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
parser.add_argument("--similarity_thr",
                    default=5,
                    type=int,
                    help="Maximal number of nucleotides that can be different from the assembled sequence "
                         "for the low abundancy kmer to be mapped to it."
                    )
parser.add_argument("--align_to",
                    default="scaled",
                    type=str,
                    help="Manner in which to add the low abundant kmers to partially assembled sequences. "
                         "Options: <scaled>: add the count to all matches respectively to the mean overall "
                         "count of the matched sequence."
                         "<max>: add everything to the most abundant (in terms of mean overall count) matched sequence."
                    )
parser.add_argument("--k", type=int)  # todo desc


def process_max(df_profiles, min_dst, alignable, kmer_count, kmer):
    ...
    # todo
    # candidates = df_profiles[alignable == min_dst]
    #if len(candidates) > 0:
    #    to_ = candidates[alignable == min_dst]["mean_count"]
    #    index = candidates.iloc[to_].name
    ##    new_val = list(df_profiles.iloc[index]["low_abund_indices"])
    #    new_val.append(index)
    #    df_profiles.loc[index, "low_abund_indices"] = new_val


def process_scaled(df_profiles, min_dst, alignable, kmer_count, kmer):
    candidates = df_profiles[alignable == min_dst].copy()
    scaled = candidates["mean_count"] / candidates["mean_count"].sum()
    scaled = np.round(scaled * kmer_count).astype(int)

    # find index where to add the count
    candidates["toadd"] = scaled
    candidates["startpos"] = candidates["sq"].apply(
        lambda longsq : edlib.align(
            kmer, longsq, task="locations",
            mode="HW", k=args.similarity_thr)["locations"][0][0]
    )

    # just adding to the start of the kmer, not minding the gaps
    for i, row in candidates.iterrows():
        c = row["counts"]
        for i in range(row["startpos"], row["startpos"]+args.k):
            c[i] += row["toadd"]
        row["counts"] = c

    # add to the dataframe
    # todo

    return df_profiles


processing_styles = {
    "scaled" : process_scaled,
    "max" : process_max
}


def main(args):
    for item in [args.assembled_profiles, args.low_abundance_kmers, args.variants_from_ec]:
        if item is None:
            print("No input given, exiting.")
            exit(1)

    if args.align_to not in processing_styles:
        raise Exception("Unsupported variant of alignment.")
    processing_func = processing_styles[args.align_to]

    df_low = pd.read_csv(args.low_abundance_kmers, header=None, sep=';')
    df_low.columns = ["kmer", "counts"]
    df_low["counts"] = df_low["counts"].str.split(",").str[0].astype(int)

    df_profiles = pd.read_csv(args.assembled_profiles, header=None, sep=';')
    df_profiles.columns = ["sq", "counts"]
    df_profiles["mean_count"] = df_profiles["counts"].str.split(",").apply(
        lambda lst: sum([int(x) for x in lst]) / len(lst)
    )
    df_profiles["low_abund_indices"] = np.empty((len(df_profiles), 0)).tolist()
    df_profiles["counts"] = df_profiles["counts"].str.split(",").apply(
        lambda lst: [int(x) for x in lst]
    )

    for i, row in df_low.iterrows():
        kmer = row['kmer']
        kmer_count = row['counts']
        alignable = df_profiles["sq"].apply(
            lambda longsq: edlib.align(kmer, longsq,
                                       mode="HW", k=args.similarity_thr)["editDistance"]
        )
        if (alignable != -1).sum() == 0:  # cannot be aligned, no sufficiently close match
            continue

        min_dst = alignable[alignable != -1].min()
        df_profiles = processing_func(df_profiles, min_dst, alignable, kmer_count, kmer)

    used_indications = df_profiles["low_abund_indices"].str.len() > 0
    used = df_profiles[used_indications]


if __name__ == '__main__':
    args = parser.parse_args([] if "__file__" not in globals() else None)
    main(args)