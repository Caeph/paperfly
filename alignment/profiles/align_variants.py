import argparse
import pandas as pd
import edlib
import numpy as np
from progress.bar import ChargingBar as bar

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
                    default="max",
                    type=str,
                    help="Manner in which to add the low abundant kmers to partially assembled sequences. "
                         "Options: <scaled>: add the count to all matches respectively to the original mean overall "
                         "count of the matched sequence."
                         "<max>: add everything to the most abundant (in terms of original mean overall count) matched "
                         "sequence."
                    )  # no mean recalculation is done
parser.add_argument("--k", type=int, default=None, help="K-mer length.")
parser.add_argument("--output_path",
                    type=str,
                    default="aligned_mapped.csv",
                    help="Path to the output: aligned sequences with mapped lowly abundant. A ; separated csv file "
                         "is produced."
                    )

chunksize = 10


def rolling(df, window, step):
    count = 0
    df_length = len(df)
    while count < (df_length - window):
        yield count, df[count:window+count]
        count += step
    yield count, df[count:]


def process_scaled(candidates, kmer, kmer_count):
    # index, startpos, scaled, kmer
    candidates["startpos"] = candidates["sq"].apply(
        lambda longsq: edlib.align(kmer, longsq,
                                   mode="HW", k=args.similarity_thr,
                                   task="locations"
                                   )["locations"][0][0]
    )
    scaled = candidates["mean_count"] / candidates["mean_count"].sum()
    candidates["scaled"] = np.round(scaled * kmer_count).astype(int)
    output = candidates[["index", "startpos", "scaled"]].copy()
    output = output[output["scaled"] > 0]
    output["kmer"] = kmer

    return output.values


def process_max(candidates, kmer, kmer_count):
    # index startpos, scaled, kmer
    index = candidates["mean_count"].argmax()
    target = candidates.iloc[index]

    startpos = edlib.align(kmer, target["sq"],
                           mode="HW",
                           k=args.similarity_thr,
                           task="locations"
                           )["locations"][0][0]
    return np.atleast_2d(np.array([target["index"], startpos, kmer_count, kmer]))


processing_styles = {
    "scaled" : process_scaled,
    "max" : process_max
}
# functions:
# input : candidates df, kmer, kmer count
# output : count,index,startpos start position


def main(args):
    for item in [args.assembled_profiles, args.low_abundance_kmers, args.variants_from_ec]:
        if item is None:
            print("No input given, exiting.")
            exit(1)

    # if args.align_to not in processing_styles:
    #     raise Exception("Unsupported variant of alignment.")
    processing_func = processing_styles[args.align_to]

    df_low = pd.read_csv(args.low_abundance_kmers, header=None, sep=';')
    df_low.columns = ["kmer", "counts"]
    df_low["counts"] = df_low["counts"].str.split(",").str[0].astype(int)
    df_low["key"] = 1

    df_profiles = pd.read_csv(args.assembled_profiles, header=None, sep=';')
    df_profiles.columns = ["sq", "counts"]
    df_profiles.reset_index(inplace=True)
    df_profiles["mean_count"] = df_profiles["counts"].str.split(",").apply(
        lambda lst: sum([int(x) for x in lst]) / len(lst)
    )
    # df_profiles["low_abund_indices"] = np.empty((len(df_profiles), 0)).tolist()
    df_profiles["counts"] = df_profiles["counts"].str.split(",").apply(
        lambda lst: [int(x) for x in lst]
    )
    df_profiles["key"] = 1
    print("Data loaded.")

    with bar("Processing low abundancy kmers...", max=np.ceil(len(df_low) / chunksize)) as b:
        overall_results = []

        for _, chunk in rolling(df_low, chunksize, chunksize):
            crossed = pd.merge(df_profiles, chunk, on="key", suffixes=("_sq", "_kmer"))
            crossed['combined'] = crossed[["kmer", "sq"]].apply(
                lambda row: row.values, axis=1)
            dsts = crossed['combined'].apply(
                lambda x: edlib.align(x[0], x[1],
                                      mode="HW", k=args.similarity_thr)["editDistance"]
            )
            candidates = crossed[dsts != -1].copy()
            candidates["distances"] = dsts[dsts != -1]

            cand_groups = candidates.groupby(by="kmer")
            for kmer, gr in cand_groups:
                min_dst = gr["distances"].min()
                kmer_count = gr["counts_kmer"].iloc[0]
                valid_candidates_gr = gr[gr["distances"] == min_dst].copy()

                aligned_result = processing_func(valid_candidates_gr, kmer, kmer_count)
                overall_results.append(aligned_result)

            b.next()

    c = np.vstack(overall_results)
    overall_results = pd.DataFrame(c)
    overall_results.columns = ["target_index", "startpos", "count_toadd", "kmer"]
    overall_results["target_index"] = overall_results["target_index"].astype(int)

    # this does not work
    results = pd.merge(df_profiles,
                       overall_results,
                       how="left", left_on="index", right_on="target_index")

    print(results.columns)

    # todo add the counts to add to general counts

    results = results.groupby(by="sq").agg(
        {
            'sq': lambda x: x[0],
            'counts': lambda x: x[0],
            'target_index': lambda x: ','.join(x),
            'kmer': lambda x: ','.join(x),
            # 'target_index': lambda x: ','.join(x),
        }
    )
    results.to_csv(args.output_path, sep=';', index=None)
    print("Finished.")


if __name__ == '__main__':
    args = parser.parse_args([] if "__file__" not in globals() else None)
    main(args)