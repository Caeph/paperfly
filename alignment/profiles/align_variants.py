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
                    default=7,
                    type=int,
                    help="Maximal number of nucleotides that can be different from the assembled sequence "
                         "for the low abundancy kmer to be mapped to it."
                    )
parser.add_argument("--align_to",
                    default="scaled",
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
        yield count, df[count:window + count]
        count += step
    yield count, df[count:]


def process_scaled(candidates, kmer, kmer_count):
    # index of the target sq, startpos of the kmer, scaled count to add, kmer, cigar of alignment
    candidates["startpos"] = candidates["sq"].apply(
        lambda longsq: edlib.align(kmer, longsq,
                                   mode="HW", k=args.similarity_thr,
                                   task="locations"
                                   )["locations"][0][0]
    )
    candidates["cigar"] = candidates["sq"].apply(
        lambda longsq: edlib.align(
            kmer, longsq, mode="HW", k=args.similarity_thr, task="path")["cigar"]
    )
    scaled = candidates["mean_count"] / candidates["mean_count"].sum()
    candidates["scaled"] = np.round(scaled * kmer_count).astype(int)
    output = candidates[["index", "startpos", "scaled", "cigar"]]
    output = output[output["scaled"] > 0]
    output["kmer"] = kmer

    return output.values


def process_max(candidates, kmer, kmer_count):
    # index of the target sq, startpos of the kmer, scaled count to add, kmer
    index = candidates["mean_count"].argmax()
    target = candidates.iloc[index]

    startpos = edlib.align(kmer, target["sq"],
                           mode="HW",
                           k=args.similarity_thr,
                           task="locations"
                           )["locations"][0][0]

    differences = edlib.align(kmer, target["sq"], mode="HW", k=args.similarity_thr, task="path")["cigar"]

    return np.atleast_2d(np.array([target["index"], startpos, kmer_count, kmer, differences]))


processing_styles = {
    "scaled": process_scaled,
    "max": process_max
}


def main(args):
    for item in [args.assembled_profiles, args.low_abundance_kmers, args.variants_from_ec]:
        if item is None:
            print("No input given, exiting.")
            exit(1)

    if args.align_to not in processing_styles:
        raise Exception("Unsupported variant of alignment.")

    processing_func = processing_styles[args.align_to]

    # load lowly abundant kmers
    df_low = pd.read_csv(args.low_abundance_kmers, header=None, sep=';')
    df_low.columns = ["kmer", "counts"]
    df_low["counts"] = df_low["counts"].str.split(",").str[0].astype(int)

    # load assembled profiles
    df_profiles = pd.read_csv(args.assembled_profiles, header=None, sep=';')
    df_profiles.columns = ["sq", "counts"]
    df_profiles.reset_index(inplace=True)
    df_profiles["mean_count"] = df_profiles["counts"].str.split(",").apply(
        lambda lst: sum([int(x) for x in lst]) / len(lst)
    )
    df_profiles["counts"] = df_profiles["counts"].str.split(",").apply(
        lambda lst: [int(x) for x in lst]
    )
    df_profiles["key"] = 1
    df_profiles["current_processed"] = ""
    df_profiles["combined"] = ""
    df_profiles["variability"] = ""
    print("Data loaded.")

    with bar("Processing low abundancy kmers...", max=np.ceil(len(df_low) / chunksize)) as b:
        overall_results = []

        for kmer, kmer_count in df_low.values:
            df_profiles["current_processed"] = kmer
            df_profiles['combined'] = df_profiles[["current_processed", "sq"]].apply(
                lambda row: row.values, axis=1)
            dsts = df_profiles['combined'].apply(
                lambda x: edlib.align(x[0], x[1],
                                      mode="HW", k=args.similarity_thr)["editDistance"]
            )
            # candidates = df_profiles[dsts != -1]  #.copy()
            applicable = dsts[dsts != -1]
            if len(applicable) == 0:
                continue

            # get best fit
            mindst = applicable.min()
            candidates = df_profiles[dsts == mindst].copy()

            aligned_result = processing_func(candidates, kmer, kmer_count)
            overall_results.append(aligned_result)

    overall_results = pd.DataFrame(np.vstack(overall_results),
                                   columns=["target_sq_index", "startpos", "count_to_add", "kmer", "cigar"]
                                   )
    for col in ["target_sq_index", "startpos", "count_to_add"]:
        overall_results[col] = overall_results[col].astype(int)
    overall_results["single_entry"] = overall_results["startpos"].astype(str) + \
                                      "," + overall_results["count_to_add"].astype(str) + "," \
                                      + overall_results["kmer"] + "," + overall_results["cigar"]
    overall_results["additions"] = overall_results[["count_to_add", "startpos"]].values.tolist()
    groups = overall_results.groupby(by="target_sq_index").agg(
        {
            "single_entry": lambda x: ";".join([str(item) for item in x]),
            "additions": lambda x: list(x)
        }
    )
    for index, addition, variability in zip(groups.index, groups["additions"], groups["single_entry"]):
        df_profiles.loc[index, "variability"] = variability
        for count, start in addition:
            current = df_profiles.iloc[index]["counts"]  # is written to df via reference
            for i in range(start, min(len(current), start+args.k)):
                current[i] += count

    df_profiles = df_profiles[["sq","counts","variability"]]
    df_profiles.to_csv(args.output_path, sep=";")
    print("Finished.")


if __name__ == '__main__':
    args = parser.parse_args([] if "__file__" not in globals() else None)
    main(args)
