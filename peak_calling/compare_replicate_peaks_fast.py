import pandas as pd
import numpy as np
from multiprocessing import Pool
import sys
from pyfastaq.sequences import file_reader as fasta_reader
import os
import edlib
from progress.bar import IncrementalBar as Bar


def get_longer_shorter(sq1, sq2):
    if len(sq1) < len(sq2):
        return sq1, sq2, 0
    else:
        return sq2, sq1, 1


def get_distance(sq1, sq2, rsq1, rsq2):
    query, target, f = get_longer_shorter(sq1, sq2)
    r_query, r_target, f = get_longer_shorter(rsq1, rsq2)
    min_dst = np.inf
    reverse_q = None
    reverse_t = None
    for q, t, rq, rt in zip([query, query, r_query], [target, r_target, r_target], [0, 0, 1], [0, 1, 1]):
        dst = edlib.align(
            q, t, mode='HW', k=len(query) * partner_threshold, task='distance'
        )["editDistance"]

        if dst == -1:
            dst = np.inf
        if dst == 0:
            return 0, f, rq, rt

        if min_dst > dst:
            min_dst = dst
            reverse_q, reverse_t = rq, rt

    return min_dst / len(query), f, reverse_q, reverse_t


def get_variability(group):
    flags = {
        0: lambda seen, seen_r, mapped, mapper_r: (mapped, mapper_r, seen, seen_r),  # mapped is query
        1: lambda seen, seen_r, mapped, mapper_r: (seen, seen_r, mapped, mapper_r),  # mapped is target
    }

    reverse = {
        0: lambda seq, rev: seq,
        1: lambda seq, rev: rev,
    }
    # mapped to sq is universal for the group
    pval = pd.concat([group["pval"], group["mapped_to_pval"]]).max()
    old_variability = np.hstack(group["mapped_to_variability"].values)

    def get_correct_qt(seen, seen_r, mapped, mapper_r, flag, rq_flag, rt_flag, blurr=5):
        q, qr, t, tr = flags[flag](seen, seen_r, mapped, mapper_r)
        corr_q = reverse[rq_flag](q, qr)
        corr_t = reverse[rt_flag](t, tr)
        A = edlib.align(corr_q, corr_t, task='path', mode='HW')
        from_, to_ = A["locations"][0]
        from_ = max(0, from_ - blurr)
        to_ += blurr

        replicated = corr_t[from_:to_]

        if A["editDistance"] == 0:
            variant = None
        else:
            variant = corr_q

        return replicated, variant

    alignments = group.apply(lambda row: get_correct_qt(
        row["seq"], row["reverse"], row["mapped_to_sq"], row["mapped_to_reverse"],
        row["qt_flag"], row["reverse_query"], row["reverse_target"]
    ), axis=1)

    # if more than one replicated -- pick the shortest from replicated, throw rest to variants
    if len(alignments) == 1:
        replicated = alignments.values[0][0]
        variants = alignments.apply(lambda x: x[1]).values
    else:
        all_replicated = alignments.apply(lambda x: x[0])
        picked = all_replicated.str.len().argmin()
        replicated = all_replicated.values[picked]
        variants = alignments.apply(lambda x: x[1]).values
        # add not picked in replicated to variants
        not_picked = [x for x in range(len(alignments)) if x != picked]
        variants = np.hstack([all_replicated.values[not_picked], variants])

    variants = [x for x in variants if x is not None]
    variants = list({*variants, *old_variability})

    d = {
        "seq": replicated,  # one for every alignment (list)
        "pval": pval,  # only one overall
        "variability": variants,
    }
    return pd.Series(d)
    # return replicated_sq, variants


def reverse_series(series):
    return series.str.replace(
        'A', 'X', regex=False).str.replace(
        'T', 'A', regex=False).str.replace(
        'X', 'T', regex=False).str.replace(
        'C', 'X', regex=False).str.replace(
        'G', 'C', regex=False).str.replace(
        'X', 'C', regex=False).str[::-1]


class SeqStorage:
    def __init__(self, df):
        self.df = df

    def calculate(self, visible):
        visible_sq, visible_rev = visible
        results = df.apply(lambda row: get_distance(visible_sq, row["seq"],
                                                    visible_rev, row["reverse"]), axis=1)
        distances = results.apply(lambda res: res[0])
        flags = results.apply(lambda res: res[1])
        reverse_q = results.apply(lambda res: res[2])
        reverse_t = results.apply(lambda res: res[3])
        return distances, flags, reverse_q, reverse_t


partner_threshold = float(sys.argv[1])
output_file = sys.argv[2]
called_peaks = sys.argv[3:]  # paths to files
threads = 8


def batcher(items, size=threads):
    seen = 0
    while seen < len(items):
        batch = items[seen:seen+threads]
        yield batch
        seen += threads


if len(called_peaks) == 1:
    print("only one file is available -- no replicate analysis can be done.")
    exit(1)

currently_visible = pd.DataFrame([[entry.id, entry.seq] for entry in fasta_reader(called_peaks[0])])
currently_visible.columns = ["id", "seq"]
currently_visible["reverse"] = currently_visible["seq"].str.replace(
    'A', 'X', regex=False).str.replace(
    'T', 'A', regex=False).str.replace(
    'X', 'T', regex=False).str.replace(
    'C', 'X', regex=False).str.replace(
    'G', 'C', regex=False).str.replace(
    'X', 'C', regex=False).str[::-1]
currently_visible["pval"] = currently_visible["id"].str.split(
    "pvalues").str[-1].str.replace("_", "").str.replace("<", "").astype(float)
currently_visible["variability"] = [[] for x in range(len(currently_visible))]

# currently_visible = currently_visible.head(111)
currently_visible.drop_duplicates("seq")


# prep done
all_lvls = len(called_peaks) - 1
counter = 1
for peaksfile in called_peaks[1:]:
    df = pd.DataFrame([[entry.id, entry.seq] for entry in fasta_reader(peaksfile)])
    df.columns = ["id", "seq"]
    df["reverse"] = reverse_series(df["seq"])
    df["pval"] = df["id"].str.split("pvalues").str[-1].str.replace("_", "").str.replace("<", "").astype(float)

    min_distances = np.zeros(len(df)) + np.inf
    min_indices_in_visible = - np.ones(len(df))
    flags_in_visible = - np.ones(len(df))
    reverse_query, reverse_target = - np.ones(len(df)), - np.ones(len(df))

    storage = SeqStorage(df)

    with Bar(f"comparing with {peaksfile}", max=len(currently_visible)) as bar:
        for batch_i, batch in enumerate(batcher(list(zip(currently_visible["seq"].values,
                                                         currently_visible["reverse"].values)))):
            with Pool(threads) as pool:
                batch_results = pool.map(storage.calculate, batch)
            for in_batch_i, item in enumerate(batch_results):
                distances, flags, reverse_q, reverse_t = item
                update = np.fmin(min_distances, distances)
                updated = min_distances != update
                min_indices_in_visible[updated] = batch_i + in_batch_i
                flags_in_visible[updated] = flags[updated]
                reverse_query[updated] = reverse_q[updated]
                reverse_target[updated] = reverse_t[updated]
                min_distances = update

                bar.next()
    partnered = np.where(min_distances != np.inf)[0]
    if len(partnered) == 0:
        print("No replicated peaks were found.")
        exit(0)
    replicated = df.iloc[partnered].copy()
    replicated["dst"] = min_distances[partnered]
    replicated["qt_flag"] = flags_in_visible[partnered].astype(int)
    replicated["reverse_query"] = reverse_query[partnered].astype(int)
    replicated["reverse_target"] = reverse_target[partnered].astype(int)

    replicated["mapped_to_index"] = min_indices_in_visible[partnered].astype(int)
    replicated["mapped_to_pval"] = replicated["mapped_to_index"].apply(lambda index:
                                                                       currently_visible.iloc[index]["pval"])
    replicated["mapped_to_sq"] = replicated["mapped_to_index"].apply(lambda index:
                                                                     currently_visible.iloc[index]["seq"])
    replicated["mapped_to_reverse"] = reverse_series(replicated["mapped_to_sq"])
    replicated["mapped_to_variability"] = replicated["mapped_to_index"].apply(lambda index:
                                                                              currently_visible.iloc[index][
                                                                                  "variability"])
    groups = replicated.groupby(by="mapped_to_index").apply(get_variability)
    currently_visible = groups.reset_index()
    currently_visible["reverse"] = reverse_series(currently_visible["seq"])

    print(f"{counter} file of {all_lvls} processed...")
    counter += 1

print(f"replicated analysis done, writing results to {output_file}")
with open(output_file, mode='w') as writer:
    for i, row in currently_visible.iterrows():
        overall_counts = 0
        pval = row["pval"]
        seq = row["seq"]
        variants = row["variability"]
        print(f">{overall_counts}_peak:{i}_pval:{pval}", file=writer)
        print(seq, file=writer)

        overall_counts += 1

        for variant_sq in variants:
            print(f">{overall_counts}_peak:{i}_pval:{pval}", file=writer)
            print(variant_sq, file=writer)
            overall_counts += 1
print("writing done")
