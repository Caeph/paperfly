import sys
import os
import pandas as pd
from pyfastaq.sequences import file_reader as fasta_reader
import edlib
import numpy as np
from progress.bar import ChargingBar
from multiprocessing import Pool
from itertools import groupby

partner_threshold = float(sys.argv[1])
output_file = sys.argv[2]
called_peaks = sys.argv[3:]  # paths to files

# friend_counter = 0
# friend_changers = {}

blurr = 1

wildcard_def = [
            (".", "A"),
            (".", "T"),
            (".", "C"),
            (".", "G"),
        ]


def get_longer_shorter(sq1, sq2):
    if len(sq1) < len(sq2):
        query = sq1
        target = sq2
    else:
        query = sq2
        target = sq1
    return query, target


def get_distance(sq1, sq2, wildcards=False):
    query, target = get_longer_shorter(sq1, sq2)
    if wildcards:
        additionalEqualities = wildcard_def
    else:
        additionalEqualities = None
    dst = edlib.align(query, target, mode='HW', task='distance', additionalEqualities=additionalEqualities)["editDistance"]
    return dst / len(query)


def get_query_target(sq1, revsq1, revsq2, sq2, how=None):
    S1, S2 = get_forrev(sq1, sq2, revsq1, revsq2, how=None)

    query, target = get_longer_shorter(S1, S2)
    return query, target


def get_forrev(sq1, sq2, revsq1, revsq2, how=None):
    if how is None:
        dsts = np.array([
            get_distance(sq1, sq2),
            get_distance(revsq1, sq2),
            get_distance(sq1, revsq2)
        ])
        how = dsts.argmin()

    if how == 0:
        S1 = sq1
        S2 = sq2
    elif how == 1:
        S1 = sq1
        S2 = revsq2
    else:
        S1 = revsq1
        S2 = sq2

    return S1, S2


class Seqstorage:
    def __init__(self, seq_array, wildcards):
        self.seq_array = seq_array
        self.wildcards = wildcards

    def get_sequences(self, seq):
        a = self.seq_array.apply(lambda sq2: get_distance(seq, sq2, wildcards=self.wildcards)).values
        return a


def get_best_matches(df1seqs, df2seqs, wildcards):  # pandas series
    with Pool(8) as pool:
        seqs = Seqstorage(df2seqs, wildcards)
        result = pool.map(seqs.get_sequences, df1seqs.values)
    result = np.vstack(result)
    min_index = np.argmin(result, axis=1)
    score = np.min(result, axis=1)

    return min_index, score


def get_consensus(nice):
    q = nice['query_aligned']
    t = nice['target_aligned']
    what = {
        '|': lambda i: q[i],
        '-': lambda i: (t[i] + q[i]).replace("-", ""),
        '.': lambda i: '.'
    }
    consensus = ''.join([what[item](i) for i, item in enumerate(nice["matched_aligned"])])
    return consensus


def get_merged_sequences(nodupsqs, dupsqs, N=blurr):
    maps = {}

    dict_md = {}
    for i1, i2 in zip(nodupsqs.index, dupsqs.index):
        if i2 in dict_md:
            dict_md[i2].append(i1)
        else:
            dict_md[i2] = [i1]

    def get_alignment(sq, rev, scaffoldsq, scaffoldrev):
        query, target = get_query_target(sq,rev,scaffoldrev,scaffoldsq)
        if target == scaffoldsq:
            reverse = False
        else:
            reverse = True
        align = edlib.align(query, target, task='path', mode='HW')
        return align["locations"][0], reverse

    merged_sqs = []
    for k,v in dict_md.items():
        df = pd.concat([dupsqs.loc[[k]], nodupsqs.loc[v]]).reset_index().drop(columns='index')
        maxlen_i = df["seq"].str.len().argmax()
        scaffold_row = df.iloc[maxlen_i]
        # return only replicated segments with reasonable widening
        others = df[df.index != maxlen_i]
        # a[align["locations"][0][0]:align["locations"][0][1] + 1] = 1
        alignments = others.apply(lambda row: get_alignment(row["seq"],
                                                            row["reverse"],
                                                            scaffold_row["seq"],
                                                            scaffold_row["reverse"],
                                                            ),
                                  axis=1)
        seen = np.zeros(len(scaffold_row["seq"]))
        for loc, reverse in alignments:
            if reverse:
                seen = seen[::-1]
                seen[loc[0]:loc[1] + 1] = 1
                seen = seen[::-1]
            else:
                seen[loc[0]:loc[1]+1] = 1
        seen = np.ceil(np.fmin(np.convolve(seen, np.ones(N) / N, mode='same'), 1))  # get sufficiently long windows
        groups = np.array([len(list(group)) for st, group in groupby(seen)])
        groups_states = np.array([st for st, group in groupby(seen)])
        starts = np.hstack([[0], groups.cumsum()])[:-1]
        ends = np.hstack([groups.cumsum(), [len(seen)]])[:-1]

        pval = "pvalues_<" + str(df["id"].str.split("__").str[-1].str.replace("pvalues_", "").str.replace("<", "").astype(float).max())
        ID = ":".join(set(df["id"].str.split("__pvalues").str[0])) + "__" + pval

        res = [[ID+f"__{from_}_{to_}", scaffold_row["seq"][from_:to_]] for from_, to_ in zip(starts[groups_states != 0], ends[groups_states != 0])]
        for from_, to_ in zip(starts[groups_states != 0], ends[groups_states != 0]):
            maps[ID + f"__{from_}_{to_}"] = [scaffold_row["seq"][from_:to_]]

        for from_, to_ in zip(starts[groups_states != 0], ends[groups_states != 0]):  # in scaffold
            for item, aligned in zip(others.iterrows(), alignments):
                i, row = item
                loc, reverse = aligned
                if reverse:
                    maps[ID + f"__{from_}_{to_}"].append(row["seq"])
                else:
                    maps[ID + f"__{from_}_{to_}"].append(row["reverse"])
        merged_sqs.extend(res)
    return pd.DataFrame(merged_sqs, columns=["id", "seq"]), maps


def merge_pair(pair):
    if len(pair) < 2:
        return pair[0]
    # essentially do cross join but in limited space
    df1, df2 = pair
    df1["sqlen"] = df1["seq"].str.len()
    df2["sqlen"] = df2["seq"].str.len()
    df1 = df1.sort_values(by="sqlen", ascending=False)
    df2 = df2.sort_values(by="sqlen", ascending=False)

    ff_index, ff_score = get_best_matches(df1["seq"], df2["seq"], True)
    rf_index, rf_score = get_best_matches(df1["reverse"], df2["seq"], True)
    fr_index, fr_score = get_best_matches(df1["seq"], df2["reverse"], True)

    S = np.vstack([ff_score, rf_score, fr_score]).T
    I = np.vstack([ff_index, rf_index, fr_index]).T

    same = (S.min(axis=1) < partner_threshold)
    Ssame, Isame = S[same,:], I[same,:]
    matchd_df1 = np.arange(len(df1))[same]
    matchd_df2 = np.array([i[ind] for i, ind in zip(Isame, Ssame.argmin(axis=1))])
    merges_sqs, maps = get_merged_sequences(df1.iloc[matchd_df1][["id", "seq", "reverse"]], # no duplicities
                         df2.iloc[matchd_df2][["id", "seq", "reverse"]], # yes duplicities
                         )
    merges_sqs["reverse"] = merges_sqs["seq"].str.replace('A', 'X', regex=False).str.replace('T', 'A', regex=False
                                            ).str.replace('X', 'T', regex=False).str.replace('C', 'X', regex=False
                                            ).str.replace('G', 'C', regex=False).str.replace('X', 'C', regex=False
                                            ).str[::-1]
    return merges_sqs, maps


def merge_maps(newmaps, oldmaps):
    new_to_old = {}
    for old_key in oldmaps:
        old_key_ident = old_key.split("__pvalues")[0]
        matches = [k for k in newmaps.keys() if old_key_ident in k]
        for match in matches:
            if match in new_to_old:
                new_to_old[match].append(old_key)
            else:
                new_to_old[match] = [old_key]

    result = {}
    for item, parts in new_to_old.items():
        seqs = []
        for x in parts:
            seqs.extend(oldmaps[x])
        result[item] = seqs

    return result


dfs = []
for i, peaksfile in enumerate(called_peaks):
    df = pd.DataFrame([[entry.id, entry.seq] for entry in fasta_reader(peaksfile)])
    df.columns = ["id", "seq"]
    df["reverse"] = df["seq"].str.replace('A', 'X', regex=False
                                          ).str.replace('T', 'A', regex=False
                                                        ).str.replace('X', 'T', regex=False
                                                                      ).str.replace('C', 'X', regex=False
                                                                                    ).str.replace('G', 'C', regex=False
                                                                                                  ).str.replace('X',
                                                                                                                'C',
                                                                                                                regex=False).str[
                    ::-1]
    df["id"] = str(i)+df["id"].str.replace("peak", "P", regex=False).str.replace("assembled_sq", "A", regex=False)
    dfs.append(df)

# map highly similar peaks on one another
# tree structure --> logarithmic number of mappings (all to all not needed)

to_merge = dfs

lvl = 0
merged_variants = None
while len(to_merge) > 1:
    # to_merge = [merge_pair(to_merge[i:i + 2], i // 2, lvl) for i in range(len(to_merge)) if i % 2 == 0]
    new_merged = []
    new_maps = {}
    for i in range(len(to_merge)):
        if i % 2 == 0:
            print(f"comparing pair {(i // 2) + 1} of {len(to_merge) // 2} on level {lvl}")
            new, maps = merge_pair(to_merge[i:i + 2])
            new_merged.append(new)
            new_maps.update(maps)

    to_merge = new_merged
    if merged_variants is None:
        merged_variants = new_maps
    else:
        merged_variants = merge_maps(new_maps, merged_variants)
    lvl += 1

replicated = to_merge[0]

print("writing to output file...")
with open(output_file, mode='w') as writer:
    i=0
    for id, seqs in merged_variants.items():
        for seq in seqs:
            print(f">{i}_{id}\n{seq}", file=writer)
            i+=1