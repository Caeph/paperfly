import sys
import os
import pandas as pd
from pyfastaq.sequences import file_reader as fasta_reader
import edlib
import numpy as np
from progress.bar import ChargingBar


partner_threshold = float(sys.argv[1])
output_file = sys.argv[2]
called_peaks = sys.argv[3:]  # paths to files


def get_distance(sq1, sq2):
    if len(sq1) < len(sq2):
        query = sq1
        target = sq2
    else:
        query = sq1
        target = sq2
    dst = edlib.align(query, target, mode='HW', task='distance')["editDistance"]
    return dst / len(query)


def get_query_target(sq1, revsq1, revsq2, sq2, how):
    if how == 0:
        S1 = sq1
        S2 = sq2
    elif how == 1:
        S1 = sq1
        S2 = revsq2
    else:
        S1 = revsq1
        S2 = sq2

    if len(S1) < len(S2):
        query = S1
        target = S2
    else:
        query = S1
        target = S2
    return query, target


def merge_pair(pair, pairno, lvlno):
    if len(pair) < 2:
        return pair[0]

    cols = ["forfor", "forrev", "revfor"]
    # essentially do cross join but in limited space
    df1, df2 = pair
    df1["sqlen"] = df1["seq"].str.len()
    df2["sqlen"] = df2["seq"].str.len()
    df1 = df1.sort_values(by="sqlen", ascending=False)  # pick a partner for longest first

    partnered = []
    with ChargingBar(f"comparing pair {pairno} on level {lvlno}", max=len(df1)) as cbar:
        for i, row in df1.iterrows():
            sq1, revsq1 = row["seq"], row["reverse"]

            # sq1 forward
            df2[f"forfor"] = df2["seq"].apply(lambda sq2: get_distance(sq1, sq2))
            df2[f"forrev"] = df2["reverse"].apply(lambda sq2: get_distance(sq1, sq2))

            # sq1 reverse
            df2[f"revfor"] = df2["seq"].apply(lambda sq2: get_distance(revsq1, sq2))

            # get best partner for
            minimal = df2[cols].values.min()

            cbar.next()

            if minimal > partner_threshold:
                continue

            # get line indexes in df2 for best partners for sq1
            possible_matches = np.unique(np.argwhere(df2[cols].values == minimal)[:, 0])

            if len(possible_matches) > 1:
                # get the longest partner
                possible_matches = possible_matches[np.argsort(df2.iloc[possible_matches]["sqlen"].values)]

            partner = possible_matches[0]
            partner_row = df2.iloc[partner]
            align_type = np.argmin(partner_row[cols].values)

            if minimal != 0:  # a variant is seen
                query, target = get_query_target(sq1, revsq1,
                                                 partner_row["seq"], partner_row["reverse"],
                                                 align_type)
                cigar = edlib.align(query, target, mode='HW', task='path')["cigar"]
            else:
                cigar = ""  # empty cigar <=> complete match

            if "variants" in df2.columns:
                variants = partner_row["variants"]
            else:
                variants = []
            if "variants" in row:
                variants.extend(row["variants"])

            if len(row["seq"]) > len(partner_row["seq"]):
                new_id = row["id"]
                sec_id = partner_row["id"]
                seq_to_use = row["seq"]
                reverse_to_use = row["reverse"]
            else:
                new_id = partner_row["id"]
                sec_id = row["id"]
                seq_to_use = partner_row["seq"]
                reverse_to_use = partner_row["reverse"]

            how = cols[align_type]

            variants.append(f"{sec_id}:{how}:{cigar}")
            partnered.append([new_id, seq_to_use, reverse_to_use, len(seq_to_use), variants])

            df2 = df2.drop(index=partner, axis=0).reset_index().drop(columns="index")

    partnered = pd.DataFrame(partnered, columns=["id", "seq", "reverse", "sqlen", "variants"])
    # print(partnered.shape)
    # remove unpartnered -- they are surely not conserved between replicates

    return partnered


dfs = []
for peaksfile in called_peaks:
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
    dfs.append(df)

# map highly similar peaks on one another
# tree structure --> logarithmic number of mappings (all to all not needed)

to_merge = dfs

lvl = 0
while len(to_merge) > 1:
    to_merge = [merge_pair(to_merge[i:i + 2], i // 2, lvl) for i in range(len(to_merge)) if i % 2 == 0]
    lvl += 1
replicated = to_merge[0]

print("writing to output file...")
with open(output_file, mode='w') as writer:
    for i, row in replicated.iterrows():
        id = row["id"] + ";".join(row["variants"])
        seq = row["seq"]
        print(f">{id}\n{seq}", file=writer)