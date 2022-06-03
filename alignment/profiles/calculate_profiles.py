import os
import re
import time
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
import math
import edlib
from progress.bar import IncrementalBar as Bar
from multiprocessing import Pool
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--pools",
                    default=4,
                    type=int,
                    help="Number of threads to use in aligning. Default 4. Optional."
                    )
parser.add_argument("--misses",
                    default=5,
                    type=float,
                    help="Number of allowed substitutions/insertions/deletions in aligning a sequence of length k. "
                         "For longer sequences, this is scaled. "
                    )
parser.add_argument("--aligned",
                    default=None,
                    type=str,
                    help="Path to the output aligned directory. Required."
                    )
parser.add_argument("--overview",
                    default=None,
                    type=str,
                    help="Path to the output description csv. Required. Pairs with <--aligned> directory."
                    )
parser.add_argument("--k",
                    default=-1,
                    type=int,
                    help="Size of the k-mer created by BCALM. Required."
                    )
parser.add_argument("--input",
                    default=None,
                    type=str,
                    help="Path to the input file."
                    )
parser.set_defaults(all_sqs_result=False)

args = parser.parse_args([] if "__file__" not in globals() else None)

bases = dict(A=0, C=1, G=2, T=3)
bases['-'] = 4
rev_bases = {v: k for k, v in bases.items()}
global_alignment_ident_no = 0


operations = {
    '.' : 0,
    '-' : 1,
    '|' : 0
}


class AlignmentProfile:
    def __init__(self, width, df, identifier):
        self.ident = identifier

        self.profile = np.zeros((5, width))
        self.repre_sq = ""
        self.seq_alignments = None  # this will be a pandas df
        self.seq_align_counter = -1

        self.calculate_profile(df)

    def calculate_profile(self, df):
        self.seq_alignments = pd.DataFrame([(index, *np.zeros(self.profile.shape[1], dtype=np.int8)) for index in df.index])

        unwrapped_sq = df['sq'].str.split('', expand=True)
        unwrapped_sq = unwrapped_sq.drop(columns=[unwrapped_sq.columns[0], unwrapped_sq.columns[-1]])

        counts = np.stack(df['count'].values)

        for base in bases:
            a = unwrapped_sq != base
            newX = np.ma.array(counts, mask=a)
            new_counts = newX.sum(axis=0)
            self.profile[bases[base], :] += new_counts

        # repre_sq
        maxs = np.argmax(self.profile, axis=0)
        self.repre_sq = "".join([rev_bases[x] for x in maxs])

    def add_sequence(self, new_sq, new_counts, nice, sq_index):
        offset = re.search(nice['target_aligned'].replace('-', ''), self.repre_sq).start(0)
        x = self.profile
        # padding with the following number of observed positions (sum of all bases)

        # pad profile with insertions
        insertions = np.where(np.array(list(nice['target_aligned'])) == '-')[0]
        for i, index in enumerate(insertions):
            if x.shape[1] >= index:
                value = 0
            else:
                value = x[:, index].sum()
            x = np.insert(x, index + offset, [0, 0, 0, 0, value], axis=1)
            self.seq_alignments.insert(loc=int(index+offset), column=self.seq_align_counter, value=1)
            self.seq_align_counter -= 1

        # pad new counts with deletions
        aligned_query = np.array(list(nice['query_aligned']))
        deletions = np.where(aligned_query == '-')[0]
        for i, index in enumerate(deletions):
            value = new_counts[index]
            new_counts = np.insert(new_counts, index, value, axis=0)

        i = offset
        for base, count in zip(aligned_query, new_counts):
            x[bases[base], i] += count
            i += 1

        self.profile = x

        # store new sequence alignment
        added_alignment = -np.ones(self.profile.shape[1])
        for i, char in enumerate(nice['target_aligned']):
            if char == '-':
                added_alignment[offset + i] = 1
            else:
                added_alignment[offset + i] = 0
        self.seq_alignments.loc[-1] = [sq_index, *added_alignment]  # adding a row
        self.seq_alignments.index = self.seq_alignments.index + 1  # shifting index

        # recalculate repre_sq -- the most probable one
        maxs = np.argmax(self.profile, axis=0)
        self.repre_sq = "".join([rev_bases[x] for x in maxs if rev_bases[x] != '-'])  # '-' is removed from the sq


def dst_func(x, y):
    return (np.array(x) != np.array(y)).sum()


def read_alignment(filename):
    for line in open(filename):
        sq, count = line.strip('\n').split(';')
        yield sq, np.array([int(x) for x in count.split(',')]), count


def cluster_group(df_group, l, dst=dst_func):
    sqs = df_group.reset_index()['sq']
    n = len(sqs)

    if n <= 1:
        return np.zeros(n)

    dst_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i):
            d = dst(sqs[i], sqs[j])
            dst_matrix[i, j] = d
            dst_matrix[j, i] = d

    model = AgglomerativeClustering(distance_threshold=threshold * l,
                                    n_clusters=None,
                                    linkage='complete',
                                    affinity='precomputed')
    clusters = model.fit_predict(dst_matrix)
    return clusters


aligned_sqs_file = args.input
k = args.k
misses = args.misses
pools = args.pools

threshold = misses / k
if args.aligned is None:
    output_profile_dir = aligned_sqs_file + "_profiles"
else:
    output_profile_dir = args.aligned

if args.overview is None:
    output_csv_file = aligned_sqs_file + "_overview.csv"
else:
    output_csv_file = args.overview

# read
df = pd.DataFrame(read_alignment(aligned_sqs_file))
df.columns = ['sq', 'count', 'str_count']
df['length'] = df['sq'].str.len()
# df['alignment'] = -1  # every aligned sq has an alignment identification
groups = df.groupby(by='length')

unique_lengths = df['length'].sort_values(ascending=False).unique()

against = []

longest = unique_lengths[0]
df_group = groups.get_group(longest).copy()

clusters = cluster_group(df_group, longest)
df_group['cluster'] = clusters

alignments = {
}

for cluster, cluster_df in df_group.groupby(by='cluster'):
    alignment = AlignmentProfile(longest, cluster_df, global_alignment_ident_no)
    alignments[global_alignment_ident_no] = alignment

    global_alignment_ident_no += 1
    against.append(alignment)

    # df.loc[df['sq'].isin(cluster_df['sq']), 'alignment'] = alignment.ident

    # to each sequence


start = time.time()

# print(df.groupby(by='length').get_group(longest))
# print("running on shorter")

with Bar("Processing length groups...", max=len(unique_lengths) - 1) as bar:
    for length in unique_lengths[1:]:
        bar.next()
        df_group = groups.get_group(length).copy()

        def getDistanceAndAlignment(sq):
            # this is a fallback, it should not happen
            maxval = np.floor(threshold * len(sq))

            min = np.inf
            min_target = None

            if maxval < 1:
                return min,min_target

            for target in against:
                align_res = edlib.align(sq, target.repre_sq, mode='HW', task='distance', k=maxval)
                if align_res['editDistance'] != -1:
                    if min > align_res['editDistance']:
                        if align_res['editDistance'] == 0:
                            return align_res['editDistance'], target.ident

                        min = align_res['editDistance']
                        min_target = target

            if min_target is not None:
                min_target = min_target.ident

            return min, min_target

        x = length * threshold
        if length * threshold >= 1:
            # try align
            with Pool(pools) as pool:
                result = pool.map(getDistanceAndAlignment, df_group['sq'])
            df_group['aligned'] = result

            # add aligned to profiles
            aligned = df_group[df_group['aligned'] != (np.inf, None)]
            for index, row in aligned.iterrows():
                to = alignments[row['aligned'][1]]
                align_res = edlib.align(row.sq, to.repre_sq, mode='HW', task='path')
                nice = edlib.getNiceAlignment(align_res, row.sq, to.repre_sq)
                to.add_sequence(row.sq, row['count'], nice, index)
                # df.loc[df['sq'] == row.sq, 'alignment'] = to.ident

            # cluster unaligned, add to against
            unaligned = df_group[df_group['aligned'] == (np.inf, None)].copy()
            clusters = cluster_group(unaligned, length)
            unaligned['cluster'] = clusters

            for cluster, cluster_df in unaligned.groupby(by='cluster'):
                alignment = AlignmentProfile(length, cluster_df, global_alignment_ident_no)
                alignments[global_alignment_ident_no] = alignment
                global_alignment_ident_no += 1
                against.append(alignment)
        else:
            # threshold is less than one, no clustering nor alignment takes place
            df_group["aligned"] = [(np.inf, None) for _ in range(len(df_group))]
            unaligned = df_group.copy()
            unaligned["cluster"] = list(range(len(unaligned)))
            # print(f"pseudoclustering elapsed: {time.time() - s}")

            s = time.time()
            for i, row in unaligned.iterrows():
                cluster_df = pd.DataFrame(row).T
                alignment = AlignmentProfile(length, cluster_df, global_alignment_ident_no)
                alignments[global_alignment_ident_no] = alignment
                global_alignment_ident_no += 1
                against.append(alignment)
            # print(f"alignment elapsed: {time.time() - s}")


print(f"{aligned_sqs_file} elapsed: {time.time() - start}")
print(f"{aligned_sqs_file} writing...")


os.makedirs(output_profile_dir, exist_ok=True)
for alignment in against:
    filename = f"{output_profile_dir}/{alignment.ident}.prf"
    np.save(filename, alignment.profile)

# get actual alignment for each sq
all_alignments = []
for alignment in against:
    itemized = alignment.seq_alignments
    num_cols = itemized.columns[1:]
    # index_col = itemized.columns[0]
    # translate to sth readable
    for col in num_cols:
        itemized[col] = itemized[col].astype(int).apply(str)

    itemized['alignment_actual'] = itemized[num_cols].agg(','.join, axis=1)  # todo maybe cigar?
    itemized = itemized.drop(columns=num_cols)
    itemized.columns = ['index_df', 'alignment_actual']
    itemized['alignment'] = alignment.ident
    all_alignments.append(itemized)

all_alignments = pd.concat(all_alignments)
merged = pd.merge(all_alignments, df, left_on='index_df', right_index=True)


# write sequences in df
merged.drop(columns=['count', 'index_df']).to_csv(output_csv_file, index=False)
print(f"{aligned_sqs_file} done")
