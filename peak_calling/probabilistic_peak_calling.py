import time
import numpy as np
from hmmlearn import hmm
import pandas as pd
from multiprocessing import Pool
from matplotlib import pyplot as plt
from matplotlib import patches
import seaborn as sns
from itertools import groupby
import os
import warnings
import argparse
from scipy.stats import mannwhitneyu

parser = argparse.ArgumentParser()
parser.add_argument("--assembled_path", default=None, type=str, help="Path to aligned profiles.")  # aligned_corrected.csv
parser.add_argument("--working_dir", default=None, type=str, help="Path to working directory")
parser.add_argument("--input_description", default=None, type=str, help="")  # input -- "ENCSR413CVQ.csv"
parser.add_argument("--k", default=None, type=int, help="K-mer size")
parser.add_argument("--draw", dest='draw', action='store_true',
                    help="Draw profiles with count and HMM-based separation. Optional.", )
parser.add_argument("--threads", default=1, type=int, help="")
parser.add_argument("--output_path", default=None, type=str,
                    help="Path to output file (will be overwritten if exists).")
parser.add_argument("--pvalue_threshold", default=0.1, type=float, help="P-value: significance threshold for a peak to "
                                                                        "be reported. Default: 10%.")

states = 5

pseudocount = 0.01
scaling_coef = 0.99
states_len_thr = 2
scale_controls = True
pval_rounding_thr = 10e-8

absolute_minimal = 10

state_colors = {
    0: 'xkcd:pale lavender',
    1: 'xkcd:pale teal',
    2: 'xkcd:ecru',
    3: 'xkcd:carolina blue',
    4: 'xkcd:ice blue'
}


def get_positions_vectors(length, last, aggr_counts, k):
    sequence_df = aggr_counts[last:last + length, :]
    positions_length = length + k - 1
    vectors = []
    for i in range(positions_length):
        start = np.fmax(0, i - k + 1)
        P = sequence_df[start:i + 1, :]
        if len(P) > 0:
            P = np.max(P, axis=0)
            vectors.append(P)
    return np.vstack(vectors)


def main(args):
    for item in [args.assembled_path, args.working_dir, args.input_description, args.output_path]:
        if item is None:
            print("Insufficient input, exiting.")
            exit(1)

    input_description = args.input_description
    working_dir = args.working_dir
    assembled_sqs_file = args.assembled_path
    processes = args.threads
    draw = args.draw
    k = args.k
    pvalue = args.pvalue_threshold

    profiles_pictures = os.path.join(working_dir, "profiles_pics")

    df_desc = pd.read_csv(input_description, sep="\t")
    # for col in ["fastq", "control"]:
    #     df_desc[col] = df_desc[col].str.replace("fastq.gz", "fasta", regex=False).str.replace("fastq", "fasta", regex=False)

    tr_to_ctrl = {rdict["fastq"]: rdict["control"] for rdict in df_desc.to_dict(orient="records")}

    treatments = df_desc["fastq"].unique()
    control = df_desc["control"].unique()
    all_files = list(treatments) + list(control)

    # bonferroni correction:
    pvalue = pvalue / len(treatments)

    counts_df = None
    for file in all_files:
        if file in control:
            label = "control"
        else:
            label = "treatment"

        counts = pd.read_csv(os.path.join(working_dir, file + ".csv"), sep=' ', header=None)
        counts.columns = ["seq", label + "_" + file]
        if counts_df is None:
            counts_df = counts
        else:
            counts_df = pd.merge(counts_df, counts, how="outer", on="seq")

    assembled = pd.read_csv(assembled_sqs_file, sep=";", header=None)
    assembled.columns = ["assembled", "counts"]
    assembled["counts"] = assembled["counts"].apply(lambda x: [int(y) for y in x.split(",")])

    treatment_cols = ["treatment" + "_" + x for x in treatments]
    control_cols = ["control" + "_" + x for x in control]

    if scale_controls:
        scale = counts_df[treatment_cols].values.sum() / counts_df[control_cols].values.sum()
        scale *= scaling_coef
        for col in control_cols:
            counts_df[col] = counts_df[col] * scale

    counts_df["reconstruction_count"] = np.fmax(
        counts_df[treatment_cols].sum(axis=1) - counts_df[control_cols].sum(axis=1), 0)

    assembled["kmer_series"] = assembled["assembled"].apply(lambda sq: [sq[i:i + k] for i in range(len(sq) - k + 1)])
    assembled["kmer_count_series"] = assembled["counts"].apply(
        lambda sq: ([min(sq[i:i + k]) for i in range(len(sq) - k + 1)]))

    kmers = np.hstack(assembled["kmer_series"].values)
    kmer_counts = np.hstack(assembled["kmer_count_series"].values)
    kmers_df = pd.DataFrame(np.vstack([kmers, kmer_counts]).T, columns=["kmer", "used_count"])

    seq_kmer_lengths = assembled["kmer_series"].str.len()

    counts_kmer_df = pd.merge(kmers_df, counts_df, how="left", left_on="kmer", right_on="seq").drop(columns=["seq"])
    counts_kmer_df["used_count"] = counts_kmer_df["used_count"].astype(int)

    counts_kmer_df["frac"] = counts_kmer_df["used_count"] / counts_kmer_df["reconstruction_count"] + pseudocount
    counts_kmer_df = counts_kmer_df.fillna(0)

    for col in [*treatment_cols, *control_cols]:
        counts_kmer_df[col + "_frac"] = counts_kmer_df[col] * counts_kmer_df["frac"]

    count_cols = [col + "_frac" for col in [*treatment_cols, *control_cols]]

    first_indices = seq_kmer_lengths.cumsum()
    first_indices = np.hstack([[0], first_indices.values[:-1]])

    aggr_counts = counts_kmer_df[count_cols].values

    start = time.time()
    with Pool(processes) as pool:
        results = pool.starmap(get_positions_vectors, iterable=zip(seq_kmer_lengths,
                                                                   first_indices,
                                                                   [aggr_counts for _ in range(len(first_indices))],
                                                                   [k for _ in range(len(first_indices))]
                                                                   )
                               )
    print(f"Position counts calculated, elapsed: {time.time() - start}")

    del first_indices, seq_kmer_lengths, counts_kmer_df, kmers_df, counts_df
    print("Prediction step initiated...")

    def draw_peaks_profile(profile, states, output_file=None, show=False):
        fig, ax = plt.subplots(1, 1, figsize=(15, 10))
        plt.subplots_adjust(bottom=0.2)

        for x in range(profile.shape[1]):
            sns.lineplot(
                x=np.arange(0, profile.shape[0]),
                y=profile[:, x],
                ax=ax,
                # color='xkcd:cyan',
                linewidth=5.0,
                label=count_cols[x]
            )
        ax.set_xlabel("position in alignment", fontsize=8)
        ax.set_ylabel("enrichment", fontsize=8)

        plt.legend(loc='upper left', fontsize=8)

        sns.despine(offset=1, trim=False)

        height = profile.max() + 5
        min = np.fmin(0, profile.min() - 5)
        for i, state in enumerate(states):
            rect = patches.Rectangle((i, min), 1, height,
                                     linewidth=1, edgecolor=None, facecolor=state_colors[state])
            ax.add_patch(rect)

        if output_file is not None:
            plt.savefig(output_file, dpi=300)
        if show:
            plt.show()
        plt.close(fig)

    def smooth_states(states, too_short_thr=states_len_thr):
        blocks = np.array([len(list(group)) for st, group in groupby(states)])
        blocks_states = np.array([st for st, group in groupby(states)])
        starts = np.hstack([[0], blocks.cumsum()])[:-1]

        candidates = np.where(blocks <= too_short_thr)[0]
        for cand in candidates:
            if (cand == 0) or (cand == len(blocks) - 1):
                if cand == 0:
                    newstate = blocks_states[1]
                else:
                    newstate = blocks_states[-2]
                states[starts[cand]:starts[cand] + blocks[cand]] = newstate

            elif blocks_states[cand - 1] == blocks_states[cand + 1]:
                start = starts[cand]
                l = blocks[cand]
                states[start:start + l] = blocks_states[cand + 1]

        return states

    def get_state_ranges(states):
        block_lengths = np.array([len(list(group)) for st, group in groupby(states)])
        # blocks_states = np.array([st for st, group in groupby(states)])
        starts = np.hstack([[0], block_lengths.cumsum()])[:-1]

        return np.vstack([starts, starts+block_lengths]).T

    if draw:
        os.makedirs(profiles_pictures, exist_ok=True)

    def get_segments(Y, states):
        if len(Y) < 1.75 * k:
            smooth_predictions = np.ones(len(Y))
            return smooth_predictions
        else:
            # seek segments
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
                    peak_identifier = hmm.GaussianHMM(n_components=states, n_iter=100)
                    peak_identifier.fit(Y)

                    # gives me segments
                    log_odds, predictions = peak_identifier.decode(Y)
                    # correction of outlier in majorities
                    smooth_predictions = smooth_states(predictions)
                    # compare control and treatment
                except ValueError:
                    # catch ValueError: transmat_ rows must sum to 1 (got [1. 0.])
                    # print(f"Skipping profile {profile_i}, value error occured")
                    return None

        return smooth_predictions

    with open(args.output_path, mode='w') as peaks_writer:
        for profile_i, Y in enumerate(results):
            # X = Y
            smooth_predictions = get_segments(Y, states)
            if smooth_predictions is None:
                trial = get_segments(Y, states-1)
                if trial is None:
                    smooth_predictions = np.zeros(len(Y))
                else:
                    smooth_predictions = trial

            # calculate peaks
            peak_segments = np.zeros(len(Y))
            tc_index_pairs = []
            for tr_col in treatment_cols:
                ctrl = tr_to_ctrl[tr_col.replace("treatment_", "")]
                ctrl_col = f"control_{ctrl}_frac"
                ci = count_cols.index(ctrl_col)
                ti = count_cols.index(f"{tr_col}_frac")
                tc_index_pairs.append((ti, ci))

            for start, stop in get_state_ranges(smooth_states(smooth_predictions)):
                segment_pvals = []
                true_treatment = 0
                for ti, ci in tc_index_pairs:
                    control_vkt = Y[:, ci]
                    treatment_vkt = Y[:, ti]

                    y_c = control_vkt[start:stop]
                    y_t = treatment_vkt[start:stop]

                    # test whether treatment is enriched sufficiently
                    if np.median(y_t) < absolute_minimal:
                        continue

                    # test whether treatment is more enriched than control
                    mwstat, mwpval = mannwhitneyu(y_t, y_c, alternative='greater')
                    if mwpval <= pvalue:
                        true_treatment += 1
                        segment_pvals.append(mwpval)

                if true_treatment >= len(treatment_cols) * 0.5:
                    # this is a peak
                    p_val_score = np.max(segment_pvals)
                    peak_segments[start:stop] = p_val_score

            peak_segments_states = np.zeros(len(Y))
            peak_segments_states[np.where(peak_segments)[0]] = 1
            if draw:
                draw_peaks_profile(Y,
                                   peak_segments_states,
                                   output_file=os.path.join(profiles_pictures, f"{profile_i}.png"),
                                   show=False)

            # print peak sequence to file
            peak_ranges = [(int(start), int(end)) for start, end in get_state_ranges(peak_segments_states) if
                           peak_segments_states[start] == 1]

            peak_i = 1
            for start, end in peak_ranges:
                if (end - start) < k:
                    # too short
                    continue
                # print(start,end)
                peak_sq = assembled.iloc[profile_i]["assembled"][start:end]
                pvals = np.round(np.max(peak_segments[start:end]), decimals=7)
                if pvals < pval_rounding_thr:
                    pvals = f"<{pval_rounding_thr}"
                else:
                    pvals = np.round(np.max(peak_segments[start:end]), decimals=7)
                print(f">assembled_sq_{profile_i}__peak_{peak_i}__pvalues_{pvals}\n{peak_sq}", file=peaks_writer)
                peak_i += 1


if __name__ == '__main__':
    args = parser.parse_args([] if "__file__" not in globals() else None)
    main(args)