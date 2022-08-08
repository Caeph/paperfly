import time
import numpy as np
from hmmlearn import hmm
import pandas as pd
from multiprocessing import Pool
from scipy.signal import find_peaks
from matplotlib import pyplot as plt
from matplotlib import patches
import seaborn as sns
from itertools import groupby
import os
import warnings
import argparse
from scipy.stats import mannwhitneyu

parser = argparse.ArgumentParser()
parser.add_argument("--assembled_path", default=None, type=str,
                    help="Path to aligned profiles.")  # aligned_corrected.csv
parser.add_argument("--working_dir", default=None, type=str, help="Path to working directory")
# parser.add_argument("--input_description", default=None, type=str, help="")  # input -- "ENCSR413CVQ.csv"
parser.add_argument("--treatment_counts", default=None, type=str,
                    help="Treatment k-mers counts file name (csv, not path).")
parser.add_argument("--control_counts", default=None, type=str, help="Control k-mers counts file name (csv, not path).")
parser.add_argument("--problems_no", default=None, type=int, help="Number of replicates for Bonferroni correction.")

parser.add_argument("--k", default=None, type=int, help="K-mer size")
parser.add_argument("--draw", dest='draw', action='store_true',
                    help="Draw profiles with count and HMM-based separation. Optional.", )
parser.add_argument("--threads", default=1, type=int, help="")
parser.add_argument("--output_path", default=None, type=str,
                    help="Path to output file (will be overwritten if exists).")
parser.add_argument("--pvalue_threshold", default=0.1, type=float, help="P-value: significance threshold for a peak to "
                                                                        "be reported. Default: 0.1.")
parser.add_argument("--window", default=120, type=int, help="Default window size.")


pseudocount = 0.01
scaling_coef = 0.99
states_len_thr = 2
scale_controls = False
pval_rounding_thr = 10e-8
max_d = 100
min_d_count = 20



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


def get_peak_ranges(smooth_profile, prominence=5, window=75):
    peaks, props = find_peaks(smooth_profile, prominence=(prominence, np.inf), distance=1, width=1)
    neg_peaks, _ = find_peaks(-smooth_profile, prominence=(prominence, np.inf), distance=1, width=1)

    profile_length = len(smooth_profile)

    next_neg_peak = 0

    ranges = []

    if len(peaks) == 0:
        return ranges

    if len(neg_peaks) == 0:
        return [[np.fmax(peak - window, 0), np.fmin(profile_length, peak + window), prom] for peak, prom in
                zip(peaks, props["prominences"])
                if np.fmin(profile_length, peak + window) - np.fmax(peak - window, 0) >  1.5 * args.k]

    for peak, prom in zip(peaks, props["prominences"]):
        if np.fmin(profile_length, peak + window) - np.fmax(peak - window, 0) <= 1.5 * args.k:
            continue
        prim_range = [np.fmax(peak - window, 0), np.fmin(profile_length, peak + window), prom]

        if next_neg_peak >= len(neg_peaks):
            ranges.append(prim_range)
            continue

        if (neg_peaks[next_neg_peak] > prim_range[0]) & (neg_peaks[next_neg_peak] < prim_range[1]) & (
                neg_peaks[next_neg_peak] < peak):
            prim_range[0] = neg_peaks[next_neg_peak] + 1

        if neg_peaks[next_neg_peak] < peak:
            next_neg_peak += 1
            if next_neg_peak >= len(neg_peaks):
                ranges.append(prim_range)
                continue
        neg_peak = neg_peaks[next_neg_peak]

        if (neg_peak > prim_range[0]) & (neg_peak < prim_range[1]):
            prim_range[1] = neg_peak - 1

        ranges.append(prim_range)

    merged_ranges = []
    last_range = None
    for i,range in enumerate(ranges):
        if last_range is None:
            last_range = range
            continue
        a = last_range
        b = range

        if a[1] >= b[0]: # ranges overlap
            last_range = [a[0], b[1], max(a[2], b[2])]
        else: # ranges dont overlap
            merged_ranges.append(last_range)
            last_range = range

    if last_range is not None:
        merged_ranges.append(last_range)

    return merged_ranges


def draw_peaks_profile(profile, smooth_profile, ranges, output_file, window):
    if len(ranges) == 0:
        return

    fig, ax = plt.subplots(1, 1, figsize=(3.5, 2))
    plt.subplots_adjust(bottom=0.2)

    sns.lineplot(
        x=np.arange(0, profile.shape[0]) - max_d // 2 + window // 2,
        y=profile,
        ax=ax,
        color='xkcd:cyan',
        label="original"
    )
    ax.set_xlabel("position in alignment", fontsize=8)
    ax.set_ylabel("enrichment", fontsize=8)

    sns.lineplot(
        x=np.arange(0, smooth_profile.shape[0]),
        y=smooth_profile,
        ax=ax,
        label="smoothed",
        color='xkcd:wine',
    )

    plt.legend(loc='upper left', fontsize=8)

    sns.despine(offset=1, trim=False)

    height = profile.max() + 10
    for begin, end, _ in ranges:
        rect = patches.Rectangle((begin, 0), end - begin, height,
                                 linewidth=1, edgecolor=None, facecolor='xkcd:pale lavender')
        ax.add_patch(rect)

    plt.savefig(output_file, dpi=300)
    plt.close(fig)


def main(args):
    for item in [args.assembled_path, args.working_dir, args.treatment_counts, args.control_counts,
                 args.problems_no, args.output_path]:
        if item is None:
            print("Insufficient input, exiting.")
            exit(1)

    # states = args.states
    # allowed_states = set([2,3,4,5])
    # if states not in allowed_states:
    #     print("Allowed state numbers are 2, 3, 4, 5.")
    #     exit(1)

    working_dir = args.working_dir
    assembled_sqs_file = args.assembled_path
    processes = args.threads
    draw = args.draw
    k = args.k
    pvalue = args.pvalue_threshold

    problems_no = args.problems_no
    treatment_filename = args.treatment_counts
    control_filename = args.control_counts

    drawing_path = os.path.join(working_dir, "profiles_pics")
    if args.draw:
        os.mkdir(drawing_path)

    # bonferroni correction:
    pvalue = pvalue / problems_no

    counts_df = None
    for file, label in zip([treatment_filename, control_filename], ["treatment", "control"]):
        counts = pd.read_csv(os.path.join(working_dir, file), sep=' ', header=None)
        counts.columns = ["seq", label + "_" + file]
        if counts_df is None:
            counts_df = counts
        else:
            counts_df = pd.merge(counts_df, counts, how="outer", on="seq")

    tr_to_ctrl = {treatment_filename : control_filename}

    assembled = pd.read_csv(assembled_sqs_file, sep=";", header=None)
    assembled.columns = ["assembled", "counts"]
    assembled["counts"] = assembled["counts"].apply(lambda x: [int(y) for y in x.split(",")])

    treatment_cols = ["treatment" + "_" + x for x in [treatment_filename]]
    control_cols = ["control" + "_" + x for x in [control_filename]]

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

    with open(args.output_path, mode='w') as peaks_writer:
        for profile_i, Y in enumerate(results):
            profile = Y[:,0]
            L = profile.shape[0]
            N = np.fmin(max_d, np.ceil(L / min_d_count).astype(int))
            smooth_profile = np.convolve(profile, np.ones(N) / N, mode='valid')

            prominence = 15

            peak_segments = np.zeros(len(Y))
            ranges = get_peak_ranges(smooth_profile, prominence=prominence, window=args.window)
            if args.draw:
                draw_peaks_profile(profile,
                                   Y[:, 1],  # control
                                   ranges,
                                   os.path.join(drawing_path, str(profile_i) + ".png"),
                                   args.window)

            for start, stop, _ in ranges:
                segment_pvals = []
                true_treatment = 0
                for ti, ci in [[0, 1]]:
                    control_vkt = Y[:, ci]
                    treatment_vkt = Y[:, ti]

                    y_c = control_vkt[start:stop]
                    y_t = treatment_vkt[start:stop]

                    # test whether treatment is more enriched than control
                    mwstat, mwpval = mannwhitneyu(y_t, y_c, alternative='greater')
                    if mwpval <= pvalue:
                        true_treatment += 1
                        segment_pvals.append(mwpval)

                if true_treatment >= len(treatment_cols) * 0.5:
                    # this is a peak
                    p_val_score = np.max(segment_pvals)
                    peak_segments[start:stop] = p_val_score

            peak_i = 1
            for start, end, _ in ranges:
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