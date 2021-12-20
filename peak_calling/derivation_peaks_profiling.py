import numpy as np
import pandas as pd
import os
from scipy.signal import find_peaks
import argparse
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib.patches as patches

parser = argparse.ArgumentParser()
parser.add_argument("--aligned",
                    default=None,
                    type=str,
                    help="Path to the aligned directory. If directory does not exists, error will be thrown. Required."
                    )
parser.add_argument("--overview",
                    default=None,
                    type=str,
                    help="Path to the description csv. Required. Pairs with <--aligned> directory."
                    )
parser.add_argument("--k",
                    default=-1,
                    type=int,
                    help="Size of the k-mer created by BCALM. Required."
                    )
parser.add_argument("--output_filename",
                    default=None,
                    type=str,
                    help="Path to the output file. Required."
                    )
parser.add_argument("--rolling_window",
                    default=100,
                    type=int,
                    help="Maximal window size for rolling mean calculation. Optional."
                    )
parser.add_argument("--rolling_window_number",
                    default=20,
                    type=int,
                    help="Minimal number of windows in rolling mean calculation. Optional."
                    )
parser.add_argument("--prominence",
                    default=20,
                    type=int,
                    help="Minimal peak prominence. Optional."
                    )
parser.add_argument("--peak_max_width",
                    default=100,
                    type=int,
                    help="Maximal width of the identified peak. Optional."
                    )
parser.add_argument("--peak_min_width",
                    default=-1,
                    type=int,
                    help="Minimal width of the identified peak. Optional."
                    )
parser.add_argument("--weighting",
                    default="consensus",
                    type=str,
                    help="Select mode of peak sequence weighting. Possibilities: \n "
                         "\t- sq_count : for every occurence of a sequence in peak, write one. Its count is stored in "
                         "the FASTA header.\n "
                         "\t- consensus : write only a consensus sequence for every peak\n"
                    )
parser.add_argument("--draw",
                    action="store_true",
                    dest='draw',
                    help="Draw the profiles."
                   )

args = parser.parse_args([] if "__file__" not in globals() else None)

directory = args.aligned
k = args.k
outfile = args.output_filename

max_d = args.rolling_window
min_d_count = args.rolling_window_number
prominence = args.prominence
window = args.peak_max_width // 2


def get_peak_ranges(smooth_profile, prominence=5, window=75):
    peaks, props = find_peaks(smooth_profile, prominence=(prominence, np.inf), distance=1, width=1)
    neg_peaks, _ = find_peaks(-smooth_profile, prominence=(prominence, np.inf), distance=1, width=1)

    next_neg_peak = 0

    ranges = []

    if len(peaks) == 0:
        return ranges

    if len(neg_peaks) == 0:
        return [[np.fmax(peak - window, 0), np.fmin(len(profile), peak + window)] for peak in peaks
                if np.fmin(len(profile), peak + window) - np.fmax(peak - window, 0) > 1.5 * window]

    for peak in peaks:
        if np.fmin(len(profile), peak + window) - np.fmax(peak - window, 0) <= args.peak_min_width:
            continue
        prim_range = [np.fmax(peak - window, 0), np.fmin(len(profile), peak + window)]

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
    skipnext = False
    for i in range(len(ranges) - 1):
        if skipnext:
            skipnext = False
            continue

        a = ranges[i]
        b = ranges[i + 1]

        if a[1] >= b[0]:
            merged_ranges.append([a[0], b[1]])
            skipnext = True
        else:
            merged_ranges.append(a)
            if i + 1 == len(ranges) - 1:
                merged_ranges.append(b)
    return merged_ranges


encoded_bases = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "AG": "R",
    "CT": "Y",
    "GT": "K",
    "AC": "M",
    "CG": "S",
    "AT": "W",
    "CGT": "B",
    "AGT": "D",
    "ACT": "H",
    "ACG": "V",
    "ACGT": "N"
}


# bases = set("ACGT")


def encode(base_set):
    base_set = base_set.intersection(bases)
    return encoded_bases["".join(sorted(list(base_set)))]


bases = dict(
    A=0, C=1, G=2, T=3
)
bases['-'] = 4
rev_bases = {v: k for k, v in bases.items()}


def get_consensus(orig_profile, ranges, ident, min_size=-np.inf, max_size=np.inf):
    for from_, to in ranges:
        current = orig_profile[from_:to, :]
        sq = "".join([rev_bases[x] for x in np.argmax(current, axis=1) if rev_bases[x] != '-']).strip(' ')

        if (len(sq) > min_size) and (len(sq) < max_size):
            yield sq, f"ALIGNMENT {ident} CONSENSUS PEAK: from={from_}, to={to}"


def get_sqs_counts(sqs_df, ranges, ident, count_agreg_func=np.min, min_size=-np.inf, max_size=np.inf):
    items = sqs_df['alignment_actual'].str.split(',')  # .apply(lambda arr: np.array([float(x) for x in arr]))
    alignments = np.vstack(items).astype(float)

    i = 0
    for _, row in sqs_df[['sq', 'str_count']].iterrows():
        sq = row.sq
        str_count = row.str_count

        index_array = alignments[i, :] == 0
        aligned_sq = []
        aligned_count = []
        j = 0
        counts = [int(x) for x in str_count.split(',')]
        for val in index_array:
            if val:
                aligned_sq.append(sq[j])
                aligned_count.append(counts[j])
                j += 1
                if j >= len(sq):
                    break
            else:
                aligned_sq.append('-')
                aligned_count.append(None)

        for k, range in enumerate(ranges):
            from_,to = range
            peak_sq = ''.join([x for x in aligned_sq[from_:to] if x != '-'])
            counts_arr = [x for x in aligned_count[from_:to] if x is not None]
            if len(counts_arr) == 0:
                continue
            count = count_agreg_func(counts_arr)
            peak_head = f"ALIGNMENT {ident} PEAK {k} SQ {i}/{alignments.shape[0]}: from={from_}, to={to}; COUNT={count}"

            if (len(peak_sq) > min_size) and (len(peak_sq) < max_size):
                yield peak_sq, peak_head

        i+=1


def draw_peaks_profile(profile, smooth_profile, ranges, output_file):
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
    for begin, end in ranges:
        rect = patches.Rectangle((begin, 0), end - begin, height,
                                 linewidth=1, edgecolor=None, facecolor='xkcd:pale lavender')
        ax.add_patch(rect)

    plt.savefig(output_file, dpi=300)
    plt.close(fig)


handle = open(outfile, mode='w')

from progress.bar import FillingCirclesBar as progresscounter

files = os.listdir(directory)
df = pd.read_csv(args.overview)

if args.draw:
    drawing_path = '.'.join(args.output_filename.split('.')[:-1])+"_plots"
    os.mkdir(drawing_path)
else:
    drawing_path = None

with progresscounter("Peak calculation: ", max=len(files) - 1) as counter:
    for file in files:
        path = os.path.join(directory, file)
        orig_profile = np.load(path).T  # in 5 options
        profile = orig_profile.sum(axis=1)
        L = profile.shape[0]
        N = np.fmin(max_d, np.ceil(L / min_d_count).astype(int))
        smooth_profile = np.convolve(profile, np.ones(N) / N, mode='valid')  # rolling mean with N-sized window
        ranges = get_peak_ranges(smooth_profile, prominence=prominence, window=window)

        if len(ranges) == 0:
            counter.next()  # we are not drawing alignments without a peak
            continue

        if args.draw:
            draw_peaks_profile(profile, smooth_profile, ranges, os.path.join(drawing_path, file+".png"))

        align_ident = int(file.split('.')[0])

        if args.weighting == 'consensus':
            for sq, head in get_consensus(orig_profile, ranges, align_ident, min_size=k):
                print(f">{head}\n{sq}", file=handle)
        elif args.weighting == 'sq_count':
            sqs_in_current = df[df['alignment'] == align_ident]
            for sq, head in get_sqs_counts(sqs_in_current, ranges, align_ident, min_size=k):
                print(f">{head}\n{sq}", file=handle)
        else:
            raise Exception("This weighting mode does not exist. See help.")

        counter.next()
