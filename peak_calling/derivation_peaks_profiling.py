import numpy as np
import pandas as pd
import os
from scipy.signal import find_peaks
import argparse

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
                         "\t- count : for every occurence of a sequence in peak, write one.\n"
                         "\t- once : write every sequence once, regardless of the occurence count\n"
                         "\t- consensus : write only a consensus sequence for every peak\n"
                    )

args = parser.parse_args([] if "__file__" not in globals() else None)

directory = args.aligned
k = args.k
outfile = args.output_filename

max_d = args.rolling_window
min_d_count = args.rolling_window_number
prominence = args.prominence
window = args.peak_max_width // 2

if args.peak_min_width == -1:
    args.peak_min_width = 1.5 * window


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


def print_sequences(sqs_df, ranges, outfile_handle, prev_path, min_segment_len=k, N=12, mode=args.weighting):
    pass
    # if mode == "consensus":
    #     ranges_results = {}
    #     for line in lines_iter:
    #         count, sq = line.strip('\n').split('\t')
    #         segments = [sq[begin:end] for begin, end in ranges if len(sq[begin:end]) > min_segment_len]
    #         for i,segment in enumerate(segments):
    #             if i in ranges_results:
    #                 ranges_results[i].append(segment)
    #             else:
    #                 ranges_results[i] = [segment]
    #     i = 0
    #     for peak in ranges_results.values():
    #         length = max(len(x.strip(' ')) for x in peak)
    #         peak_sqs = set(x for x in peak if len(x.strip(' ')) == length)
    #         for peak_sq in peak_sqs:
    #             print(
    #                 f">{''.join(random.choices(string.ascii_uppercase + string.digits, k=N))}: {prev_path} {i}",
    #                 file=outfile_handle)
    #             print(f"{peak_sq}", file=outfile_handle)
    #         i+=1
    # else:
    #     peaks = Counter()
    #     for line in lines_iter:
    #         count, sq = line.strip('\n').split('\t')
    #         for begin, end in ranges:
    #             segment = sq[begin:end].strip(' ')
    #             if len(segment) > min_segment_len:
    #                 peaks[segment] += int(count)
    #     for segment, count in peaks.items():
    #         if mode == "count":
    #             for i in range(count):
    #                 print(f">{''.join(random.choices(string.ascii_uppercase + string.digits, k=N))}: {prev_path}, {i + 1}/{count}", file=outfile_handle)
    #                 print(f"{segment}", file=outfile_handle)
    #                 pass
    #         elif mode == "once":
    #             print(
    #                 f">{''.join(random.choices(string.ascii_uppercase + string.digits, k=N))}: {prev_path}, count={count}",
    #                 file=outfile_handle)
    #             print(f"{segment}", file=outfile_handle)
    #         else:
    #             raise Exception("Mode unrecognized.")


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


def get_sqs_counts(sqs_df, ranges):
    # todo
    pass


handle = open(outfile, mode='w')

from progress.bar import FillingCirclesBar as progresscounter

files = os.listdir(directory)
df = pd.read_csv(args.overview)

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
            counter.next()
            continue

        align_ident = int(file.split('.')[0])

        if args.weighting == 'consensus':
            for sq, head in get_consensus(orig_profile, ranges, align_ident, min_size=k):
                print(f">{head}\n{sq}", file=handle)
        elif args.weighting == 'sq_count':
            raise NotImplementedError("To do.")

            # sqs_in_current = df[df['alignment'] == align_ident]
            # for sq, head in get_sqs_counts(sqs_in_current, ranges):
            #     print(f">{head}\n{sq}", file=handle)
        else:
            raise Exception("This weighting mode does not exist. See help.")

        counter.next()
