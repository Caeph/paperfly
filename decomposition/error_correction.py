import os
import time

from progress.bar import ChargingBar as bar
import networkx as nx
import preprocessing as prep
import edlib
from itertools import groupby
import argparse
from multiprocessing import Pool
import heapq

# this is done as described in Velvet genome assembler -- Tour bus algorithm

parser = argparse.ArgumentParser()
parser.add_argument("--raw_components",
                    default=None,
                    type=str,
                    help="Path to input: directory of fastas of a weakly connected component."
                    )
parser.add_argument("--corrected_components",
                    default=None,
                    type=str,
                    help="Path to output: directory of fastas of a corrected weakly connected component."
                    )
parser.add_argument("--k",
                    default=31,
                    type=int,
                    help="Size of the k-mer created by BCALM. Defaults at 31."
                    )
parser.add_argument("--replaced",
                    default=None,
                    type=str,
                    help="Path to a tsv with information on errors"
                    )
parser.add_argument("--max_bubble_sequence",
                    default=None,
                    type=int,
                    help="Maximal length of a sequence in bubble to be merged."
                    )


def archive_ungapped(sq, mapped_to, alignment_score, counts, contracted_info):
    counts = ",".join([str(x) for x in counts])
    print(f"count included\t{sq}\t{mapped_to}\t{counts}\t{alignment_score}", file=contracted_info)
    # count included in mapped


def archive_gapped(sq, mapped_to, alignment_score, counts, contracted_info):
    counts = ",".join([str(x) for x in counts])
    print(f"count excluded\t{sq}\t{mapped_to}\t{counts}\t{alignment_score}", file=contracted_info)


def sequence_from_path(G, k, path):
    if len(path) == 0:
        return "", None
    endings = {}
    sqs = [G.nodes[path[0]]["sq"]]
    current_size = len(G.nodes[path[0]]["sq"])
    endings[path[0]] = current_size  # where the part first seen in the node ends in the alignment
    for node in path[1:]:
        additional = G.nodes[node]["sq"][k - 1:]
        current_size += len(additional)
        endings[node] = current_size
        sqs.append(
            additional
        )
    return "".join(sqs), endings


def compare_sqs(s1, s2):
    align = edlib.align(s1, s2, mode='NW', task='path')

    similarity = 1 - (2 * align["editDistance"] / (len(s1) + len(s2)))  # frakce stejnych
    return similarity, align


def get_endings_in_alignment(s1, s2, e1, e2, align):
    nice = edlib.getNiceAlignment(align, s1, s2)
    alignment = {
        "s1": nice["query_aligned"],
        "s2": nice["target_aligned"],
        "score": align["editDistance"]
    }

    endings_in_align = {}  # vertex: end index
    for label, sq, endings in zip(["s1", "s2"], [s1, s2], [e1, e2]):
        aligned = alignment[label]
        # find indexes of gaps
        gap_indices = [i for i, x in enumerate(aligned) if x == '-']
        for index in gap_indices:
            for k in endings:
                if endings[k] >= index:
                    endings[k] += 1

        for v in endings:
            if v not in endings_in_align:
                endings_in_align[v] = endings[v]
    return alignment, endings_in_align


def get_new_node_info(aligned, endindex, last_endindex,
                      more_abundant, counts1, counts2, k, contracted_info):
    if last_endindex is None:
        init = 0
    else:
        init = last_endindex - k + 1

    represented1 = aligned["s1"][init:endindex]
    represented2 = aligned["s2"][init:endindex]
    kmers = endindex - init - k + 1

    # no gaps expected
    c1, c2 = counts1[init:init + kmers], counts2[init:init + kmers]

    if represented1 == represented2:
        return represented1, c1

    new_counts = [x + y for x, y in zip(c1, c2)]
    new_sq = [represented1, represented2][more_abundant]

    if more_abundant == 0:
        archive_ungapped(represented2,
                         mapped_to=represented1,
                         alignment_score=aligned["score"],
                         counts=c2,
                         contracted_info=contracted_info
                         )
    else:
        archive_ungapped(represented1,
                         mapped_to=represented2,
                         alignment_score=aligned["score"],
                         counts=c1,
                         contracted_info=contracted_info
                         )

    return new_sq, new_counts


def get_bubble_startpoint(subG, bubble_endpoint):
    # endpoint has always two incoming edges

    e1, e2 = subG.in_edges(bubble_endpoint)

    # probe first
    seen = set()
    walkthrough1 = []
    queue = [e1[0]]
    while len(queue) > 0:
        current = queue.pop(0)
        walkthrough1.append(current)
        seen.add(current)

        # if subG.in_degree(current) > 1:
        #     print(current)

        for prev_, _ in subG.in_edges(current):
            if prev_ in seen:
                continue
            # stack.insert(0, prev_)
            queue.append(prev_)

    # probe second
    seen = set()
    queue = [e2[0]]
    while len(queue) > 0:
        current = queue.pop(0)
        if current in walkthrough1:
            return current
        seen.add(current)
        for prev_, _ in subG.in_edges(current):
            if prev_ in seen:
                continue
            # stack.insert(0, prev_)
            queue.append(prev_)

    return None


def get_more_abundant(G, path1, path2):
    cs1, l1 = 0, 0
    for item in path1:
        cs1 += sum([x for x in G.nodes[item]["counts"]])
        l1 += len(G.nodes[item]["counts"])

    if l1 == 0:
        count_s1 = 0
    else:
        count_s1 = cs1 / l1

    cs2, l2 = 0, 0
    for item in path2:
        cs2 += sum([x for x in G.nodes[item]["counts"]])
        l2 += len(G.nodes[item]["counts"])

    if l2 == 0:
        count_s2 = 0
    else:
        count_s2 = cs2 / l2

    if count_s1 >= count_s2:
        return 0
    return 1


def contract_bubble(G, k, bubble_endpoint, seen_vertices, max_bubble_length, contracted_info, seen_edges):
    seen_vertices.add(bubble_endpoint)

    subG = nx.DiGraph()
    subG.add_edges_from(seen_edges)

    if len(subG.in_edges(bubble_endpoint)) != 2:  # a cycle is involved
        return False

    bubble_startpoint = get_bubble_startpoint(subG, bubble_endpoint)

    # test for cycles
    if bubble_startpoint == bubble_endpoint:
        return False

    paths_to_contract = list(nx.all_simple_paths(subG, bubble_startpoint, bubble_endpoint))
    if len(paths_to_contract) < 2:
        # this was a cycle
        return False

    if len(paths_to_contract) != 2:
        print("ERROR INFO:")
        print(paths_to_contract)
        print("subG:")
        print(subG.edges())
        raise Exception("Wrong number of paths to contract.")

    # generate sequences, compare them
    path1, path2 = paths_to_contract[0], paths_to_contract[1]

    sq1, endings1 = sequence_from_path(G, k, path1)
    sq2, endings2 = sequence_from_path(G, k, path2)

    if (len(sq1) > max_bubble_length) and (len(sq2) > max_bubble_length):
        return None

    counts1, counts2 = [], []
    for v in path1:
        counts1.extend(G.nodes[v]["counts"])
    for v in path2:
        counts2.extend(G.nodes[v]["counts"])

    # if one of the paths is gapped, keep the whole one and archive the other
    # if both are, keep the more abundant and keep it, archive the other
    similarity, alignment = compare_sqs(sq1, sq2)
    if ("I" in alignment["cigar"]) and ("D" in alignment["cigar"]):
        # both are gapped: pick the more abundant
        more_abundant = get_more_abundant(G, path1, path2)
        vertices_to_remove = set([path1, path2][more_abundant])
        archived_sq = [sq1, sq2][more_abundant]
        archived_counts = [counts1, counts2][more_abundant]
        if more_abundant == 0:
            archive_gapped(archived_sq,
                           mapped_to=sq2,
                           alignment_score=alignment["editDistance"],
                           counts=archived_counts,
                           contracted_info=contracted_info
                           )
        else:
            archive_gapped(archived_sq,
                           mapped_to=sq1,
                           alignment_score=alignment["editDistance"],
                           counts=archived_counts,
                           contracted_info=contracted_info
                           )
        unchanged = [x for x in G.nodes if x not in vertices_to_remove]
        newG = G.subgraph(unchanged).copy()

        return newG
    elif "D" in alignment["cigar"]:
        # path1 is gapped
        vertices_to_remove = path1
        archive_gapped(sq1,
                       mapped_to=sq2,
                       alignment_score=alignment["editDistance"],
                       counts=counts1,
                       contracted_info=contracted_info
                       )
        unchanged = [x for x in G.nodes if x not in vertices_to_remove]
        newG = G.subgraph(unchanged).copy()

        return newG
    elif "I" in alignment["cigar"]:
        # path2 is gapped
        vertices_to_remove = path2
        archive_gapped(sq2,
                       mapped_to=sq1,
                       alignment_score=alignment["editDistance"],
                       counts=counts2,
                       contracted_info=contracted_info
                       )
        unchanged = [x for x in G.nodes if x not in vertices_to_remove]
        newG = G.subgraph(unchanged).copy()

        return newG

    # generate new graph with contracted vertices
    changed = set(path1).union(set(path2))
    unchanged = [x for x in G.nodes if x not in changed]
    newG = G.subgraph(unchanged).copy()

    # managing ungapped alignments
    aligned_sqs, endings_aligned = get_endings_in_alignment(sq1, sq2, endings1, endings2, alignment)

    reversed_endings = {}
    for vertex, index in endings_aligned.items():
        if index in reversed_endings:
            reversed_endings[index].append(vertex)
        else:
            reversed_endings[index] = [vertex]

    # find vertices, after which there must be a break in alignment]
    branchings_after = []
    for path in [path1, path2]:
        # find vertices, after which there is a branch out of the bubble
        for n in path[:-1]:  # not interested in bubble endpoint
            if G.out_degree(n) > 1:  # makes unnecessary break after bubble startpoint but whatever
                branchings_after.append(n)
        # find vertices, before which there is a branch into the bubble
        for i, n in enumerate(path[1:]):  # not interested in bubble startpoint
            if G.in_degree(n) > subG.in_degree(n):
                branchings_after.append(path[i])  # append preceding vertex

    # add changed to new graph
    branchings_after = list(set(branchings_after))
    breaks = [(endings_aligned[v], v) for v in branchings_after]
    breaks.append((len(aligned_sqs["s1"]), bubble_endpoint))
    sorted_breaks = sorted(breaks, key=lambda x: x[0])
    breakings_in_aligned = groupby(sorted_breaks, key=lambda x: x[0])

    # pick the more abundant path as a whole
    more_abundant = get_more_abundant(G, path1, path2)

    last_vertex, last_index = None, None
    for index, group in breakings_in_aligned:
        items = list(group)
        # these vertices end on the same point in the alignment AND after them a vertex must end
        # for each break find all vertices that end there (that do not necessarily require break)
        ending_vertices = reversed_endings[index]
        new_vertex = min([v for v in ending_vertices])

        new_sq, new_count = get_new_node_info(aligned_sqs, index, last_index,
                                              more_abundant, counts1, counts2, k, contracted_info)
        newG.add_node(new_vertex, sq=new_sq, counts=new_count, length=len(new_sq))
        # link with the node added in previous step
        if last_vertex is not None:
            newG.add_edge(last_vertex, new_vertex)
        # get other edges
        for _, v in items:
            for in_edge_in, _ in G.in_edges(v):
                if in_edge_in not in changed:
                    newG.add_edge(in_edge_in, new_vertex)
            for _, out_edge_out in G.out_edges(v):
                if out_edge_out not in changed:
                    newG.add_edge(new_vertex, out_edge_out)
        last_vertex = new_vertex
        last_index = index

    # create the path and add it to the new graph
    return newG


class dijkstraer:
    def __init__(self, G, max_bubble):
        self.graph = G
        self.max_bubble = max_bubble + 1

    def calc_distances(self, s):
        dists = nx.single_source_dijkstra_path_length(self.graph,
                                                      s,
                                                      weight="extension",
                                                      cutoff=self.max_bubble)
        # dists = sorted(dists.items(), key=lambda x: x[1])
        return dists


def identify_contract_bubble(G, k, max_bubble_length, contracted_info, all_dists):
    # returns a new graph with contracted
    discovered = {}  # node : set
    closed = {}  # node : set
    closed_edges = {}  # node : list of edges

    for from_, to_ in G.edges:
        if G.nodes[to_]["length"] >= k:
            G.edges[(from_, to_)]["extension"] = G.nodes[to_]["length"] - k + 1
        else:
            raise Exception(f"Insertion of node with less than k={k} bases.")
            # print()
            # G.edges[(from_, to_)]["extension"] = None

    reasonable = [n for n, deg in G.out_degree() if deg > 1]
    # reasonable = [n for n in reasonable if G.in_degree(n) > 0]  # sufficient starting points

    if len(reasonable) == 0:
        return None, None

    # start = time.time()
    to_dijkstra = []
    for source in reasonable:
        discovered[source] = set()
        closed[source] = set()
        closed_edges[source] = []
        if source not in all_dists:
            to_dijkstra.append(source)

    d = dijkstraer(G, max_bubble_length)

    with Pool(8) as pool:
        distances = pool.map(d.calc_distances, to_dijkstra)
    for source, dists in zip(to_dijkstra, distances):
        all_dists[source] = dists

    reasonable = set(reasonable)

    to_process = []
    for item in all_dists:
        to_process.extend(
            [(item, to, dst) for to, dst in all_dists[item].items()]
        )
    to_process = sorted(to_process, key=lambda x: x[2])

    for node, current, distance in to_process:
        closed[node].add(current)
        # print(f"DST:{distance}, {node}")
        for _, next_ in G.edges(current):
            # print(f"\t {next_}")
            closed_edges[node].append((current, next_))
            if next_ in discovered[node]:
                contracted_graph = contract_bubble(G, k, next_, closed[node], max_bubble_length, contracted_info, closed_edges[node])
                if contracted_graph is None:
                    return contracted_graph, None
                if contracted_graph:
                    # remove affected vertices
                    all_dists.pop(next_, None)
                    for affected in closed[node]:  # these are definitely affected
                        all_dists.pop(affected, None)
                    to_remove = set()
                    for n in all_dists:
                        for m, dst in all_dists[n].items():
                            if m in closed[node]:
                                to_remove.add(n)
                                break
                    for n in to_remove:
                        all_dists.pop(n, None)

                    # print(f"identification done: {time.time() - start}")
                    return contracted_graph, all_dists
            else:
                discovered[node].add(next_)

    return None, None


def signum(no):
    if no < 0:
        return "-"
    elif no > 0:
        return "+"
    else:
        return ""


emptyline = "".join([" " for i in range(200)])


def process_component(filename, contracted_info, output_file, args, bub_done, comp_done):
    # print(filename)
    try:
        k = args.k
        if args.max_bubble_sequence is None:
            max_bubble_length = 2 * k
        else:
            max_bubble_length = args.max_bubble_sequence
            # max length of sequence represented in a bubble

        orig_graph = prep.debruijn_fasta_to_digraph(filename, extended_data=True)
        # stuff is being added to this with contractions of paths -- info on errors and variability

        graph = orig_graph
        graph.remove_edges_from(nx.selfloop_edges(graph))
        output_graph = graph
        all_dists = {}

        bubbles = 0
        while graph is not None:
            new_graph, all_dists = identify_contract_bubble(graph, args.k, max_bubble_length,
                                                            contracted_info, all_dists)
            print(f"now processing {filename}; already processed components: {comp_done}, corrected bubbles {bub_done + bubbles} ({bubbles} in current)", end='\r')
            # if graph is none, keep the old graph as the output one
            if new_graph is not None:
                output_graph = new_graph
                bubbles += 1
            graph = new_graph
    except Exception as e:
        print(f"Error occurred at {filename}. Exception is being raised.")
        raise e

    # write output_graph to file
    with open(output_file, mode='w') as writer:
        for node in output_graph.nodes:
            sq = output_graph.nodes[node]["sq"]
            counts = output_graph.nodes[node]["counts"]
            total_count = sum(counts)
            str_counts = ",".join([str(x) for x in counts])
            edges = []
            for _, to_ in output_graph.out_edges(node):
                sign_from = signum(node)
                sign_to = signum(to_)
                abs_to = abs(to_)
                edges.append(
                    f"L:{sign_from}:{abs_to}:{sign_to}"
                )

            edges = " ".join(edges)

            id = f"{node} LN:i:{len(sq)} KC:i:{total_count} km:f:{total_count / len(counts)} {str_counts} {edges}"

            print(f">{id}", file=writer)
            print(f"{sq}", file=writer)

    return bubbles


def main(args):
    input_directory = args.raw_components
    output_directory = args.corrected_components
    os.makedirs(output_directory)
    contracted_info_filename = args.replaced
    contracted_info = open(contracted_info_filename, mode="w")
    print(f"about\tarchived_sq\tmapped_to_sq\tcounts_from_total\tedit distance", file=contracted_info)

    files = os.listdir(input_directory)

    repaired_bubble = 0
    repaired_comps = 0

    # with bar("correcting errors...", max=len(files)) as progressbar:
    for comp in files:
        filename = f"{input_directory}/{comp}"
        output_file = f"{output_directory}/corr_{comp}"
        repaired_bubble += process_component(filename, contracted_info, output_file, args, repaired_bubble, repaired_comps)
        repaired_comps += 1

    print("Error correction finished.")


if __name__ == '__main__':
    args = parser.parse_args([] if "__file__" not in globals() else None)
    main(args)
