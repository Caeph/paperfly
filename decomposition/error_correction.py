import networkx as nx
import preprocessing as prep
import edlib
from itertools import groupby
import argparse

# this is done as described in Velvet genome assembler -- Tour bus algorithm

parser = argparse.ArgumentParser()
parser.add_argument("--component",
                    default=None,
                    type=str,
                    help="Path to input: fasta of a weakly connected component."
                    )
parser.add_argument("--output",
                    default=None,
                    type=str,
                    help="Path to output: fasta of a corrected weakly connected component."
                    )
parser.add_argument("--k",
                    default=31,
                    type=int,
                    help="Size of the k-mer created by BCALM. Defaults at 31."
                    )
parser.add_argument("--replaced",
                    default=None,
                    type=str,
                    help="Path to a tsv with information on errore"
                    )


def main(args):
    filename = args.component
    k = args.k
    contracted_info_filename = args.replaced
    output_file = args.output

    # todo maybe variable threshold for seq similarity

    def archive(sq, mapped_to, counts):
        counts = ",".join([str(x) for x in counts])
        print(f"{sq}\t{mapped_to}\t{counts}", file=contracted_info)

    def sequence_from_path(G, path):
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

        # no gaps -- if best alignment contains gaps (I -- insertion, D -- deletion), dont contract
        if ("I" in align["cigar"]) or ("D" in align["cigar"]):
            print("gapped")
            print(edlib.getNiceAlignment(align, s1, s2))
            return 0, align

        similarity = 1 - (2 * align["editDistance"] / (len(s1) + len(s2)))  # frakce stejnych
        return similarity, align

    def get_endings_in_alignment(s1, s2, e1, e2, align):
        nice = edlib.getNiceAlignment(align, s1, s2)
        alignment = {
            "s1": nice["query_aligned"],
            "s2": nice["target_aligned"]
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

    def get_new_node_info(aligned, endindex, last_endindex, more_abundant, counts1, counts2):
        if last_endindex is None:
            init = 0
        else:
            init = last_endindex - k + 1

        represented1 = aligned["s1"][init:endindex]
        represented2 = aligned["s2"][init:endindex]

        kmers = endindex - init - k + 1
        a = len(represented1) - k + 1

        # no gaps expected
        c1, c2 = counts1[init:init + kmers], counts2[init:init + kmers]

        if represented1 == represented2:
            return represented1, c1

        new_counts = [x + y for x, y in zip(c1, c2)]
        new_sq = [represented1, represented2][more_abundant]

        if more_abundant == 0:
            archive(represented2, mapped_to=represented1, counts=c2)
        else:
            archive(represented1, mapped_to=represented2, counts=c1)

        return new_sq, new_counts

    def get_bubble_startpoint(subG, bubble_endpoint):
        stack = [bubble_endpoint]

        seen = set()
        # bubble_startpoint = None

        while len(stack) > 0:
            current = stack.pop(0)
            seen.add(current)
            for prev_, _ in subG.in_edges(current):
                if prev_ in seen:
                    bubble_startpoint = prev_
                    return bubble_startpoint
                stack.insert(0, prev_)
        return None

    def get_more_abundant(G, path1, path2):
        cs1, l1 = 0, 0
        for item in path1:
            cs1 += sum([x for x in G.nodes[item]["counts"]])
            l1 += len(G.nodes[item]["counts"])
        count_s1 = cs1 / l1
        cs2, l2 = 0, 0
        for item in path2:
            cs2 += sum([x for x in G.nodes[item]["counts"]])
            l2 += len(G.nodes[item]["counts"])
        count_s2 = cs2 / l2

        if count_s1 >= count_s2:
            return 0
        return 1

    def contract_bubble(G, bubble_endpoint, seen_vertices, thr=(1 - 1 / k)):
        seen_vertices.add(bubble_endpoint)
        subG = G.subgraph(seen_vertices)
        bubble_startpoint = get_bubble_startpoint(subG, bubble_endpoint)

        # print(bubble_startpoint, bubble_endpoint)
        paths_to_contract = list(nx.all_simple_paths(subG, bubble_startpoint, bubble_endpoint))
        # todo osetrit -- nemergnute bubbles delaji trochu bugr a mam tu vic cest
        if len(paths_to_contract) != 2:
            raise Exception("Wrong number of paths to contract.")

        # generate sequences, compare them
        path1, path2 = paths_to_contract[0], paths_to_contract[1]

        print("path1", path1)
        print("path2", path2)

        sq1, endings1 = sequence_from_path(G, path1)
        sq2, endings2 = sequence_from_path(G, path2)
        similarity, alignment = compare_sqs(sq1, sq2)
        if similarity < thr:
            return None

        print("contracted")
        print()

        aligned_sqs, endings_aligned = get_endings_in_alignment(sq1, sq2, endings1, endings2, alignment)
        # reversed endings
        reversed_endings = {}
        for vertex, index in endings_aligned.items():
            if index in reversed_endings:
                reversed_endings[index].append(vertex)
            else:
                reversed_endings[index] = [vertex]

        # generate new graph with contracted vertices
        changed = set(path1).union(set(path2))
        unchanged = [x for x in G.nodes if x not in changed]
        newG = G.subgraph(unchanged).copy()

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

        counts1, counts2 = [], []
        for v in path1:
            counts1.extend(G.nodes[v]["counts"])
        for v in path2:
            counts2.extend(G.nodes[v]["counts"])

        last_vertex, last_index = None, None
        for index, group in breakings_in_aligned:
            items = list(
                group)  # these vertices end on the same point in the alignment AND after them a vertex must end
            # for each break find all vertices that end there (that do not necessarily require break)
            ending_vertices = reversed_endings[index]

            new_vertex = min([v for v in ending_vertices])

            new_sq, new_count = get_new_node_info(aligned_sqs, index, last_index, more_abundant, counts1, counts2)
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

    def identify_contract_bubble(G):
        # returns a new graph with contracted
        def node_distance(from_, to_, edge_attr):
            # length of an extension in a path
            return G.nodes[to_]["length"] - k + 1

        all_dists = {}
        discovered = {}  # node : set
        closed = {}  # node : set
        for source in G.nodes:
            dists = nx.single_source_dijkstra_path_length(G, source, weight=node_distance)
            dists = sorted(dists.items(), key=lambda x: x[1])
            all_dists[source] = dists
            discovered[source] = set()
            closed[source] = set()

        for step in range(4 * len(G.nodes())):  # top estimate on no of steps -- from degree
            for node in G.nodes():
                if len(all_dists[node]) <= step:
                    continue
                # do step and check for bubble
                current, distance = all_dists[node][step]
                closed[node].add(current)

                for _, next_ in G.edges(current):
                    if (next_ in discovered[node]) and (next_ not in forbidden[node]):
                        contracted_graph = contract_bubble(G, next_, closed[node])
                        if contracted_graph is None:
                            forbidden[node].add(next_)
                        return contracted_graph, False
                    else:
                        discovered[node].add(next_)

        return None, True

    def signum(no):
        if no < 0:
            return "-"
        elif no > 0:
            return "+"
        else:
            return ""



    orig_graph = prep.debruijn_fasta_to_digraph(filename, extended_data=True)
    # stuff is being added to this with contractions of paths -- info on errors and variability
    contracted_info = open(contracted_info_filename, mode="w")
    print(f"archived_sq\tmapped_to_sq\tcounts_from_total", file=contracted_info)

    # get source
    graph = orig_graph

    # get source:
    # if component is a dag, get first
    output_graph = graph
    forbidden = {}  # storing endpoints of too dissimilar bubbles
    for item in graph.nodes:
        forbidden[item] = set()

    while graph is not None:
        new_graph, end = identify_contract_bubble(graph)
        if new_graph is not None:
            output_graph = new_graph
            graph = new_graph
        elif end:
            graph = new_graph
        # if graph is none and not end, keep the old graph

    # print(output_graph.edges)
    contracted_info.close()

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
                    f"L:{sign_from}:{to_}:{sign_to}"
                )

            id = f"{node} LN:i:{len(sq)} KC:i:{total_count} km:f:{total_count / len(counts)} {str_counts}"

            print(f">{id}", file=writer)
            print(f"{sq}", file=writer)


if __name__ == '__main__':
    args = parser.parse_args([] if "__file__" not in globals() else None)
    main(args)
