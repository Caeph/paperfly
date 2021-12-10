import numpy as np
import os
from pyfastaq import sequences as pyfs
from pyfastaq import tasks as pyftasks
import pandas as pd
from collections import Counter
import igraph


def debruijn_fasta_to_igraph(fastafile):
    unitigs = pyfs.file_reader(fastafile)
    G = igraph.Graph().as_directed()

    ids = []
    min_counts = []
    all_edges = []

    for entry in unitigs:
        items = [x for x in entry.id.split(' ') if len(x) > 0]
        id_ = items[0]
        min_counts.append(min([int(x) for x in items[4].split(',')]))
        edges = [x for x in items if x[0:2] == 'L:']

        G.add_vertex(name=id_,
                     )

        edges_to_add = [(id_, "%s%s" % (e.split(':')[3].replace('+', ''), e.split(':')[2])) for e in edges]
        all_edges.extend(edges_to_add)

    G.add_edges(all_edges)
    G.vs["min_count"] = min_counts
    return G


def strongly_connected_components_description(dir):
    connected = 0
    single = 0

    for file in os.listdir(dir):
        comp_graph = os.path.join(dir, file)
        G = debruijn_fasta_to_igraph(comp_graph)
        scc = G.clusters(mode="strong")
        sizes = np.array(scc.sizes())
        big_scc = set(np.where(sizes > 1)[0])
        if len(big_scc) > 0:
            #print(comp_graph, [(x, sizes[x]) for x in big_scc])
            big_sizes = [str(sizes[x]) for x in big_scc]
            print(f"{comp_graph} has {len(big_sizes)} cycle(s) of sizes", ", ".join(big_sizes))
            connected += 1
        else:
            single += 1
    print(f"{connected} weakly connected components in {dir} contain cycle.")


def covered_base_number(assembly_dir,
                        orig_data="../../dataset_foxo_19_11/Selex2_S1_L001_R1_001.fasta"):
    fasta = pyfs.file_reader(orig_data)
    base_number_orig = 0
    for entry in fasta:
        base_number_orig += len(entry.seq)
    print(f"Original file: {base_number_orig}")

    base_number_assem = 0
    for csv in os.listdir(assembly_dir):
        df = pd.read_csv(os.path.join(assembly_dir, csv), sep=';', header=None)
        df.columns = ['sq', 'count']
        lengths = df['sq'].apply(lambda row: len(row))
        counts = lengths * df["count"]
        base_number_assem += counts.sum()
    print(f"Assembled files: {base_number_assem}")
    print(f"Coverage of data: {base_number_assem / base_number_orig}%")


def component_size_histogram(dir):
    sq_nos = Counter()
    for file in os.listdir(dir):
        comp_file = os.path.join(dir, file)
        no = pyftasks.count_sequences(comp_file)
        sq_nos[no] += 1
    return sq_nos


def get_candidates_for_separate_learning(dir, min_node_count=100, min_forking_part=0.3, for_view=True):
    candidates = []
    for file in os.listdir(dir):
        comp_graph = os.path.join(dir, file)
        G = debruijn_fasta_to_igraph(comp_graph)

        if G.vcount() >= min_node_count:
            a = Counter(G.outdegree())
            forks = sum([a[i] for i in [2, 3, 4]])  # at least two sufficiently abundant variants follow

            # candidates for separate learning
            if forks / G.vcount() >= min_forking_part:
                if for_view:
                    candidates.append(os.path.join(dir.strip('/')+'_pictures', file+'.pdf'))
                else:
                    candidates.append(comp_graph)
    return candidates