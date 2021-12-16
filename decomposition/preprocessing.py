import networkx as nx
import subprocess
import pandas as pd
import pyfastaq.sequences as pyfs
import shutil
import os
import numpy as np

bcalm_path = "bcalm"
jellyfish_path = "jellyfish"
sfdp_path = "sfdp"
drawing_timeout = 20
jellyfish_settings = "-s 1G -t 10 -C"


def store_low_abundace_kmers(counts_csv, min_abundance, k, directory=None):
    # identify kmers present in data
    # jellyfish count -m $k -s 100M -t 10 -o mers.unitigs.$k.jf $unitigs
    if (directory is not None) and (len(directory) > 0):
        dir = directory + '/'
    else:
        dir = ""

    counts_df = pd.read_csv(counts_csv, sep=' ', header=None)
    counts_df.columns = ['kmer', 'count']
    counts_df = counts_df[(counts_df['count'] < min_abundance) & (counts_df['count'] > 0)]
    counts_df['count'] = counts_df['count'].map(lambda x : ','.join([str(1) for x in range(k)]))

    low_abund_name = f"{dir}low_abund_kmers.csv"
    counts_df.to_csv(low_abund_name, header=False, index=False, sep=';')
    return low_abund_name


def get_control_counts(control, k, flag="unitigs", directory=None):
    print(f"running jellyfish on {control}")
    if (directory is not None) and (len(directory) > 0):
        dir = directory + '/'
    else:
        dir = ""

    control_mers = f"{dir}control.{k}.fa"

    # jellyfish count -m $k -s 1G -t 10 -o control_counts.tmp --if mers.unitigs.$k.fa $control
    process = subprocess.run(
        [
            f"{jellyfish_path} count -m {k} -s 1G -t 10 -o {control_mers}.tmp {control}"],
        shell=True)
    evaluate_process(process)

    # jellyfish dump -c counts.tmp > $outfile
    process = subprocess.run(
        [f"{jellyfish_path} dump -c {control_mers}.tmp > {control_mers}"],
        shell=True)
    evaluate_process(process)

    return control_mers


def evaluate_process(completed, continue_if_bad=False):
    if not continue_if_bad:
        completed.check_returncode()
        print("process finished succesfully.")
    else:
        if completed.returncode == 0:
            print("process finished succesfully.")
        else:
            print("error in process, skipping and continuing")


def run_sfdp(dot, dot_pof, timeout=drawing_timeout):
    try:
        # process = subprocess.run([f"sfdp_path -x -Goverlap=scale -Tpdf {dot} >{dot_pof}"],
        #                         shell=False, timeout=timeout)
        process = subprocess.run([sfdp_path, "-x", "-Goverlap=scale", "-Tpdf", f"{dot}", "-o", f"{dot_pof}"],
                                 shell=False, timeout=timeout)
        evaluate_process(process, continue_if_bad=True)
    except subprocess.TimeoutExpired:
        print("Timeout expired, skipping...")


def run_bcalm(input, k, abu, stub=False):
    print(f"running BCALM on {input}")
    outputname = input.split('/')[-1].split('.')[0] + ".unitigs.fa"
    if stub:
        return outputname
    process = subprocess.run([f"{bcalm_path} -in {input} -kmer-size {k} -abundance-min {abu}"],
                             shell=True)
    evaluate_process(process)
    return outputname


def get_kmers_for_initials(mersfile):
    df = pd.DataFrame(
        [entry.seq for entry in pyfs.file_reader(mersfile)]
    )
    df.columns = ["present"]

    for base in ['A', 'C', 'G', "T"]:
        df["next_" + base] = df["present"].str[1:] + base
        df["prev_" + base] = df["present"].str[:-1] + base

    orig = df["present"]
    inits = df.drop(columns=["present"]).melt()["value"].drop_duplicates()

    df = pd.merge(orig, inits, left_on="present", right_on="value", how="right")
    return df[df["present"].isna()]["value"].values


def run_jellyfish_bcalm_kmers_canonical(fasta, kmers_file, k, flag="unitigs_canonical", directory=None, stub=False,
                                        cleanup=True):
    print(f"running jellyfish on {fasta}")
    outputname = kmers_file + ".counts.csv"

    # jellyfish count -m $k -s 100M -t 10 -o mers.unitigs.$k.jf $unitigs
    if (directory is not None) and (len(directory) > 0):
        dir = directory + '/'
    else:
        dir = ""

    process = subprocess.run([f"{jellyfish_path} count -m {k} -s 100M -t 10 -o {dir}mers.{flag}.{k}.jf {kmers_file}"],
                             shell=True)
    evaluate_process(process)

    # jellyfish dump mers.unitigs.$k.jf > mers.unitigs.$k.fa
    process = subprocess.run([f"{jellyfish_path} dump {dir}mers.{flag}.{k}.jf > {dir}mers.{flag}.{k}.fa"],
                             shell=True)
    evaluate_process(process)

    # jellyfish count -m $k -s 1G -t 10 -o counts.tmp --if mers.unitigs.$k.fa $infile
    process = subprocess.run(
        [
            f"{jellyfish_path} count -m {k} {jellyfish_settings} -o {dir}counts.{flag}.tmp --if {dir}mers.{flag}.{k}.fa {fasta}"],
        shell=True)
    evaluate_process(process)

    # jellyfish dump -c counts.tmp > $outfile
    process = subprocess.run(
        [f"{jellyfish_path} dump -c {dir}counts.{flag}.tmp > {outputname}"],
        shell=True)
    evaluate_process(process)
    # cleanup

    if cleanup:
        for fl in [f"{dir}counts.{flag}.tmp", f"{dir}mers.{flag}.{k}.fa", f"{dir}mers.{flag}.{k}.jf"]:
            os.remove(fl)

    return outputname


def run_jellyfish_bcalm_kmers(fasta, unitigs, k, flag="unitigs", directory=None, stub=False, cleanup=True, control=None, control_flag="control"):
    print(f"running jellyfish on {fasta}")
    outputname = unitigs + ".counts.csv"
    if stub:
        return outputname

    # get k-mers that are found in the unitigs file
    # jellyfish count -m $k -s 100M -t 10 -o mers.unitigs.$k.jf $unitigs
    if (directory is not None) and (len(directory) > 0):
        dir = directory + '/'
    else:
        dir = ""
    process = subprocess.run([f"{jellyfish_path} count -m {k} -s 100M -t 10 -o {dir}mers.{flag}.{k}.jf {unitigs}"],
                             shell=True)
    evaluate_process(process)

    # jellyfish dump mers.unitigs.$k.jf > mers.unitigs.$k.fa
    process = subprocess.run([f"{jellyfish_path} dump {dir}mers.{flag}.{k}.jf > {dir}mers.{flag}.{k}.fa"],
                             shell=True)
    evaluate_process(process)

    # count the number of k-mers of interest
    # jellyfish count -m $k -s 1G -t 10 -o counts.tmp --if mers.unitigs.$k.fa $infile
    process = subprocess.run(
        [
            f"{jellyfish_path} count -m {k} -s 1G -t 10 -o {dir}counts.{flag}.tmp --if {dir}mers.{flag}.{k}.fa {fasta}"],
        shell=True)
    evaluate_process(process)

    # jellyfish dump -c counts.tmp > $outfile
    process = subprocess.run(
        [f"{jellyfish_path} dump -c {dir}counts.{flag}.tmp > {outputname}"],
        shell=True)
    evaluate_process(process)

    if cleanup:
        for fl in [f"{dir}counts.{flag}.tmp", f"{dir}mers.{flag}.{k}.fa", f"{dir}mers.{flag}.{k}.jf"]:
            os.remove(fl)

    return outputname


def divide_to_components(filtered_unitigs, comp_dir="components", canonical=False):
    print(f'dividing {filtered_unitigs} to components, dir {comp_dir}')
    remove_dir_if_exists(comp_dir)
    G = debruijn_fasta_to_digraph(filtered_unitigs, canonical=canonical)

    components = nx.algorithms.components.weakly_connected_components(G)

    dict_ = {}
    i = 0
    singles = 0
    for comp in components:
        if len(comp) <= 1:
            singles += 1
            continue

        for node in comp:
            if node in dict_:
                print("DUPLICATION!", node)
                raise Exception("Node duplication.")
            # should not happen, weakly connected components should be disjuct
            dict_[node] = i
        i += 1

    print(f"found {singles} isolated nodes and ~{i} weakly connected components")  # the i is not necessarily precise
    del G  # memory efficiency
    del components

    ##########################################################################################
    # write components into respective files

    os.makedirs(comp_dir, exist_ok=False)
    print("Directory now exists.")

    unitigs = pyfs.file_reader(filtered_unitigs)
    for entry in unitigs:
        items = [x for x in entry.id.split(' ') if len(x) > 0]
        id_ = int(items[0])

        if id_ in dict_:
            flname = "%s/%i.comp.fasta" % (comp_dir, dict_[id_])
            file = open(flname, mode='a')
        else:
            flname = "%s/singles.comp.fasta" % (comp_dir)
            file = open(flname, mode='a')
        file.write(">%s\n" % entry.id)
        file.write("%s\n" % entry.seq)
        file.close()
    return comp_dir


def reverse(seq):
    return seq.replace('A', 'x'
                       ).replace('T', 'A'
                                 ).replace('x', 'T'
                                           ).replace('C', 'x'
                                                     ).replace('G', 'C'
                                                               ).replace('x', 'G')[::-1]


def replace_zero_index(unitigs):
    print(f"replacing zero index on {unitigs}")
    outfilename = unitigs + ".moved.fa"
    fasta = pyfs.file_reader(unitigs)
    outfile = open(outfilename, mode='w')

    for entry in fasta:
        items = entry.id.split(' ')
        id_ = int(items[0])
        new_id = id_ + 1  # zero cannot be an id
        edges = [x for x in items[6:] if len(x) > 0]
        # print('edges', edges)

        new_edges = ["L:%s:%i:%s" % (e.split(':')[1], int(e.split(':')[2]) + 1, e.split(':')[3]) for e in edges]
        outfile.write(">%i %s %s\n" % (new_id, " ".join(items[1:5]), " ".join(new_edges)))
        outfile.write("%s\n" % entry.seq)
    outfile.close()
    return outfilename


def expand(input):
    print(f"expanding {input}")
    outfilename = input + ".expanded.fa"
    fasta = pyfs.file_reader(input)

    # edges_buffer = open("EDGES.tmp", mode='w')
    outfile = open(outfilename, mode='w')

    for entry in fasta:
        items = entry.id.split(' ')
        id_ = int(items[0])
        edges = [x for x in items[5:] if len(x) > 0]

        for sign in ['+', '-']:
            if sign == '-':
                new_id = -id_
                new_sq = reverse(entry.seq)
            else:
                new_id = id_
                new_sq = entry.seq
            edges_subset = [e for e in edges if e[2] == sign]
            outfile.write(">%i %s %s\n" % (new_id, " ".join(items[1:5]), " ".join(edges_subset)))
            outfile.write(new_sq + "\n")
    return outfilename


def link_jellyfish_counts(input, jf_input, k):
    print(f"linking {jf_input} to {input}")
    outputname = input + ".expanded_counts.fa"

    # load description
    desc = pd.read_csv(jf_input, sep=' ', header=None)
    desc.columns = ['kmer', 'count']

    fasta = pyfs.file_reader(input)
    seqs = []

    for entry in fasta:
        items = entry.id.split(' ')
        ln = int(items[1].split(':')[-1])
        id_ = items[0]
        edges = " ".join(items[4:])
        del items

        if ln == k:
            seqs.append([id_, ln, edges, entry.seq, entry.seq])
        else:
            for i in range(ln - k + 1):
                seqs.append([id_, ln, edges, entry.seq[i:i + k], entry.seq])

    df = pd.DataFrame(seqs, columns=['unitig_id', 'len', 'edges', 'kmer', 'seq'])

    # find and calculate noncanonical counts
    m = pd.merge(
        desc, df, how='right', on='kmer'
    ).fillna(0)
    m = m.groupby(
        by=['unitig_id', 'len', 'edges', 'seq']).agg(
        {'count': ['sum', 'mean', lambda x: ','.join(["%.0f" % i for i in x])]})
    # del df
    del desc
    m.columns = ['_'.join(col).strip().replace('<lambda_0>', 'all') for col in m.columns.values]
    m.reset_index(inplace=True)

    # write output
    outfl = open(outputname, mode='w')
    for x in m.values:
        outfl.write('>%s LN:i:%i KC:i:%i km:f:%.1f %s %s\n' % (x[0], x[1], x[4], x[5], x[6], x[2]))
        outfl.write(x[3] + '\n')
    outfl.close()
    return outputname


def filter_abundance_simple(input, tmp_outfile, req_abu):
    removing = set()

    fasta = pyfs.file_reader(input)
    for entry in fasta:
        items = entry.id.split(' ')
        unitig_abu = int(items[2].split(':')[2])
        if unitig_abu < req_abu:
            removing.add(int(items[0]))
            continue
        counts = [int(x) for x in items[4].split(',')]
        if len(counts) == len([1 for x in counts if x < req_abu]):  # every unitig is underabundant
            removing.add(int(items[0]))
    # print(removing)

    outfile = open(tmp_outfile, mode='w')
    fasta = pyfs.file_reader(input)
    for entry in fasta:
        id_ = int(entry.id.split(' ')[0])
        if id_ in removing:
            continue
        # print(entry.id)

        edges = [x for x in entry.id.split(' ')[6:] if len(x) > 0]
        edges = [x for x in edges if int("%s%s" % (x.split(':')[3], x.split(':')[2])) not in removing]
        outfile.write(">%s\n" % (" ".join([*entry.id.split(' ')[:6], *edges])))
        outfile.write("%s\n" % entry.seq)

    outfile.close()
    ### \first walkthrough


def filter_abundance(input, req_abu, k, tmp_outfile='filtered_unitigs.tmp',
                     buffer_name="filtering_buffer.tmp"):
    print(f"filtering abundance on {input}")
    outputname = input + ".filtered.fa"

    filter_abundance_simple(input, tmp_outfile, req_abu)

    ### unitig filtering
    fasta = pyfs.file_reader(tmp_outfile)
    abs_max_id = -1

    i = 0

    remove_before = set()  # original id, TO is being deleted
    split_unitigs = set()  # original id, FROM is being deleted

    buffer = open(buffer_name, mode='w')

    for entry in fasta:
        items = entry.id.split(' ')
        counts = [int(x) for x in items[4].split(',')]
        id_ = int(items[0])

        if abs(id_) > abs_max_id:
            abs_max_id = abs(id_)

        if len([x for x in counts if x < req_abu]) <= 0:
            continue

        # budeme splitovat
        split_unitigs.add(id_)
        edges = [x for x in items[4:] if x[0:2] == "L:"]

        prev = []
        new_tigs = []

        for c, j in zip(counts, range(len(counts))):
            # print(counts[j])
            if c < req_abu:
                if len(prev) > 0:
                    new_tigs.append(prev)
                prev = []
            else:

                prev.append(j)
            # print(prev)
        if len(prev) > 0:
            new_tigs.append(prev)
        # print(new_tigs)

        remove_after = False
        if 0 not in new_tigs[0]:
            remove_before.add(id_)
        if (len(counts) - 1) not in new_tigs[-1]:
            remove_after = True
        ##
        for tig, j in zip(new_tigs, range(len(new_tigs))):
            if (j == 0) and not remove_after:
                edges = " ".join(edges)
            else:
                edges = ''

            counts_tig = counts[tig[0]:tig[-1] + 1]
            tig_sq = entry.seq[tig[0]:(k + tig[-1])]

            buffer.write(">%i_%i LN:i:%i KC:i:%i km:f:%.1f %s %s\n" % (id_,
                                                                       j,
                                                                       len(tig_sq),
                                                                       sum(counts_tig),
                                                                       sum(counts_tig) / len(tig_sq),
                                                                       ",".join([str(x) for x in counts_tig]),
                                                                       edges
                                                                       ))
            buffer.write("%s\n" % tig_sq)
        i += 1
    print("no of unitigs split: %i" % i)
    buffer.close()

    output = open(outputname, mode='w')
    fasta = pyfs.file_reader(tmp_outfile)
    for entry in fasta:
        items = entry.id.split(' ')
        id_ = int(items[0])
        if id_ in split_unitigs:
            continue

        edges = [x for x in items[4:] if x[0:2] == "L:"]
        new_edges = [e for e in edges if int('%s%s' % (e.split(':')[3], e.split(':')[2])) not in remove_before]
        output.write(">%i %s %s\n" % (id_,
                                      ' '.join(items[1:5]),
                                      ' '.join(new_edges)
                                      ))
        output.write("%s\n" % entry.seq)

    fasta = pyfs.file_reader(buffer_name)
    from collections import Counter
    c = Counter()
    for entry in fasta:
        items = entry.id.split(' ')
        orig_id = int(items[0].split('_')[0])
        unitig_num = int(items[0].split('_')[1])

        if unitig_num == 0:
            new_id = orig_id
        else:
            new_id = abs_max_id + 1
            abs_max_id = new_id
            c[new_id] += 1

        edges = [x for x in items[4:] if x[0:2] == "L:"]
        new_edges = [e for e in edges if int('%s%s' % (e.split(':')[3], e.split(':')[2])) not in remove_before]

        output.write(">%i %s %s\n" % (new_id,
                                      ' '.join(items[1:5]),
                                      ' '.join(new_edges)
                                      ))
        output.write("%s\n" % entry.seq)
    for fl in [tmp_outfile, buffer_name]:
        os.remove(fl)
    return outputname

def draw_components(compdir, req_abu, picdir=None, canonical=False, skip_linear=False):
    if picdir is None:
        picdir = compdir
        if picdir[-1] == '/':
            picdir = picdir[:-1]
        picdir = picdir + "_pictures"
    remove_dir_if_exists(picdir)
    os.makedirs(picdir, exist_ok=False)
    for file in os.listdir(compdir):
        print(f"drawing {file}")
        if file == "singles.comp.fasta":
            print(file)
        debruijn_fasta_to_dot(f"{compdir}/{file}", f"{picdir}/{file}.dot", f"{picdir}/{file}.pdf", req_abu,
                              emptylabel=True, canonical=canonical, skip_linear=skip_linear)


def remove_dir_if_exists(dirname):
    if os.path.isdir(dirname):
        print(f'removing existing directory {dirname}')
        shutil.rmtree(dirname)


def debruijn_fasta_to_dot(debruijn_fl, dotname, pdfname, req_abu, draw=True, remove_dot=True, write_node_labels=True,
                          ranksep=6, nodesep=1.2, splines="false",
                          emptylabel=False, shape='point', nodecolor='blue',
                          arrowsize=0.3, manipulate_edgewidth=True, fontsize=12, canonical=False, skip_linear=False,
                          ):
    dot = open(dotname, mode='w')
    dot.write("digraph {\n")
    dot.write(f"graph[ranksep={ranksep}, nodesep={nodesep}, splines={splines}]\n")

    if emptylabel:
        labelattr = "label=\"\","
    else:
        labelattr = ''

    dot.write(f"node[{labelattr} shape={shape}, color={nodecolor}, fontsize={fontsize}]\n")
    dot.write(f"edge[arrowsize={arrowsize}]\n")
    written = infer_graph_layout(debruijn_fl, dot, write_nodes=write_node_labels, manipulate_edgewidth=manipulate_edgewidth,
                       req_abu=req_abu, canonical=canonical, skip_linear=skip_linear)
    dot.write('}\n')
    dot.close()

    if (draw & written):
        run_sfdp(dotname, pdfname)
    if remove_dot:
        os.remove(dotname)


def split_unitigs(unitigs_counts, split_unitigs_filename, k):
    max_ident = - np.inf

    for entry in pyfs.file_reader(unitigs_counts):
        ident = int(entry.id.split(' ')[0])
        if ident > max_ident:
            max_ident = ident

    id = max_ident + 1

    with open(split_unitigs_filename, mode='w') as writer:
        for entry in pyfs.file_reader(unitigs_counts):
            desc = entry.id.split(' ')
            orig_id = desc[0]
            counts = [int(x) for x in desc[4].split(',')]

            edges = " ".join(desc[5:])
            if len(counts) == 1:
                count = counts[0]
                kmer = entry.seq
                writer.write(f">{orig_id} LN:i:{k} KC:i:{count} km:f:{float(count)} {count} {edges}\n")
                writer.write(f"{kmer}\n")
                continue

            for i, count in enumerate(counts):
                kmer = entry.seq[i:i + k]
                if i == 0:
                    if orig_id[0] == '-':
                        sign = '-'
                    else:
                        sign = '+'
                    to_next = f"L:{sign}:{id}:+"  # edge to next
                    writer.write(f">{orig_id} LN:i:{k} KC:i:{count} km:f:{float(count)} {count} {to_next}\n")
                    writer.write(f"{kmer}\n")
                    continue
                if i == (len(counts) - 1):
                    writer.write(f">{id} LN:i:{k} KC:i:{count} km:f:{float(count)} {count} {edges}\n")
                    writer.write(f"{kmer}\n")
                    id += 1
                    continue
                # middles
                to_next = f"L:+:{id + 1}:+"
                writer.write(f">{id} LN:i:{k} KC:i:{count} km:f:{float(count)} {count} {to_next}\n")
                writer.write(f"{kmer}\n")
                id += 1


def infer_graph_layout(debruijn_fl, dothandle, write_nodes=False, manipulate_edgewidth=False, req_abu=100,
                       max_edge_width=15, canonical=False, skip_linear=False):
    fasta = pyfs.file_reader(debruijn_fl)

    if skip_linear:
        graph = debruijn_fasta_to_digraph(debruijn_fl)
        a = max(graph.in_degree(graph.nodes), key=lambda x: x[1])[1]
        b = max(graph.out_degree(graph.nodes), key=lambda x: x[1])[1]

        if ((a <= 1) & (b <= 1)):
            return False

    for entry in fasta:
        items = [x for x in entry.id.split(' ') if len(x) > 0]
        id_ = items[0].replace('-', 'rev_')  # if canonical, nothing is replaced
        edges = [x for x in items if x[0:2] == 'L:']
        if write_nodes:
            if len(items) > 4:
                dothandle.write(f"{id_} [label=\"mincount={items[4]},id={items[0]}\\n{entry.seq}\"]\n")
            else:
                dothandle.write(f"{id_} [label=\"mincount={items[3].split(':')[-1]},id={items[0]}\\n{entry.seq}\"]\n")

        for e in edges:
            if manipulate_edgewidth:
                edgewidth = "[penwidth=%.2f]" % min(float(items[3][5:]) / req_abu, max_edge_width)  # mean count
            else:
                edgewidth = ''
            spl = e.split(':')
            if canonical:
                to = spl[2]
            else:
                to = ("%s%s" % (spl[3], spl[2])).replace('-', 'rev_').replace('+', '')
            dothandle.write("%s -> %s %s;\n" % (id_, to, edgewidth))

    return True


def debruijn_fasta_to_digraph(fastafile, extended_data=False, min_count=False, canonical=False):
    unitigs = pyfs.file_reader(fastafile)
    G = nx.DiGraph()

    for entry in unitigs:
        items = [x for x in entry.id.split(' ') if len(x) > 0]
        id_ = int(items[0])
        # print(items)

        if extended_data:
            G.add_node(id_,
                       sq=entry.seq,
                       length=int(items[1][5:]),
                       counts=[int(x) for x in items[4].split(',')]
                       )
        else:
            G.add_node(id_,
                       sq=entry.seq,
                       )
        if (not extended_data) and min_count:
            G.nodes[id_]['min_count'] = min([int(x) for x in items[4].split(',')])

        edges = [x for x in items if x[0:2] == 'L:']

        if not canonical:
            for edge in edges:
                # print(edge)
                spl = edge.split(':')
                to = int("%s%s" % (spl[3], spl[2]))
                G.add_edge(id_, to)
        else:
            for edge in edges:
                # print(edge)
                spl = edge.split(':')
                to = int(spl[2])
                G.add_edge(id_, to)

    return G
