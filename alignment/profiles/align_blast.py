from Bio.Blast.Applications import NcbiblastxCommandline
import os
import pandas as pd

working_dir = "../.."
profiles_file = "aligned"
low_abund_file = "low_abund_kmers.fasta"

# low_abund to fasta --
#df_low = pd.read_csv(os.path.join(working_dir, low_abund_file), header=None, sep=';')
#df_low.columns = ["kmer", "counts"]
#df_low["counts"] = df_low["counts"].str.split(",").str[0].astype(int)

# with open("low_abund_kmers.fasta", mode="w") as writer:
#    for i,row in df_low.iterrows():
#        print(f">{row.counts}", file=writer)
#        print(row.kmer, file=writer)


# Define paths to input and output directories -- the same: the working directory

# Define paths to input and output files
#query = os.path.join(working_dir, 'k_sp_CB01950_penicillin.fasta')           # query sequence(s)
#db = os.path.join(working_dir, 'kitasatospora_proteins.faa')                 # BLAST database
#blastout = os.path.join(working_dir, 'AMK19_00175_blastx_kitasatospora.tab')  # BLAST output

#cmd_blastx = NcbiblastxCommandline(query=query, out=blastout, outfmt=6, db=db)
#results = pd.read_csv(blastout, sep="\t", header=None)