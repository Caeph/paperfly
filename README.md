# Paperfly: a software for analysis of ChIP-seq (or similar) experiments without the reference genome

The PseudoAssembly based Peak Finder (PAPerFly) assembles the sequencing reads seen during a ChIP-seq or similar experiment. 
This repository is still under construction and some additions will be made to the code.

# Requirements
- python 3.8 or newer, python3.8-venv to create a virtual environment
- jellyfish -- can be downloaded from Linux repositories or compiled from https://github.com/gmarcais/Jellyfish
- bcalm -- can be downloaded from Linux repositories or compiled from https://github.com/GATB/bcalm
- mono, msbuild and nuget for building and running the c# applications (msbuild can usually be downloaded as a part mono-complete packages)
- graphviz for drawing (only sfdp is used)

All of these are available from the Linux repositories. All of these must be in your path for the program to function properly.

# Instalation
Clone the repository. Use the included makefile to 
- create a python virtual environment and install the required packages via the venv pip
- build the C# binaries

Just run 
```
make
```
and then update your path. 
```
export PATH=$PATH:<paperfly directory>
```

The core of paperfly is a bash script.

# Running
To run the program when it is in your path, use
```
paperfly --input_fast<q|a> <input file> --k <k>
```
The input file can be in FASTA, FASTQ or gzipped FASTQ format (then use the fastq option). The value of K corresponds to the lenght of k-mers in the de Bruijn graph of the reads.

Other parameters are:
- ```--control_file <filename>``` 
- ```--minimal_abundance <value>```
- ```--minimal_abundance_percentile <value>```
- ```--working_directory <path>```
- ```--draw```
- ```--allowed_misses <value>```
- ```--distance_metric <hamming|levenstein>```
- ```--window <value>```
- ```--prominence <value>```
- ```--peak_format <consensus|sq_count>```
- ```--peak_min_width <value>```
- ```--peak_max_width <value>```

# Reference
The article describing the method was not published yet, but we are working on that. Both BCALM and Jellyfish will be properly cited in the publication.
