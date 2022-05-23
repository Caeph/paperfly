# Paperfly: a software for analysis of ChIP-seq (or similar) experiments without the reference genome

The PseudoAssembly based Peak Finder (PAPerFly) assembles the sequencing reads seen during a ChIP-seq or similar experiment. 
This repository is still under construction and some additions may be made to the code. The first release is an alpha version, a new release will be made when the work is published.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6379332.svg)](https://doi.org/10.5281/zenodo.6379332)


# Requirements
- python 3.8 or newer, python3.8-venv to create a virtual environment
- jellyfish -- can be downloaded from Linux repositories or compiled from https://github.com/gmarcais/Jellyfish
- bcalm -- can be downloaded from Linux repositories or compiled from https://github.com/GATB/bcalm
- mono, msbuild and nuget for building and running the C# applications (msbuild can usually be downloaded as a part mono-complete packages)
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

A Dockerfile is also available. If you are using it, please note that the recommended RAM is at least 8 GB. 
I also recommend you check it out if you are using Ubuntu 20.04, as you could run into problems when installing mono and msbuild.

# Running
To run the program when it is in your path, use
```
paperfly --input_fast<q|a> <input filename> --k <k> --control_file <control filename>
```
The input and the control file can be in FASTA, FASTQ or gzipped FASTQ format (then use the fastq option). The value of K corresponds to the lenght of k-mers in the de Bruijn graph of the reads. We recommend you start with k=+-30.
The order of the parameters is not fixed.

Other parameters are:
- ```--omit_control```: logical switch. If you use this switch instead of specifying a control file, the program will not perform any control normalization.
- ```--minimal_abundance <value>```: absolute number denoting the minimal abundance of a k-mer. 
- ```--minimal_abundance_percentile <value>```: percentile of k-mer abundancy to use as minimal abundance threshold. If absolute value is defined, this parameter is not used. Default: 75.
- ```--working_directory <path>```: name of the directory with results. Default: "output_paperfly_datetime".
- ```--draw```: logical switch. If used, the layout of the weakly connected components is drawn. There is a timeout threshold, therefore, too big components may not be drawn.
- ```--exclude_low```: logical switch. If used, the insufficiently abundant k-mers will be completely deleted. Otherwise, these k-mers are included in the clustering alignment step.
- ```--allowed_misses <value>```: distance threshold for clustering relative to k ("how many errors can be in a k-mer for two sequences to be mapped together"). Default: 1. On human data, we got much better results with setting 0.01.
- ```--window <value>```: rolling window width for curve smoothing during peak finding. Default: 100.
- ```--prominence <value>```: minimal prominence of a peak. Default: 20.
- ```--peak_format <consensus|sq_count|meme>```: whether to print only a consensus sequence of a peak (default), to print all the sequences with their count, or to print a single sequence for every occurence (recommended if you plan to follow up with MEME).
- ```--peak_min_width <value>```: minimal peak width. Default: k/2.
- ```--peak_max_width <value>```: maximal peak width. Default: infinity.

# Reference
The article describing the method was not published yet, but we are working on that. Both BCALM and Jellyfish will be properly cited in the publication.
