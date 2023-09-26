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

These parameters are valid for the latest package. We are currently reworking parts of the software. Once ready, some changes to parameter setting will be done.

To run the program when it is in your path, use
```
paperfly --input_description <tsv file> --input_directory <path to directory with all input files> --k <k>
```
The input and the control file can be in FASTA, FASTQ or gzipped FASTQ format (then use the fastq option). The value of K corresponds to the lenght of k-mers in the de Bruijn graph of the reads. We recommend you start with k=+-30.
The order of the parameters is not fixed.

The required parameters are:
```--input_description <tsv file>```: path to a tsv file describing treatment vs. control relationships for replicates. 
It should look like this: 

```
fastq   control
<input file 1>    <respective control file>
<input file 2>    <respective control file>
```

```--input_directory <path to directory with ALL files in input_description>```: path to a directory with all input files. Cannot be current directory."
```--k <K>```: k-mer length for sliding window probing of the sequencing data. Cannot be a multiple of 4 (BCALM property).

Other parameters are:

```--working_dir <path>```: path to a directory to store results. Must not exist, otherwise exception is thrown.

```--minimal_abundance <N>```: minimal abundance of a k-mer. Should be lower for higher k-mers. Note that lower abundance numbers lead to higher data complexity and longer runtime.If not specified, it is set as 90th percentile in the (non-unique) k-mer counts. Percentile can
                        be adjusted by the --minimal_abundance_percentile parameter.

```--minimal_abundance_percentile <0-100>```: defined percentile of sufficiently abundant kmers. This calculation is time-consuming, but comes in handy if you don\'t know much about the input data size.
  BCALM also makes a number of temporary files during the calculation. These will be removed afterwards. Default option, 95th percentile.

```--minimal_abundance_mapping <N>```: abundance threshold of a low abundance kmer count to be considered a sequencing error. Default: 10.

```--assembly_report_step <N>```: number of steps after which time and graph state are reported during pseudoassembly. Default: 250.

```--draw```: option to draw the components graphs and graphs of assembled profile enrichments.

```--no_store_low```: option to throw away low abundance kmers.

```--miss_percentage <0-100>```: identity percentage for a low abundance kmer to be mapped to assembled profile. Values from 0 to 100, default 90"

# Reference
The article describing the method was not published yet, but we are working on that. Both BCALM and Jellyfish will be properly cited in the publication.
