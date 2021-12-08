# Paperfly: a software for analysis of ChIP-seq (or similar) experiments without the reference genome

The PseudoAssembly based Peak Finder (PAPerFly) assembles the sequencing reads seen during a ChIP-seq or similar experiment. 
This repository is still under construction and some additions will be made to the code.

# Requirements
- python 3.8 or newer, python3.8-venv 
- jellyfish -- can be downloaded from repositories or compiled from https://github.com/gmarcais/Jellyfish
- bcalm -- can be downloaded from repositories or compiled from https://github.com/GATB/bcalm
- mono, msbuild and nuget for building and running the c# applications (msbuild can usually be downloaded as a part mono-complete packages)
- graphviz for drawing (only sfdp is used)

All of these are available from the Linux repositories. All of these should be in your path.

# Instalation
Clone the repository. Use the included makefile to 
- create a python virtual environment and install the required packages
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
Still to do. A compplete list of required and optional parameters will be added as soon as possible, as well as some example files.

# Reference
The article describing the method was not published yet, but we are working on that. Both BCALM and Jellyfish will be properly cited in the publication.
