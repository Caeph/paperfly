# Paperfly: a software for analysis of ChIP-seq (or similar) experiments without the reference genome

# Requirements:
- python 3.8 or newer
- jellyfish -- can be downloaded from repositories
- bcalm -- can be downloaded from repositories
- mono-complete for building the c# 
- graphviz for drawing

All of these are available from the Linux repositories.

# Instalation:
Clone the repository. Use the included makefile to 
- create a python virtual environment and install the required packages
- build the C# binaries

Just run 
```
make
```
and then update your path if needed.
```
export PATH=$PATH:<paperfly directory>
```

# Reference
The article describing hte method was not published yet.
