# Fplib
Fplib is is a library for calculating fingerprint and measuring similarites of materials (e.g. crystals, clusters, and molecules). The library is written in C, with python interface. 

The detailed algorithm of fplib is describted the following paper:

- "[A fingerprint based metric for measuring similarities of crystaline structures](http://scitation.aip.org/content/aip/journal/jcp/144/3/10.1063/1.4940026)",
  Li Zhu, Maximilian Amsler, Tobias Fuhrer, Bastian Schaefer, Somayeh Fareji, Alireza Ghasemi, Migle Grauzinyte, Chris Wolverton, and Stefan Goedecker
  **J. Chem. Phys. 144**, 034203 (2016)



# Installation

To install python-fplib using `setup.py`, python header files (python-dev), C-compiler (e.g., gcc, clang), numpy, and lapack (openblas, atlas, Intel-MKL) are required before the build. The installation steps are shown as follows:

1. Go to the fppy directory

2. Type the command:

   `% pip install .

   

# Example

Examples are found in `examples` directory.

