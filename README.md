# libfp

libfp is a library for calculating crystalline fingerprints and measuring similarities of materials (e.g., crystals, clusters, and molecules). The library is written in C with a Python interface.

Research papers where this code has been used include:

- [EOSnet: Embedded Overlap Structures for Graph Neural Networks in Predicting Material Properties, **J. Phys. Chem. Lett. 16**, 717 (2025)](https://pubs.acs.org/doi/full/10.1021/acs.jpclett.4c03179)
- [Accelerating Structural Optimization through Fingerprinting Space Integration on the Potential Energy Surface, **J. Phys. Chem. Lett. 15**, 3185 (2024)](https://pubs.acs.org/doi/10.1021/acs.jpclett.4c00275)
- [Quantum structural fluxion in superconducting lanthanum polyhydride, **Nature Commun. 14**, 1674 (2023)](https://www.nature.com/articles/s41467-023-37295-1)
- [Phase Transition Pathway Sampling via Swarm Intelligence and Graph Theory, **J. Phys. Chem. Lett. 10**, 5019 (2019)](https://pubs.acs.org/doi/10.1021/acs.jpclett.9b01715)



## Installation

### Prerequisites

Before installing libfp, ensure you have the following:

- Python header files (python-dev)
- C compiler (e.g., gcc, clang)
- NumPy
- OpenBLAS, LAPACK, or MKL

### Installation Steps

You can install libfp using one of the following methods:

#### Option 1: Install from PyPI

To install the latest stable version from PyPI, simply run:

```
pip install libfp
```

#### Option 2: Install from source

1. Clone the repository:
   ```
   git clone https://github.com/Rutgers-ZRG/libfp.git
   cd libfp
   ```

2. Install using pip:
   ```
   pip install .
   ```

## Usage

To use libfp in your Python project:

```python
import libfp

# Your code here
```

## Examples

Example
Examples are found in `examples` directory.

## Algorithm

The detailed algorithm of libfp is described in the following paper:

- "[A fingerprint based metric for measuring similarities of crystalline structures](http://scitation.aip.org/content/aip/journal/jcp/144/3/10.1063/1.4940026)",
  Li Zhu, Maximilian Amsler, Tobias Fuhrer, Bastian Schaefer, Somayeh Fareji, Alireza Ghasemi, Migle Grauzinyte, Chris Wolverton, and Stefan Goedecker
  **J. Chem. Phys. 144**, 034203 (2016)

## Acknowledgements
This library is derived from an original Fortran version, which I wrote during my postdoctoral work in Prof. Goedecker's research group at the University of Basel. I am grateful to Prof. Goedecker for his guidance during my postdoc period. I also acknowledge support from Rutgers startup funding and the National Science Foundation Division of Materials Research (NSF-DMR), under Grant No. 2226700.