# fplib

Fplib is a library for calculating fingerprints and measuring similarities of materials (e.g., crystals, clusters, and molecules). The library is written in C with a Python interface.


## Installation

### Prerequisites

Before installing python-fplib, ensure you have the following:

- Python header files (python-dev)
- C compiler (e.g., gcc, clang)
- NumPy
- OpenBLAS

### Installation Steps

1. Clone the repository:
   ```
   git clone https://github.com/Rutgers-ZRG/fplib.git
   cd fplib
   ```

2. Install using pip:
   ```
   pip install .
   ```

## Usage

To use fplib in your Python project:

```python
import fplib

# Your code here
```

## Examples

Example
Examples are found in `examples` directory.

## Algorithm

The detailed algorithm of fplib is described in the following paper:

- "[A fingerprint based metric for measuring similarities of crystalline structures](http://scitation.aip.org/content/aip/journal/jcp/144/3/10.1063/1.4940026)",
  Li Zhu, Maximilian Amsler, Tobias Fuhrer, Bastian Schaefer, Somayeh Fareji, Alireza Ghasemi, Migle Grauzinyte, Chris Wolverton, and Stefan Goedecker
  **J. Chem. Phys. 144**, 034203 (2016)