# TopModel

A python script to check structure models.

## Installation

```bash
pip install git+https://github.com/liedllab/topmodel.git
```

## Test Installation

```bash
topmodel --help
```

# Command-Line Interface (CLI)
```bash
topmodel path/to/file*.pdb
```

TopModel is a command line tool that checks the chiralities, the amide bonds and overall clashes of 
the amino acids in a structure model. Optionally, the structure can be opened in PyMOL to visualize 
the issues found. For this PyMOL needs to be in PATH.

## Chirality

TopModel can assign `L` and `D` to aminoacids. In the command-line interface (CLI) only the `D`
aminoacids are shown.

## Amide bonds

TopModel can assign `CIS` and `TRANS` to the amide bonds depending on the dihedral angle between
CA-C-N-CA. Amide bonds that could not be assigned to either `CIS` or `TRANS` are labeled as
`NON_PLANAR`.
Cis amide bonds to prolines are labelled separately as they occur more frequently.

## Clashes

TopModel detects Van der Waals clashes once the distance between any atom of the sidechain is
closer to any other atom than their combined Van der Waals radii.

# As a package

TopModel can be imported as a package. The package contains the small modules `chirality`,
`amide_bond` and `clashes` that provide functions to calculate the respective property.

```python
import topmodel # or
from topmodel import chirality, amide_bond, clashes
```

# Dependencies

- biopython
- scipy
- mendeleev
- Click
- colorama
- PyMOL (to open the structure in PyMOL)
