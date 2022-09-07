# PDBear

A python script to check PDB structures.

## Installation

```bash
pip install git+https://github.com/liedllab/PDBear.git
```

## Test Installation

```bash
pdbear --help
```

# Command-Line Interface (CLI)
```bash
pdbear path/to/file(s).pdb
```

PDBear is a small script that checks the chiralities and the amide bonds of the amino acids in a 
given PDB structure. Optionally, the PDBear can open the structure in PyMOL. The problems that were
found in the PDB structure by the PDBear are highlighted accordingly. For this PyMOL needs to be in
the PATH.

## Chirality

The PDBear can assign `L` and `D` to aminoacids. In the command-line interface (CLI) only the `D`
aminoacids are printed.

## Amide bonds

The PDBear can assign `CIS` and `TRANS` to the amide bonds depending on the dihedral angle between
CA-C-N-CA. Amide bonds that could not be assigned to either `CIS` or `TRANS` are labeled as
`NON_PLANAR`.
Cis amide bonds to prolines are labelled separately as they occur more frequently.

## Clashes

The PDBear detects Van der Waals clashes once the distance between any atom of the sidechain is
closer to any atom of any other residue than their combined Van der Waals radii.

# As a package

PDBear can be imported as a package. The package contains the small modules `chirality`,
`amide_bond` and `clashes` that provide functions to calculate the respective property.

```python
import pdbear # or
from pdbear import chirality, amide_bond
```

# Dependencies

- biopython
- Click
- colorama
- PyMOL (to open the structure in PyMOL)
