# PDBear

A small python script to check PDB structures.

## Installation

```bash
pip install git+https://github.com/liedllab/PDBear.git
```

## Test Installation

```bash
pdbear --help
```

# Usage

PDBear is a small script that checks the chiralities and the amide bonds of the amino acids in a 
given PDB structure. Optionally, the PDBear can open the structure in PyMOL. The problems that were
found in the PDB structure by the PDBear are highlighted accordingly. For this pymol needs to be in
the PATH.

```bash
pdbear -f path/to/file.pdb
```

# Dependencies

- Biopython
- Click
- PyMOL (to open the structure in PyMOL)
