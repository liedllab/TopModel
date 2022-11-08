# TopModel

A python script to quickly inspect and highlight issues in structure models.

![Tests](https://github.com/liedllab/TopModel/actions/workflows/tests.yml/badge.svg)

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

TopModel is a command line tool that checks the chiralities, the amide bonds and overall VDW clashes
in a structure model. Optionally, the structure can be opened in PyMOL to visualize 
the issues found. For this PyMOL needs to be in PATH.

## Chirality

TopModel can assign `L` and `D` to aminoacids. In the command-line interface (CLI) only the `D`
aminoacids are shown.
The chirality is computed by computing a normal vector $\vec{n}$ to the plane defined by the triangle 
of the three highest priority atoms $A$ around the chiral center where the suffix denotes the
priority. 

$$\vec{n} =  \overrightarrow{A_3A_1} \times \overrightarrow{A_3A_2}$$

The direction of the normal vector is determined by the orientation of the side chains. By 
calculating the dot product of the normal vector and a vector from the chiral center to the plane the 
relative orientation of the three atoms around the chiral center, and thus the chirality, can be 
determined from the sign of the result.

$$
chirality = \vec{n} \cdot \overrightarrow{A_3A_{center}} = 
        \begin{cases} 
        \mathbf{D}\, \text{if } > 0 \\
        \mathbf{L}\, \text{if } < 0 
        \end{cases}
$$

## Amide bonds

TopModel can assign `CIS` and `TRANS` to the amide bonds depending on the dihedral angle defined by
$C_{\alpha}CNC_\alpha$.
Amide bonds that could not be assigned to either `CIS` or `TRANS` are labeled as `NON_PLANAR`.

Cis amide bonds to prolines are labelled separately as they occur more frequently.

## Clashes

TopModel detects Van der Waals clashes by calculating the distance between all pair of atoms that
are within 5 Å of each other. A clash is defined by:
$$d_{AB} < r_A + r_B - 0.5Å$$

# As a package

TopModel can be imported as a package. The package contains the small modules `chirality`,
`amide_bond` and `clashes` that provide functions to calculate the respective property.

```python
import topmodel # or
from topmodel.check import chirality, amide_bond, clashes
```

# Dependencies

- biopython
- scipy
- mendeleev
- Click
- colorama
- PyMOL (to open the structure in PyMOL)
