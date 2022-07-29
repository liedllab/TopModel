from pathlib import Path
import warnings
import yaml

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from hypothesis import given, example, settings, strategies as st

from pdbear import chirality


BASE = Path('./data/alanine_dipeptide')
with open(BASE / "config.yaml") as f:
    CONFIG = yaml.safe_load(f)

parser = PDBParser()


@given(
        n = st.integers(min_value=1, max_value=8),
        res_number = st.integers(min_value=0, max_value=CONFIG["n_res"]),
        )
def test_chirality_alanine_dipeptide(n: int, res_number: int):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=PDBConstructionWarning)
        pdb = parser.get_structure(n, BASE / f"{n}.pdb")

    residue = list(pdb.get_residues())[res_number]

    assert chirality.assign_chirality(residue) == CONFIG[n]["chirality"][res_number]

if __name__ == '__main__':
    test_alanine()
