from pathlib import Path
import warnings
import yaml

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from hypothesis import given, example, settings, strategies as st

from pdbear import amide_bond

parser = PDBParser()
BASE = Path('./data/alanine_dipeptide')
with open(BASE / "config.yaml") as f:
    CONFIG = yaml.safe_load(f)

@given(st.integers(min_value=1, max_value=8))
def test_amide_alanine_dipeptide(n: int):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=PDBConstructionWarning)
        pdb = parser.get_structure(n, BASE / f"{n}.pdb")

    stereo = amide_bond.get_stereo(pdb)
    target = CONFIG[n]["amide"]
    assert stereo == target

if __name__ == '__main__':
    test_alanine()
