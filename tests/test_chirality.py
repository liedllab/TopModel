from hypothesis import given, example, settings, strategies as st

from pdbear import chirality
from pathlib import Path
import yaml
from Bio.PDB.PDBParser import PDBParser

parser = PDBParser()
BASE = Path('./data/alanine_dipeptide')
with open(BASE / "config.yaml") as f:
    CONFIG = yaml.safe_load(f)

@given(st.integers(min_value=1, max_value=8))
def test_alanine(n: int):
    pdb = parser.get_structure(n, BASE / f"{n}.pdb")
    chiralities = chirality.get_chiralities(pdb)
    target = CONFIG[n]["stereo"]
    print(target)
    assert chiralities == target

if __name__ == '__main__':
    test_alanine()
