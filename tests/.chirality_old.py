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


@given(st.integers(min_value=1, max_value=8))
def test_chirality_alanine_dipeptide(n: int):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=PDBConstructionWarning)
        pdb = parser.get_structure(n, BASE / f"{n}.pdb")

    chiralities = chirality.get_chiralities(pdb)
    target = CONFIG[n]["stereo"]
    assert [chirality for _, chirality in chiralities.values()] == target

if __name__ == '__main__':
    test_alanine()
