from pathlib import Path
import warnings
import yaml
from hypothesis import given, assume
from hypothesis.strategies import integers, sampled_from, composite
import pytest
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from topmodel.check import chirality

BASE = Path('./data/modelled')
with open(BASE / "config.yaml") as f:
    CONFIG = yaml.safe_load(f)

parser = PDBParser()


@composite
def example_data(draw, examples: dict):
    key = draw(sampled_from(list(examples)))
    res = draw(integers(min_value=examples[key]['min_res'], max_value=examples[key]['max_res']))

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=PDBConstructionWarning)
        pdb = parser.get_structure(key, BASE / examples[key]['pdb'])

    residue = list(pdb.get_residues())[res]

    return key, residue, res


@given(example_data(CONFIG))
def test_chirality_modelled(data):
    key, residue, res_number  = data
    
    if residue.get_resname() == 'GLY':
        with pytest.raises(chirality.GlycineException):
            chirality.assign_chirality_amino_acid(residue)
    else:
        chiral_info = chirality.assign_chirality_amino_acid(residue) 
        if chiral_info == 'L':
            assert res_number not in CONFIG[key]['d_amino']
        else:
            assert res_number in CONFIG[key]['d_amino']


