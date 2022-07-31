from pathlib import Path
import warnings
import yaml

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from hypothesis import given, assume
from hypothesis.strategies import integers, sampled_from, composite

from pdbear import chirality


BASE = Path('./data/modelled')
with open(BASE / "config.yaml") as f:
    CONFIG = yaml.safe_load(f)

parser = PDBParser()


@composite
def example_data(draw, examples: dict):
    key = draw(sampled_from(list(examples)))
    res = draw(integers(min_value=examples[key]['min_res'], max_value=examples[key]['max_res']))

    return key, examples[key]['pdb'], res


@given(example_data(CONFIG))
def test_chirality_alanine_dipeptide(data):
    key, path, res_number = data
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=PDBConstructionWarning)
        pdb = parser.get_structure("", BASE / path)

    residue = list(pdb.get_residues())[res_number]
    
    assume(residue.get_resname() != 'GLY')

    chiral_info = chirality.assign_chirality_amino_acid(residue) 
    if chiral_info == 'L':
        assert res_number not in CONFIG[key]['d_amino']
    else:
        assert res_number in CONFIG[key]['d_amino']

