from itertools import chain
import warnings
import pytest
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from topmodel.check import amide_bond
from topmodel.util.utils import AmideBonds, StructuralIrregularity


parser = PDBParser()

def get_structure(file):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=PDBConstructionWarning)
        struc = parser.get_structure('', file)
    return struc


class TestGet:
    struc = get_structure('./data/modelled/output.pdb')
    stereo = amide_bond.get_amide_stereo(struc)

    def test_output_type(self):
        """make sure output format is correct."""
        assert isinstance(self.stereo, dict)

    def test_output_keys(self):
        assert all(isinstance(x, AmideBonds) for x in self.stereo)

    def test_output_values(self):
        assert all(isinstance(x, StructuralIrregularity) for x in chain(*self.stereo.values()))

    def test_get_output_length(self):
        len_values = sum(len(x) for x in self.stereo.values())
        len_structure = len(list(self.struc.get_residues()))
        # 1 less bond than residues
        assert len_values == (len_structure - 1)


@pytest.mark.parametrize("pdb_file,expected", [
    ('./data/alanine_dipeptide/1.pdb', AmideBonds.TRANS),
    ('./data/alanine_dipeptide/2.pdb', AmideBonds.TRANS),
    ('./data/alanine_dipeptide/3.pdb', AmideBonds.TRANS),
    ('./data/alanine_dipeptide/4.pdb', AmideBonds.TRANS),
    ('./data/alanine_dipeptide/5.pdb', AmideBonds.CIS),
    ('./data/alanine_dipeptide/6.pdb', AmideBonds.CIS),
    ('./data/alanine_dipeptide/7.pdb', AmideBonds.CIS),
    ('./data/alanine_dipeptide/8.pdb', AmideBonds.CIS),
    ])
def test_assign_alanine_dipeptide(pdb_file, expected):
    """test alanine dipeptide."""
    struc = get_structure(pdb_file)
    a, b = struc.get_residues()
    assert amide_bond.assign_stereo(a, b) is expected
