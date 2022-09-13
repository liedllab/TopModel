from itertools import chain
import warnings
import pytest
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from topmodel.check import amide_bond
from topmodel.util.utils import AmideBonds, StructuralIrregularity
from topmodel.util.errors import MissingInformationError
from topmodel.util.parser import get_structure as special_get_structure


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


class Test7SG5:
    file = './data/fileformats/7SG5_model.pdb'

    def test_7SG5_failure(self):
        struc = get_structure(self.file)
        *_, a, b = struc.get_residues()
        with pytest.raises(MissingInformationError):
            amide_bond.assign_stereo(a, b)

    def test_7SG5_and_solution(self):
        struc = special_get_structure(self.file)
        *_, a, b = struc.get_residues()
        amide_bond.assign_stereo(a, b)

    def test_7SG5_cis_amides(self):
        struc = special_get_structure(self.file)
        stereo = amide_bond.get_amide_stereo(struc)
        assert len(stereo[AmideBonds.CIS]) == 3
