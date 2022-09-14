from itertools import chain
import warnings
import pytest
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from topmodel.check import amide_bond
from topmodel.util.utils import AmideBonds, StructuralIrregularity
from topmodel.util.errors import MissingInformationError
from topmodel.util.parser import get_structure as get_clean_structure


parser = PDBParser()


def get_plain_structure(file):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=PDBConstructionWarning)
        struc = parser.get_structure('', file)
    return struc


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
@pytest.mark.parametrize("get_structure", [get_plain_structure, get_clean_structure])
def test_assign_alanine_dipeptide(pdb_file, expected, get_structure):
    """test alanine dipeptide."""
    struc = get_structure(pdb_file)
    a, b = struc.get_residues()
    assert amide_bond.assign_stereo(a, b) is expected


@pytest.mark.parametrize("get_structure", [get_plain_structure, get_clean_structure])
class TestGet:
    file = './data/modelled/output.pdb'

    def test_output_type(self, get_structure):
        struc = get_structure(self.file)
        stereo = amide_bond.get_amide_stereo(struc)
        """make sure output format is correct."""
        assert isinstance(stereo, dict)

    def test_output_keys(self, get_structure):
        struc = get_structure(self.file)
        stereo = amide_bond.get_amide_stereo(struc)
        assert all(isinstance(x, AmideBonds) for x in stereo)

    def test_output_values(self, get_structure):
        struc = get_structure(self.file)
        stereo = amide_bond.get_amide_stereo(struc)
        assert all(isinstance(x, StructuralIrregularity) for x in chain(*stereo.values()))

    def test_get_output_length(self, get_structure):
        struc = get_structure(self.file)
        stereo = amide_bond.get_amide_stereo(struc)
        len_values = sum(len(x) for x in stereo.values())
        len_structure = len(list(struc.get_residues()))
        # 1 less bond than residues
        assert len_values == (len_structure - 1)


class Test7SG5:
    file = './data/fileformats/7SG5_model.pdb'

    def test_7SG5_failure(self):
        struc = get_plain_structure(self.file)
        *_, a, b = struc.get_residues()
        with pytest.raises(MissingInformationError):
            amide_bond.assign_stereo(a, b)

    def test_7SG5_and_solution(self):
        struc = get_clean_structure(self.file)
        *_, a, b = struc.get_residues()
        amide_bond.assign_stereo(a, b)

    def test_7SG5_number_cis_amides(self):
        struc = get_clean_structure(self.file)
        stereo = amide_bond.get_amide_stereo(struc)
        assert len(stereo[AmideBonds.CIS]) == 3
