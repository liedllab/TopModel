from itertools import chain
import warnings
import pytest
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from topmodel.check import chirality
from topmodel.util.utils import ChiralCenters, StructuralIrregularity
from topmodel.util.errors import MissingInformationError
from topmodel.util.parser import get_structure as get_clean_structure


parser = PDBParser()


def get_plain_structure(file):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=PDBConstructionWarning)
        struc = parser.get_structure('', file)
    return struc


@pytest.mark.parametrize("pdb_file,expected", [
    ('./data/alanine_dipeptide/1.pdb', (ChiralCenters.L, ChiralCenters.L)),
    ('./data/alanine_dipeptide/2.pdb', (ChiralCenters.L, ChiralCenters.D)),
    ('./data/alanine_dipeptide/3.pdb', (ChiralCenters.D, ChiralCenters.L)),
    ('./data/alanine_dipeptide/4.pdb', (ChiralCenters.D, ChiralCenters.D)),
    ('./data/alanine_dipeptide/5.pdb', (ChiralCenters.L, ChiralCenters.L)),
    ('./data/alanine_dipeptide/6.pdb', (ChiralCenters.L, ChiralCenters.D)),
    ('./data/alanine_dipeptide/7.pdb', (ChiralCenters.D, ChiralCenters.L)),
    ('./data/alanine_dipeptide/8.pdb', (ChiralCenters.D, ChiralCenters.D)),
    ])
@pytest.mark.parametrize("get_structure", [get_plain_structure, get_clean_structure])
def test_assign_alanine_dipeptide(pdb_file, expected, get_structure):
    """test alanine dipeptide."""
    struc = get_structure(pdb_file)
    a, b = struc.get_residues()
    chiral_info = chirality.assign_chirality_amino_acid(a), chirality.assign_chirality_amino_acid(b)
    assert chiral_info == expected


@pytest.mark.parametrize("get_structure", [get_plain_structure, get_clean_structure])
class TestSGMT:
    file = './data/modelled/sg_imgt_protonated.pdb'

    def test_output_type(self, get_structure):
        struc = get_structure(self.file)
        chiral = chirality.get_chirality(struc)
        assert isinstance(chiral, dict)

    def test_output_keys(self, get_structure):
        struc = get_structure(self.file)
        chiral = chirality.get_chirality(struc)
        assert all(isinstance(key, ChiralCenters) for key in chiral)

    def test_output_values(self, get_structure):
        struc = get_structure(self.file)
        chiral = chirality.get_chirality(struc)
        assert all(isinstance(val, StructuralIrregularity) for val in chain(*chiral.values()))

    def test_output_length(self, get_structure):
        struc = get_structure(self.file)
        chiral = chirality.get_chirality(struc)
        len_values = sum(len(x) for x in chiral.values())
        len_structure = len(list(struc.get_residues()))
        assert len_values == len_structure

    def test_found_d_aminoacids_length(self, get_structure):
        """right answer was checked using moe and cpptraj."""
        struc = get_structure(self.file)
        chiral = chirality.get_chirality(struc)
        assert len(chiral[ChiralCenters.D]) == 3

    def test_found_d_aminacids(self, get_structure):
        """correct answer was checked with moe."""
        # cpptraj assigns 3 different residues?
        struc = get_structure(self.file)
        chiral = chirality.get_chirality(struc)
        found = {109, 112, 115}
        expected = {calculated.number for calculated in chiral[ChiralCenters.D]}
        assert found == expected


class Test7SG5:
    file = './data/fileformats/7SG5_model_complex.pdb'
    struc = get_clean_structure(file)
    chiral = chirality.get_chirality(struc)

    def test_found_d_aminoacids_length(self):
        """right answer was checked using moe and cpptraj."""
        assert len(self.chiral[ChiralCenters.D]) == 4

    def test_found_d_aminacids(self):
        """correct answer was checked with moe and cpptraj."""
        found = {216, 218, 221, 222}
        expected = {calculated.number for calculated in self.chiral[ChiralCenters.D]}
        assert found == expected

    def test_error_nme(self):
        """CA needed to determine the chirality."""
        struc = get_plain_structure(self.file)
        with pytest.raises(MissingInformationError):
            chirality.get_chirality(struc)


# old implementation raised an error, which now no longer works
@pytest.mark.xfail
@pytest.mark.parametrize("get_structure", [get_plain_structure, get_clean_structure])
def test_failure_missing_ha(get_structure):
    struc = get_structure('./data/modelled/sg_imgt.pdb')
    with pytest.raises(MissingInformationError):
        chirality.get_chirality(struc)
