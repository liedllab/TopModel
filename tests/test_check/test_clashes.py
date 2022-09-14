"""only small test set to check the overall return format"""
from itertools import chain
from enum import Enum
from topmodel.util.utils import StructuralIrregularity
from topmodel.util.parser import get_structure
from topmodel.check import clashes


structure = get_structure('./data/clashes/tryptophan_helix.pdb')
found = clashes.get_clashes(structure)


def test_type():
    assert isinstance(found, dict)


def test_keys():
    assert all(isinstance(key, Enum) for key in found)


def test_values():
    assert all(isinstance(value, StructuralIrregularity)
               for value in chain(*found.values()))
