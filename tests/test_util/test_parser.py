import pytest
from hypothesis import given, strategies as st
from topmodel.util import parser as tm_parser
from topmodel.util.errors import PDBCodeError


@pytest.mark.parametrize('filename', [
    './data/fileformats/7lka.cif',
    './data/fileformats/4HHB.mmtf',
    './data/fileformats/7SG5_model.pdb'])
def test_parser_protocol(filename):
    parser = tm_parser.get_parser(filename)
    assert isinstance(parser, tm_parser.Parser)


def test_parser_failure():
    with pytest.raises(ValueError):
        tm_parser.get_parser('test.xyz')


@given(code=st.from_regex('[0-9]{1}[0-9a-zA-Z]{3}', fullmatch=True))
def test_valid_pdb_code(code):
    assert tm_parser.PDBCode(code).valid


def test_db_failure(tmpdir):
    with pytest.raises(KeyError):
        tm_parser.from_database('9zzz', tmpdir)


def test_db(tmpdir):
    tm_parser.from_database('1igy', tmpdir)


def test_get_structure_pdb_retrieve():
    struc = tm_parser.get_structure('7lka')
    struc2 = tm_parser.get_structure('./data/fileformats/7lka.cif')
    assert struc == struc2


def test_get_structure_failure():
    with pytest.raises(PDBCodeError):
        tm_parser.get_structure('test.xyz')
