"""Handle the loading of the Structure from User input."""
from __future__ import annotations
from pathlib import Path
import re
import tempfile
from typing import Protocol, runtime_checkable
import warnings
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBParser, PDBList, PDBExceptions
from Bio.PDB.Structure import Structure

from .errors import PDBCodeError
from .utils import BlockSTDOUT, SelectiveStructure
from .MMTFParser import TopModelMMTFParser


@runtime_checkable
class Parser(Protocol):
    """Parser Protocol which the other functions depend on."""
    def get_structure(self, identifier: str, filename: str) -> Structure:
        """get structure from filename."""


def get_parser(filename: Path | str) -> Parser:
    """return appropriate parser depending on filetype.
:param filename: str

:returns: Parser

:raises: ValueError
    if no appropriate parser could be assigned."""
    filename = str(filename)
    suffix = filename.rsplit('.', maxsplit=1)[-1].lower()
    parser: Parser
    if suffix in {'mmcif', 'cif'}:
        parser = MMCIFParser()
    elif suffix == 'mmtf':
        parser = TopModelMMTFParser()
    elif suffix == 'pdb':
        parser = PDBParser()
    else:
        raise ValueError(f'No parser available for {suffix} format')
    return parser
#    match filename.rsplit('.', maxsplit=1)[-1]:
#        case 'mmcif' | 'cif':
#            parser = MMCIFParser()
#        case 'mmtf':
#            parser = MMTFParser()
#        case 'pdb':
#            parser = PDBParser()
#        case ending:
#            raise ValueError(f"No parser available for {ending} format")


class PDBCode(str):
    _PATTERN = re.compile('[0-9]{1}[0-9a-zA-Z]{3}')

    @property
    def valid(self) -> bool:
        return bool(self._PATTERN.fullmatch(self))


def from_database(code: Path | str, directory: Path | str) -> Path:
    """Retrieve a code from the PDB database.

:param code: str
    PDB code
:param directory: str
    path to directory where the file should be downloaded to.

:returns: str
    path to downloaded file

:raises: KeyError if no PDB was downloaded."""
    code = PDBCode(code)
    if not code.valid:
        raise PDBCodeError('Invalid PDB access code')
    directory = Path(directory)

    with BlockSTDOUT():
        database = PDBList()
        database.retrieve_pdb_file(code, pdir=str(directory))
    try:
        file = next(directory.glob('*'))
    except StopIteration:
        raise KeyError(f'Desired structure \"{code}\" does not exist or could not be retrieved')\
                from None
    return file


def struc_from_file(filename: Path | str, parser: Parser) -> Structure:
    """Construct a Structure given the file and a parser.

:param filename: str
    path to file
:param parser: Parser
    Parser that can construct the Structure object.

:returns: Structure"""
    stem = Path(filename).stem
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=PDBExceptions.PDBConstructionWarning)
        structure = parser.get_structure(stem, str(filename))
    return structure


def get_structure(filename: Path | str) -> Structure:
    """Get Structure from filename or PDB Code. Downloads the files from the PDB into a temporary
directory.

:param filename: str
    path to file or PDB code.

:returns: Structure"""
    with tempfile.TemporaryDirectory() as directory:
        try:
            parser: Parser = get_parser(filename)
        except ValueError:
            filename = from_database(filename, directory)
            parser = get_parser(filename)

        structure: Structure = struc_from_file(filename, parser)
        # replace get_residues() function
        structure.__class__ = SelectiveStructure
    return structure


def main():
    """Main."""
    pass


if __name__ == '__main__':
    raise SystemExit(main())
