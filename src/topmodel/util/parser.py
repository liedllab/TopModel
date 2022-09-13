"""Handle the loading of the Structure from User input."""
from __future__ import annotations
from pathlib import Path
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
        ...


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


def validate_pdb_code(code: str) -> None:
    """check if the code is a valid PDB code.

:param code: str
:raises: PDBCodeError
    if the provided code does not conform to the required format."""
    code_length = 4
    if len(code) != code_length:
        raise PDBCodeError(f"{code} is does not match the required length of {code_length}")
    if not code.isalnum():
        raise PDBCodeError(f"Invalid characters detected in {code}. Only 0-9 and A-Z is allowed")
    if not code[0].isdecimal():
        raise PDBCodeError("First character is required to be numeric.")


def from_database(code: Path | str, directory: Path | str) -> Path:
    """Retrieve a code from the PDB database.

:param code: str
    PDB code
:param directory: str
    path to directory where the file should be downloaded to.

:returns: str
    path to downloaded file

:raises: KeyError if no PDB was downloaded."""
    code = str(code)
    directory = Path(directory)

    with BlockSTDOUT():
        database = PDBList()
        database.retrieve_pdb_file(code, pdir=str(directory))
    try:
        file = next(directory.glob(f'*{code}*'))
    except StopIteration:
        raise KeyError(f'Desired structure \"{code}\" does not exist') from None
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
            validate_pdb_code(str(filename))
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
