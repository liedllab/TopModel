"""Handle the loading of the Structure from User input."""
from __future__ import annotations
import sys
import os
from pathlib import Path
import tempfile
from typing import Protocol
import warnings
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmtf import MMTFParser
from Bio.PDB import PDBParser, PDBList, PDBExceptions
from Bio.PDB.Structure import Structure

from .errors import PDBCodeError


class Parser(Protocol):
    """Parser Protocol which the other functions depend on."""
    def get_structure(self, identifier: str, filename: str) -> Structure:
        """get structure from filename."""
        ...


def get_parser(filename: str) -> Parser:
    """return appropriate parser depending on filetype.
:param filename: str

:returns: Parser

:raises: ValueError
    if no appropriate parser could be assigned."""

    match filename.split('.')[-1]:
        case 'mmcif' | 'cif':
            parser = MMCIFParser()
        case 'mmtf':
            parser = MMTFParser()
        case 'pdb':
            parser = PDBParser()
        case ending:
            raise ValueError(f"No parser available for {ending} format")
    return parser


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


class BlockSTDOUT:
    """Context manager to block all stdout and stderr calls.
Usefull for library methods that print information instead of using warnings."""
    def __init__(self):
        self._original_stdout = sys.stdout
        self._original_stderr = sys.stderr

    def __enter__(self):
        sys.stdout = open(os.devnull, 'w', encoding='utf8')
        sys.stderr = open(os.devnull, 'w', encoding='utf8')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stderr.close()
        sys.stdout = self._original_stdout
        sys.stderr = self._original_stderr


def from_database(code: str, directory: str) -> str:
    """Retrieve a code from the PDB database.

:param code: str
    PDB code
:param directory: str
    path to directory where the file should be downloaded to.

:returns: str
    path to downloaded file

:raises: KeyError if no PDB was downloaded."""

    with BlockSTDOUT():
        database = PDBList()
        database.retrieve_pdb_file(code, pdir=directory)
        directory = Path(directory)
    try:
        file = next(directory.glob(f'*{code}*'))
    except StopIteration:
        raise KeyError(f'Desired structure \"{code}\" does not exist') from None
    return str(file)


def struc_from_file(filename: str, parser: Parser) -> Structure:
    """Construct a Structure given the file and a parser.

:param filename: str
    path to file
:param parser: Parser
    Parser that can construct the Structure object.

:returns: Structure"""

    stem = Path(filename).stem
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=PDBExceptions.PDBConstructionWarning)
        structure = parser.get_structure(stem, filename)
    return structure


def get_structure(filename: str) -> Structure:
    """Get Structure from filename or PDB Code. Downloads the files from the PDB into a temporary
directory.

:param filename: str
    path to file or PDB code.

:returns: Structure"""

    with tempfile.TemporaryDirectory() as directory:
        try:
            parser: Parser = get_parser(filename)
        except ValueError:
            validate_pdb_code(filename)
            filename: str = from_database(filename, directory)
            parser: Parser = get_parser(filename)

        structure: Structure = struc_from_file(filename, parser)
    return structure

def main():
    pass

if __name__ == '__main__':
    raise SystemExit(main())
