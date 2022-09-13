"""Enumerated Labels and custom exception and errors."""
from __future__ import annotations
import os
import sys
from typing import Protocol, runtime_checkable
from enum import Enum
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.SeqUtils import seq1


class ChiralCenters(Enum):
    """Enumeration of L/D Chirality"""
    L = 0
    D = 10


class AmideBonds(Enum):
    """Enumeration to handle trans, cis stereo isomerism"""
    TRANS = 0
    CIS_PROLINE = 1
    NON_PLANAR = 2
    CIS = 10


class Clashes(Enum):
    """Enumeration for clashes."""
    VDW = 1


@runtime_checkable
class StructuralIrregularity(Protocol):
    """Irregularity interface."""
    score: float

    def to_pymol(self):
        ...

    def to_cli(self):
        ...


class CoupleIrregularity:
    """Handles irregularities that depend on two residues."""
    def __init__(self, res_a: Residue, res_b: Residue, score):
        self.res_a = SingleIrregularity(res_a, 0)
        self.res_b = SingleIrregularity(res_b, 0)
        self.score = score

    def to_pymol(self) -> str:
        return f'{self.res_a.to_pymol()} or {self.res_b.to_pymol()}'

    def to_cli(self) -> str:
        return f'{self.res_a.to_cli()}-{self.res_b.to_cli()}'

    def __repr__(self):
        cls = self.__class__
        return f'{cls.__name__}({self.res_a}, {self.res_b}, {self.score!r})'


class SingleIrregularity:
    """Handles irregularities that only depend on one residue."""
    def __init__(self, residue: Residue, score):
        self.code = seq1(residue.get_resname())
        self.number = residue.get_id()[1]
        self.score = score

    def to_pymol(self) -> str:
        """Convert to pymol selectable string."""
        return f'resid {self.number}'

    def to_cli(self) -> str:
        """Convert to displayable string in CLI."""
        return f'{self.code}{self.number:03}'

    def __repr__(self):
        cls = self.__class__
        return f'{cls.__name__}({self.code!r}, {self.number!r}, {self.score!r})'

    def __str__(self):
        cls = self.__class__
        return f'{cls.__name__}({self.code!r}, {self.number!r})'


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


class SelectiveStructure(Structure):
    """Wrap Structure functionality with selective output of residues."""
    def get_residues(self):
        """return only non heteroatoms."""
        for chain in self.get_chains():
            for residue in chain:
                if residue.resname in {'NME', 'HOH', 'ACE', 'WAT'}:
                    continue
                hetero, *_ = residue.get_id()
                if hetero.strip() != '':
                    continue
                if residue.resname in {'NME', 'HOH', 'ACE', 'WAT'}:
                    continue
                yield residue
