"""Enumerated Labels and custom exception and errors."""
from __future__ import annotations
from typing import Protocol
from enum import Enum
from Bio.PDB.Residue import Residue
from Bio.SeqUtils import seq1


class ChiralCenters(Enum):
    """Enumeration of L/D Chirality"""
    L = 0
    D = 2


class AmideBonds(Enum):
    """Enumeration to handle trans, cis stereo isomerism"""
    TRANS = 0
    CIS_PROLINE = 1
    NON_PLANAR = 2
    CIS = 3


class Clashes(Enum):
    """Enumeration for clashes."""
    VDW = 1


class StructuralIrregularity(Protocol):
    """Irregularity interface."""
    score: float
    def to_pymol(self):
        ...

    def to_cli(self):
        ...


class CoupleIrregularity:
    """handles irregularities that depend on two residues."""
    def __init__(self, res_a: Residue, res_b: Residue, score) -> CoupleIrregularity:
        self.res_a = SingleIrregularity(res_a, 0)
        self.res_b = SingleIrregularity(res_b, 0)
        self.score = score

    def to_pymol(self) -> str:
        return f'{self.res_a.to_pymol()} or {self.res_b.to_pymol()}'

    def to_cli(self) -> str:
        return f'{self.res_a.to_cli()}-{self.res_b.to_cli()}'


class SingleIrregularity:
    """handles irregularities that only depend on one residue."""
    def __init__(self, residue: Residue, score) -> StructuralIrregularity:
        self.code = seq1(residue.get_resname())
        self.number = residue.get_id()[1]
        self.score = score

    def to_pymol(self) -> str:
        """convert to pymol selectable string."""
        return f'resid {self.number}'

    def to_cli(self) -> str:
        """convert to displayable string in CLI."""
        return f'{self.code}{self.number:03}'
