"""Enumerated Labels and custom exception and errors."""

from enum import Enum, auto


class ChiralCenters(Enum):
    """Enumeration of L/D Chirality"""
    L = auto()
    D = auto()


class StereoIsomers(Enum):
    """Enumeration to handle trans, cis stereo isomerism"""
    TRANS = auto()
    CIS = auto()
    INBETWEEN = auto()
    CIS_PROLINE = auto()


class Clashes(Enum):
    VDW = auto()


class Color(Enum):
    RED = 'red'
    YELLOW = 'yellow'
    MAGENTA = 'magenta'
    BLUE = 'blue'


class ProlineException(Exception):
    """Prolines are more likely to be in cis-conformation, especially if they are preceded by
    Glycine or an aromatic residue."""


class GlycineException(Exception):
    """Glycine has no chiral centre"""


class PDBError(Exception):
    """Raised when needed information is missing in Structure"""
