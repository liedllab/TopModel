class ProlineException(Exception):
    """Prolines are more likely to be in cis-conformation, especially if they are preceded by
    Glycine or an aromatic residue."""


class GlycineException(Exception):
    """Glycine has no chiral centre"""


class PDBError(Exception):
    """Raised when needed information is missing in Structure"""
