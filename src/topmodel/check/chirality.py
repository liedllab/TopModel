"""Calculate chirality of an aminoacid."""

from __future__ import annotations

import warnings
from Bio.PDB import vectors, Residue, Structure
import numpy as np
import numpy.typing as npt
from topmodel.util.errors import GlycineException, MissingInformationError
from topmodel.util.utils import ChiralCenters, SingleIrregularity


def assign_chirality_amino_acid(residue: Residue.Residue) -> ChiralCenters:
    """Assign L/D label to residue."""
    try:
        c_alpha = residue['CA'].coord
    except KeyError as error:
        raise MissingInformationError from error
    try:
        c_beta = residue['CB'].coord
    except KeyError as error:
        if residue.resname == 'GLY':
            return ChiralCenters.L
        else:
            raise MissingInformationError from error
    n = residue['N'].coord
    c = residue['C'].coord

    assignment = np.dot(c_alpha - c_beta, np.cross(n-c_beta, c-c_beta))
    if assignment > 0:
        return ChiralCenters.D
    elif assignment < 0:
        return ChiralCenters.L
    raise ValueError


def get_chirality(
        struc: Structure.Structure
        ) -> dict[ChiralCenters, list[SingleIrregularity]]:
    """Iterate over structure to yield a defaultdict with `L` and `D` as keys and lists of Residues
as corresponding values."""

    chirality: dict[ChiralCenters, list[SingleIrregularity]] = {e: [] for e in ChiralCenters}
    for residue in struc.get_residues():
        label = assign_chirality_amino_acid(residue)
        chirality[label].append(
                SingleIrregularity(residue, label.value)
                )
    return chirality


def _assign_chirality_amino_acid(residue: Residue.Residue) -> ChiralCenters:
    """deprecated version.
    Assign L/D label to residue."""
    warnings.warn("Deprecated.", DeprecationWarning)
    try:
        c_alpha = residue['CA'].get_vector().copy()
    except KeyError as error:
        raise MissingInformationError from error
    atoms = {atom.name: (atom.get_vector() - c_alpha) for atom in residue}

    try:
        transformer = vectors.rotmat(atoms['HA'], vectors.Vector(0, 0, -1))
    except KeyError as error:
        if residue.resname == 'GLY':
            raise GlycineException from error
        raise MissingInformationError(
                "No H-alphas are in the structure which is needed to determine the chiralities"
                ) from error

    theta_zero = _get_theta(atoms['N'], transformer)
    rotations = {name: np.mod(_get_theta(vec, transformer) - theta_zero, 2*np.pi)
                 for name, vec in atoms.items()}

    chirality = ChiralCenters.L if rotations['C'] < rotations['CB'] else ChiralCenters.D
    return chirality


def _get_theta(vec: vectors.Vector, rotation_matrix: npt.NDArray[np.float64]) -> float:
    """Deprecated. No longer needed.
Rotate vector so that HA lies at x=0, y=0 and then calculate angle of newly obtained vector."""
    warnings.warn("Deprecated.", DeprecationWarning)
    x, y, _ = vec.left_multiply(rotation_matrix)
    theta: float = np.arctan2(y, x)
    return theta
