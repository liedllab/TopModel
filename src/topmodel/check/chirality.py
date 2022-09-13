"""Calculate chirality of an aminoacid."""

from __future__ import annotations

from Bio.PDB import vectors, Residue, Structure
import numpy as np
import numpy.typing as npt
from topmodel.util.errors import GlycineException, MissingInformationError
from topmodel.util.utils import ChiralCenters, SingleIrregularity


def get_chirality(
        struc: Structure.Structure
        ) -> dict[ChiralCenters, list[SingleIrregularity]]:
    """Iterate over structure to yield a defaultdict with `L` and `D` as keys and lists of Residues
as corresponding values."""

    chirality: dict[ChiralCenters, list[SingleIrregularity]] = {e: [] for e in ChiralCenters}
    for residue in struc.get_residues():
        try:
            label = assign_chirality_amino_acid(residue)
        except GlycineException:
            label = ChiralCenters.L

        chirality[label].append(
                SingleIrregularity(residue, label.value)
                )
    return chirality


def assign_chirality_amino_acid(residue: Residue.Residue) -> ChiralCenters:
    """Assign L/D label to residue."""
    c_alpha = residue['CA'].get_vector().copy()
    atoms = {atom.name: (atom.get_vector() - c_alpha) for atom in residue}

    try:
        transformer = vectors.rotmat(atoms['HA'], vectors.Vector(0, 0, -1))
    except KeyError as error:
        if residue.resname == 'GLY':
            raise GlycineException from error
        raise MissingInformationError(
                "No H-alphas are in the structure which is needed to determine the chiralities"
                ) from error

    theta_zero = get_theta(atoms['N'], transformer)
    rotations = {name: np.mod(get_theta(vec, transformer) - theta_zero, 2*np.pi)
                 for name, vec in atoms.items()}

    chirality = ChiralCenters.L if rotations['C'] < rotations['CB'] else ChiralCenters.D
    return chirality


def get_theta(vec: vectors.Vector, rotation_matrix: npt.NDArray[np.float64]) -> float:
    """Rotate vector so that HA lies at x=0, y=0 and then calculate angle of newly obtained
vector."""
    x, y, _ = vec.left_multiply(rotation_matrix)
    theta: float = np.arctan2(y, x)
    return theta
