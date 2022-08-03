"""Calculate chirality of an aminoacid."""

from __future__ import annotations

from collections import defaultdict
from Bio.PDB import vectors, Residue, Structure
import numpy as np
from .utils import GlycineException, PDBError, ChiralCenter


def get_chirality(pdb: Structure.Structure) -> dict[str, list[Residue.Residue]]:
    """Iterate over structure to yield a defaultdict with `L` and `D` as keys and lists of Residues
as corresponding values."""

    chirality = defaultdict(list)
    for residue in pdb.get_residues():
        res_number = residue.get_id()[1]
        try:
            label = assign_chirality_amino_acid(residue)
        except GlycineException:
            label = ChiralCenter.L
        chirality[label].append(res_number)
    return chirality


def assign_chirality_amino_acid(residue: Residue.Residue) -> str:
    """Assign L/D label to residue."""

    CA = residue['CA'].get_vector().copy()
    atoms = {atom.name: atom.get_vector() - CA  for atom in residue}

    try:
        transformer = vectors.rotmat(atoms['HA'], vectors.Vector(0,0,-1))
    except KeyError as error:
        if residue.resname == 'GLY':
            raise GlycineException from error
        raise PDBError(
        "No H-alphas are in the PDB structure which is needed to \
determine the chiralities") from error

    theta_zero = get_theta(atoms['N'], transformer)
    rotations = {name: np.mod(get_theta(vec, transformer) - theta_zero, 2*np.pi) \
                    for name, vec in atoms.items()}

    chirality = ChiralCenter.L if rotations['C'] < rotations['CB'] else ChiralCenter.D
    return chirality


def get_theta(vec: vectors.Vector, rotation_matrix: np.ndarray) -> float:
    """Rotate vector so that HA lies at x=0, y=0 and then calculate angle of newly obtained
vector."""

    x, y, _ = vec.left_multiply(rotation_matrix)
    theta = np.arctan2(y, x)
    return theta
