from __future__ import annotations
import numpy as np
from Bio.PDB.Residue import Residue


def set_zero_point(values: np.ndarray, zero_point: np.ndarray) -> np.ndarray:
    return values - zero_point

def get_transform(vector: np.ndarray) -> np.ndarray:
    A = np.array(
            [vector,
             [0, 1, 0],
             [0, 0, 1]]
            )


    q, r = np.linalg.qr(A.T)
    if np.sign(q.T @ vector)[0] == 1:
        q = -1 * q
    return q.T


def get_theta(position: np.ndarray) -> float:
    theta = np.arctan2(*reversed(position))
    normalized_theta = np.mod(theta, 2*np.pi)
    return normalized_theta

def assign_chirality_amino_acid(residue: Residue) -> str:

    atoms = {atom.name : atom.coord.copy() for atom in residue.get_atoms() \
                if atom.name in ('C', 'CA', 'CB', 'N', 'HA')}
    zero = atoms.pop('CA')
    for atom, position in atoms.items():
        atoms[atom] = set_zero_point(position, zero)

    transformer = get_transform(atoms.pop('HA'))[1:]
    
    theta_zero = get_theta(transformer @ atoms['N'])
    rotations = {atom : get_theta(transformer @ position) - theta_zero \
                    for atom, position in atoms.items()}

    chirality = 'L' if rotations['C'] < rotations['CB'] else 'D'

    return chirality


