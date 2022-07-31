from __future__ import annotations
import numpy as np
from Bio.PDB.Residue import Residue
from .logger import logger

def set_zero_point(values: np.ndarray, zero_point: np.ndarray) -> np.ndarray:
    return values - zero_point

def get_transform(vector: np.ndarray) -> np.ndarray:
    A = np.array(
            [vector,
             [0, 1, 0],
             [0, 0, 1]]
            )

    q, r = np.linalg.qr(A.T)
    transform = q.T
    if np.sign(transform @ vector)[0] == 1:
        rot = np.array([[-1,0,0],[0,1,0],[0,0,1]])
        transform = rot @ transform
    if np.sign(np.linalg.det(transform)) == -1:
        transform[[1,2]] = transform[[2,1]]
    return transform


def get_theta(position: np.ndarray) -> float:
    theta = np.arctan2(*reversed(position))
    normalized_theta = np.mod(theta, 2*np.pi)
    return normalized_theta


def assign_chirality_amino_acid(residue: Residue) -> str:
    atoms = {atom.name: atom.coord.copy() for atom in residue.get_atoms() \
                if atom.name in ('C', 'CA', 'CB', 'N', 'HA')}
    zero = atoms['CA'].copy()
    for atom, position in atoms.items():
        atoms[atom] = set_zero_point(position, zero)

    transformer = get_transform(atoms['HA'])
    factor = np.sign(np.linalg.det(transformer))
    transformed_pos = {atom: transformer @ position for atom, position in atoms.items()}

    theta_zero = get_theta(transformed_pos['N'][1:])
    rotations = {atom : get_theta(position[1:]) - theta_zero \
                    for atom, position in transformed_pos.items()}
    
    if rotations['C'] < rotations['CB']:
        chirality = 'L'
    else:
        chirality = 'D'
    
    logger.warning(
        f"AminoAcid {residue} has been asigned to be {chirality}\
            {np.linalg.det(transformer)=}")
    return chirality


