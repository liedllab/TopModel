from __future__ import annotations
import numpy as np
from Bio.PDB.Residue import Residue


def set_zero_point(values: np.ndarray, zero_point: np.ndarray) -> np.ndarray:
    return values - zero_point

def get_projector(vector: np.ndarray) -> np.ndarray:
    A = np.array(
            [vector,
             [0, 1, 0],
             [0, 0, 1]]
            )


    q, r = np.linalg.qr(A.T)
    if np.sign(q.T @ vector)[0] == 1:
        q = -1 * q
    return q.T

def assign_chirality(residue: Residue) -> str:

    atoms = {atom.name : atom.coord.copy() for atom in residue.get_atoms() \
             if atom.name in ('C', 'CA', 'CB', 'N', 'HA')}
    zero = atoms.pop('CA')
    for atom, position in atoms.items():
        atoms[atom] = set_zero_point(position, zero)

    transformer = get_projector(atoms.pop('HA'))[1:]
    
    get_theta = lambda pos: np.mod(np.arctan2(*reversed(transformer @ pos)), 2*np.pi)
    theta_zero = get_theta(atoms['N'])
    rotations = {atom : get_theta(position) - theta_zero for atom, position in atoms.items()}

    label = 'L' if rotations['C'] < rotations['CB'] else 'D'

    return label


def get_chiralities(pdb: Structure) -> dict[int, tuple[str, str]]:
    result = dict()
    for residue in pdb.get_residues():
        atoms = {atom.name : atom.coord for atom in residue.get_atoms()}
        try: 
            h_vec = atoms['HA'] - atoms['CA']
        except KeyError:
            if residue.get_resname() == "GLY":
                label = "None"
            else:
                raise KeyError
        else:
            rotator = get_rotator(h_vec)
            rotations = dict()
            for atom in ['N', 'C', 'CB']:
                rotated_xyz = rotator @ (atoms[atom] - atoms['CA'])
                rotation = np.arctan2(rotated_xyz[1], rotated_xyz[0])
                rotations[atom] = rotation
                rotations[atom] -= rotations['N']
                rotations[atom] = np.mod(rotations[atom], 2*np.pi)


            label = "D" if rotations["C"] < rotations['CB'] else "L"
        result[residue.get_id()[1]] = residue.get_resname(), label
            
    return result


def get_theta(vec: np.array) -> float:
    return np.arctan2(vec[1], vec[2])


def get_phi(vec: np.array) -> float:
    return np.arctan2(-vec[0], vec[2])


def rotate_x(theta: float) -> np.array:
    sin = np.sin(theta)
    cos = np.cos(theta)

    M = np.array(
            [[1,   0,    0],
             [0, cos, -sin],
             [0, sin,  cos]]
            )
    return M


def rotate_y(phi: float) -> np.array:
    sin = np.sin(phi)
    cos = np.cos(phi)

    M = np.array(
            [[ cos,  0, sin],
             [   0,  1,   0],
             [-sin,  0, cos]]
            )
    return M

def get_rotator(vec: np.array) -> np.array:
    theta = get_theta(vec)
    intermediate = rotate_x(theta) @ vec
    phi = get_phi(intermediate)
    rotated = rotate_y(phi) @ intermediate

    factor = -1 if rotated[2] >= 0 else 1

    transformed = factor * rotate_y(phi) @ rotate_x(theta)
    return transformed
