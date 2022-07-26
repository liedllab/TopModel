import numpy as np
from Bio.PDB.Structure import Structure


def get_chiralities(pdb: Structure) -> list[str]:
    result = list()
    for residue in pdb.get_residues():
        atoms = {atom.name : atom.coord for atom in residue.get_atoms()}

        h_vec = atoms['H'] - atoms['CA']
        rotate_side = get_rotator(h_vec) @ ( atoms['CB'] - atoms['CA'])

        if rotate_side[1] < 0:
            result.append("L")
        elif rotate_side[1] > 0:
            result.append("D")
        else:
            raise ValueError("No chiral center found")
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
