from __future__ import annotations
from collections import defaultdict
import numpy as np
from Bio.PDB.Residue import Residue
from Bio.PDB import vectors

class GlycineException(Exception):
    """Glycine has no chiral centre"""

def get_theta(vec, rotation_matrix) -> float:
    x, y, _ = vec.left_multiply(rotation_matrix)
    theta = np.arctan2(y, x)
    return theta


def assign_chirality_amino_acid(residue: Residue) -> str:
    CA = residue['CA'].get_vector().copy()
    atoms = {atom.name: atom.get_vector() - CA  for atom in residue}
    
    try:    
        transformer = vectors.rotmat(atoms['HA'], vectors.Vector(0,0,-1))
    except KeyError:
        if residue.resname == 'GLY':
            raise GlycineException
        else:
            raise KeyError

    theta_zero = get_theta(atoms['N'], transformer)
    rotations = {name: np.mod(get_theta(vec, transformer) - theta_zero, 2*np.pi) \
                    for name, vec in atoms.items()}
    
    chirality = 'L' if rotations['C'] < rotations['CB'] else 'D'
    return chirality

def get_chirality(pdb: Structure) -> defaultdict:
    chirality = defaultdict(list)
    for residue in pdb.get_residues():
        res_number = residue.get_id()[1]
        try:
            label = assign_chirality_amino_acid(residue)
        except GlycineException:
            label = "L"
        chirality[label].append(res_number)
    return chirality
