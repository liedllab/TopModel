from collections import defaultdict
from Bio.PDB.Structure import Structure
from Bio.PDB.vectors import calc_dihedral
import numpy as np


class ProlineException(Exception):
    """Prolines are more likely to be in cis-conformation, especially if they are preceded by
    Glycine or an aromatic residue."""


def get_stereo(pdb: Structure) -> dict[int, tuple[str, str, float]]:
    stereo = defaultdict(list)
    
    header = pdb.get_residues()
    tailer = pdb.get_residues()
    next(tailer)
    for head, tail in zip(header, tailer):
        try:
            label = assign_stereo(head, tail)
        except ProlineException:
            label = '?'
        res_number = head.get_id()[1]
        stereo[label].append(res_number) 
    return stereo


def assign_stereo(head, tail):
    angle = calc_dihedral(
            head['CA'].get_vector(), head['C'].get_vector(),
            tail['N'].get_vector(), tail['CA'].get_vector(),
            )
    angle = np.mod(angle, 2*np.pi)
    if 5*np.pi/6 <= angle <= 7*np.pi/6:
        label = 'trans'
    elif (0 <= angle <= np.pi/6) or (11*np.pi/6 <= angle <= 2*np.pi):
        if tail.resname == 'PRO':
            raise ProlineException
        label = 'cis'
    else:
        label = '?'
    return label

