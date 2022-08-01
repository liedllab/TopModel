from collections import defaultdict
from Bio.PDB.Structure import Structure
from Bio.PDB.vectors import calc_dihedral
import numpy as np

def get_stereo(pdb: Structure) -> dict[int, tuple[str, str, float]]:
    mask = {"N", "CA", "C", "O"}
    bb_iter = (atom for atom in pdb.get_atoms() if atom.name in mask)
    stereo = defaultdict(list)
    for atom in bb_iter:
        while atom.name != 'C':
            atom = next(bb_iter)

        C = atom.get_vector()
        try:
            O = next(bb_iter).get_vector()
            N = next(bb_iter).get_vector()
            CA = next(bb_iter).get_vector()
            # between = {p.get_id()[1] for p in (c, o, n, ca)}
        except StopIteration:
            break
        
        angle = calc_dihedral(C, O, N, CA)
        label = assign_stereo(np.mod(angle, 2*np.pi))
        res_number = atom.get_parent().get_id()[1]
        stereo[label].append(res_number) 
    return stereo

def assign_stereo(angle: float) -> str:
    if 5*np.pi/6 <= angle <= 7*np.pi/6:
        return "trans"
    elif (0 <= angle <= np.pi/6) or (11*np.pi/6 <= angle <= 2*np.pi):
        return "cis"
    return "?"
