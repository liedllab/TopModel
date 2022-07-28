from Bio.PDB.Structure import Structure
import numpy as np

def get_stereo(pdb: Structure) -> dict[int, tuple[str, str, float]]:
    mask = {"N", "CA", "C", "O"}
    atoms = (atom for atom in pdb.get_atoms() if atom.name in mask)

    stereo = dict()
    for atom in atoms:
        if atom.name == 'C':
            try:
                c = atom
                o = next(atoms)
                n = next(atoms)
                ca = next(atoms)
                # between = {p.get_id()[1] for p in (c, o, n, ca)}
            except StopIteration:
                break
            else:
                phi = calc_angle(
                        c.coord - o.coord,
                        n.coord - c.coord,
                        ca.coord - n.coord,
                        )
                if 11*np.pi/6 <= phi or phi <= np.pi/6:
                    label = "c"
                    label = 't'
                else:
                    label = 't'

                parent = c.get_parent()
                stereo[parent.get_id()[1]] = (parent, ca.get_parent(), label, phi)
    return stereo


def calc_angle(u1: np.array, u2: np.array, u3: np.array) -> float:
    cl = np.cross(u1, u2)
    cr = np.cross(u2, u3)
    n_cl = np.linalg.norm(cl)
    n_cr = np.linalg.norm(cr)

    cos_phi = np.dot(cl, cr) / np.dot(n_cl, n_cr)
    phi = np.arccos(cos_phi)
    while phi < 0:
        phi += 2*np.pi

    return phi