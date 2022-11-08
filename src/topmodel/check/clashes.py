"""Compute VDW clashes of a Structure"""
from __future__ import annotations
from collections import defaultdict
from scipy.spatial import KDTree
from mendeleev import fetch
from Bio.PDB import Structure, Entity

from topmodel.util.utils import Clashes, CoupleIrregularity


# perfoms this on import
_df = fetch.fetch_table('elements')[['symbol', 'vdw_radius']]
_df = _df.set_index('symbol')
_df = _df / 100  # conversion from pm to angstrom
VDW_RADII = _df.to_dict()['vdw_radius']
del _df  # so it cant be imported anymore


def get_clashes(struc: Structure.Structure) -> dict[Clashes, list[CoupleIrregularity]]:
    """Iterates over structure and yields a set of clashes."""

    atoms = []              # contains references to the atom object
    residues = []           # contains references to the residue object
    atomic_coordinates = []
    atom_to_residue = {}    # will have atom index as key and residue number as value
    clashes: dict[int, set[int]] = defaultdict(set)
    backbone = {'C', 'CA', 'N', 'O', 'CD'}  # CD in prolines

    atoms_so_far = 0
    for residue_index, residue in enumerate(struc.get_residues()):
        residues.append(residue)
        for atom_index, atom in enumerate(heavy(residue), start=atoms_so_far):
            atoms.append(atom)                     # reference to atom object
            atomic_coordinates.append(atom.coord)  # coordinates of atom object

            atom_to_residue[atom_index] = residue_index  # mapping of atom index to res index
        atoms_so_far += sum(1 for _ in heavy(residue))   # keep count for proper atom indexing

    tree = KDTree(atomic_coordinates)
    close_points = tree.query_pairs(5)  # within 5 angstroms

    for a, b in close_points:
        res_a, res_b = atom_to_residue[a], atom_to_residue[b]
        if res_a != res_b and \
           (res_a not in clashes or res_b not in clashes[res_a]) and \
           (res_b - res_a != 1 or ((atoms[a].name in backbone) != (atoms[b].name in backbone))):
            # distance of pair is only calculated if a and b:
            # not same residue
            # and not already computed
            # if next to each other, make sure not both are in the backbone
            distance = atoms[a] - atoms[b]
            combined_radii = VDW_RADII[atoms[a].element] + VDW_RADII[atoms[b].element]

            if distance < combined_radii - 0.5 and not \
               (residues[res_a].resname == 'CYX' and residues[res_b].resname == 'CYX'):
                # ignore 'clash' if its a 'clash' between cysteine disulfide bridges
                clashes[res_a].add(res_b)

    return {Clashes.VDW: [CoupleIrregularity(residues[a], residues[b], Clashes.VDW.value)
                          for a, values in sorted(clashes.items())
                          for b in sorted(values)]
            }


def heavy(entity: Entity.Entity):
    """Generator that excludes Hydrogens."""
    for atom in entity.get_atoms():
        if atom.element != 'H':
            yield atom
