"""Compute VDW clashes of a Structure"""
from __future__ import annotations

import itertools
from scipy import spatial
from mendeleev import fetch
from Bio.PDB import Structure, Residue, Entity

from topmodel.util.utils import Clashes, CoupleIrregularity


# perfoms this on import
_df = fetch.fetch_table('elements')[['symbol', 'vdw_radius']]
_df = _df.set_index('symbol')
VDW_RADII = _df / 100  # conversion from pm to angstrom
del _df  # so it cant be imported anymore


def get_clashes(struc: Structure.Structure) -> dict[Clashes, list[CoupleIrregularity]]:
    """Iterates over structure and yields a set of clashes."""
    all_clashes = set()
    _info = ((atom.coord, atom.parent) for atom in struc.get_atoms())
    coord, labels = zip(*_info)
    all_points = spatial.KDTree(coord)

    for residue in struc.get_residues():
        clashes = compute_clash(residue, all_points, labels)
        for clash in clashes:
            all_clashes.add(frozenset([residue, clash]))

    return {Clashes.VDW: [CoupleIrregularity(a, b, Clashes.VDW.value)
                          for a, b in all_clashes]}


def compute_clash(residue: Residue.Residue,
                  tree: spatial.KDTree,
                  labels: list[Residue.Residue],
                  ) -> set[Residue.Residue]:
    """Return clashes of a residue as a set using a KDTree and corresponding labels."""
    info = ((atom.coord, VDW_RADII.loc[atom.element.capitalize()].values[0])
            for atom in sidechains(residue))
    coords, radii = zip(*info)
    nearby_points = tree.query_ball_point(coords, radii)
    nearby_labelled = {labels[index] for index in itertools.chain(*nearby_points)}

    return nearby_labelled.difference([residue])


def sidechains(entity: Entity.Entity):
    """Generator that yields sidechain atoms of Entity."""
    for atom in entity.get_atoms():
        if atom.name not in {'C', 'CA', 'N', 'HA'}:
            yield atom
