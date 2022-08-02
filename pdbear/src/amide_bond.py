"""Calculate and assign `cis` and `trans` label to amide bonds. Bonds that are neither cis nor trans
are assigned the `strange` label."""

from __future__ import annotations

from collections import defaultdict
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.vectors import calc_dihedral
import numpy as np
from .pdb_errors import ProlineException


def get_stereo(pdb: Structure) -> dict[str, list[Residue]]:
    """Iterates over structure and yields result in a dictionary that maps the label to a list of
Residues."""

    stereo = defaultdict(list)
    header = pdb.get_residues()
    tailer = pdb.get_residues()
    next(tailer)
    for head, tail in zip(header, tailer):
        try:
            label = assign_stereo(head, tail)
        except ProlineException:
            label = 'strange'
        res_number = head.get_id()[1]
        stereo[label].append(res_number)
    return stereo


def assign_stereo(head: Residue, tail: Residue) -> str:
    """Calculates dihedral angle of amide bond between Residue `head` and `tail`. The angle then is
mapped to the labels `cis`, `trans` or `strange`."""

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
        label = 'strange'
    return label
