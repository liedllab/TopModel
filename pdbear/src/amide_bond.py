"""Calculate and assign `cis` and `trans` label to amide bonds. Bonds that are neither cis nor trans
are assigned the `strange` label."""

from __future__ import annotations

from collections import defaultdict
from Bio.PDB import Structure, Residue, vectors
import numpy as np
from .utils import ProlineException, PDBError, StereoIsomer


def get_stereo(pdb: Structure.Structure) -> dict[str, list[Residue.Residue]]:
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
            label = StereoIsomer.CIS_PROLINE
        res_number = head.get_id()[1]
        stereo[label].append(res_number)
    return stereo


def assign_stereo(head: Residue.Residue, tail: Residue.Residue) -> str:
    """Calculates dihedral angle of amide bond between Residue `head` and `tail`. The angle then is
mapped to the labels `cis`, `trans` or `strange`."""
    try:
        angle = vectors.calc_dihedral(
                head['CA'].get_vector(), head['C'].get_vector(),
                tail['N'].get_vector(), tail['CA'].get_vector(),
                )
    except KeyError as error:
        head_atoms = {atom.get_name() for atom in head.get_atoms()}
        tail_atoms = {atom.get_name() for atom in tail.get_atoms()}
        diff = {'C', 'CA', 'N'}.difference(head_atoms.union(tail_atoms))
        diff_str = ", ".join(diff)
        raise PDBError(f"{diff_str} needed to determine the amide bond orientation") from error

    angle = np.mod(angle, 2*np.pi)
    if 5*np.pi/6 <= angle <= 7*np.pi/6:
        label = StereoIsomer.TRANS
    elif (0 <= angle <= np.pi/6) or (11*np.pi/6 <= angle <= 2*np.pi):
        if tail.resname == 'PRO':
            raise ProlineException
        label = StereoIsomer.CIS
    else:
        label = StereoIsomer.INBETWEEN
    return label
