"""Calculate and assign `cis` and `trans` label to amide bonds. Bonds that are neither cis nor trans
are assigned the `strange` label."""

from __future__ import annotations

from collections import defaultdict
from Bio.PDB import Structure, Residue, vectors
import numpy as np
from .errors import ProlineException, MissingInformationError
from .utils import AmideBonds, StructuralIrregularity, CoupleIrregularity


def get_amid_stereo(struc: Structure.Structure) -> dict[str, list[StructuralIrregularity]]:
    """Iterates over structure and yields result in a dictionary that maps the label to a list of
Residues."""

    stereo = defaultdict(list)
    header = struc.get_residues()
    tailer = struc.get_residues()
    next(tailer) # advance second iterator so tailer is always one step ahead of header
    for head, tail in zip(header, tailer):
        try:
            label = assign_stereo(head, tail)
        except ProlineException:
            label = AmideBonds.CIS_PROLINE

        stereo[label].append(CoupleIrregularity(head, tail, label.value))
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
        raise MissingInformationError(f"{diff_str} needed to determine the amide bond orientation")\
                from error

    if head['C'] - tail['N'] > 2:
        return AmideBonds.TRANS
    angle = np.mod(angle, 2*np.pi)
    if 5*np.pi/6 <= angle <= 7*np.pi/6:
        return AmideBonds.TRANS
    if (0 <= angle <= np.pi/6) or (11*np.pi/6 <= angle <= 2*np.pi):
        if tail.resname == 'PRO':
            raise ProlineException
        return AmideBonds.CIS
    return AmideBonds.NON_PLANAR
