from __future__ import annotations

from collections import defaultdict
import textwrap
import os
import sys
import warnings

from Bio.SeqUtils import seq1
from Bio.PDB.Structure import Structure
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import click

from pdbear import amide_bond, chirality

class PDBError(Exception):
    """Raised when needed information is missing in Structure"""

class App:
    def __init__(self, amides, chiralities):
        self.width = os.get_terminal_size().columns
        self.amides = amides
        self.chiralities = chiralities

    def load_structure(self, path: str) -> Structure:
        parser = PDBParser()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=PDBConstructionWarning)
            structure = parser.get_structure("", path)

        all_atom_types = {atom.name for atom in structure.get_atoms()}

        if self.chiralities and "HA" not in all_atom_types:
            raise PDBError(
                    "No H-alphas are in the PDB structure which is needed to \
determine the chiralities"
                    )
        if self.amides:
            diff = all_atom_types.difference({'CB', 'C', 'CA', 'N'})
            if diff:
                diff_str = ", ".join(diff)
                raise PDBError(f"{diff_str} needed to determine the amide bond orientation")

        return structure

    def process_structure(self, structure) -> None:
        header, tailer = structure.get_residues(), structure.get_residues()
        next(tailer) 
        
        structure_data = defaultdict(set)
        for structure_len, (head, tail) in enumerate(zip(header, tailer), start=1):
            if self.amides:
                try:
                    stereo = amide_bond.assign_stereo(head, tail)
                except amide_bond.ProlineException:
                    stereo = 'cis_proline'
                structure_data[stereo].add(head)
            if self.chiralities:
                try:
                    chiral = chirality.assign_chirality_amino_acid(head)
                except chirality.GlycineException:
                    chiral = 'L'
                structure_data[chiral].add(head)

        if self.chiralities:
            try:
                chiral = chirality.assign_chirality_amino_acid(tail)
            except chirality.GlycineException:
                chiral = 'L'
            structure_data[chiral].add(tail)

        structure_data['len'] = structure_len

        return structure_data

    def output_to_terminal(self, structure_data):
        click.echo("")
        if self.chiralities:
            self.display_information(
                    structure_data['D'],
                    "D amino acids",
                    structure_data['len'],
                    )
        # amide
        if self.amides:
            self.display_information(
                    structure_data['cis'],
                    "Cis amide bonds",
                    structure_data['len'],
                    )
            self.display_information(
                    structure_data['cis_proline'],
                    "Bonds to cis prolines",
                    structure_data['len'],
                    )
            self.display_information(
                    structure_data['?'],
                    "Strange torsion angles for the amide bond",
                    structure_data['len'],
                    )

        if not structure_data['D'] and not structure_data['cis'] and not structure_data['?']:
            click.echo(click.style("No irregularities were found.", bold=True, fg='green'))
            click.echo("")

        return None
    
    def display_information(
            self, 
            residue_list: list[Residue], msg: str, 
            structure_len: int) -> None:
        
        if not residue_list:
            return None

        click.echo(click.style(msg, bold=True, fg='red') + " have been detected at:")
        raw = ", ".join(
                [f"{seq1(residue.resname)}{residue.get_id()[1]:0{len(str(structure_len))}}" \
                        for residue in residue_list]
                )
        wrapped_text = textwrap.wrap(raw, width=self.width)
        for line in wrapped_text:
            click.echo(line)
        click.echo("")



@click.command()
@click.option(
        "--file", "-f",
        required=True,
        type=click.Path(exists=True),
        help="Path to the PDB file",
        )
@click.option("--amides", "-a", is_flag=True, default=True, show_default=True)
@click.option("--chiralities", "-c", is_flag=True, default=True, show_default=True)
@click.option("--pymol", "-p", is_flag=True, default=False, show_default=True)
def main(file, amides, chiralities, pymol) -> None:
    app = App(amides, chiralities)
    try:
        struc = app.load_structure(file)
    except PDBError as e:
        click.echo(click.style(e, fg='white', bg='red', bold=True))
        sys.exit(1)
    data = app.process_structure(struc)
    app.output_to_terminal(data)


if __name__ == '__main__':
    sys.tracebacklimit=0
    main()
