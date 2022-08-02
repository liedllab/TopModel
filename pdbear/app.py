"""Entry point for CLI"""
from __future__ import annotations

from collections import defaultdict
import textwrap
import os
import sys
import warnings
import pickle
import subprocess
from pathlib import Path

from Bio.SeqUtils import seq1
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import click

from pdbear.src import amide_bond, chirality
from pdbear.src.pdb_errors import PDBError

class App:
    def __init__(self, amides, chiralities):
        self.width = os.get_terminal_size().columns
        self.amides = amides
        self.chiralities = chiralities
        self.colors = {
                'D': 'magenta',
                'cis': 'red',
                'cis_proline': 'red',
                'strange': 'yellow',
                }

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
            diff = {'CB', 'C', 'CA', 'N'}.difference(all_atom_types)
            if diff:
                diff_str = ", ".join(diff)
                raise PDBError(f"{diff_str} needed to determine the amide bond orientation")

        return structure

    def process_structure(self, structure: Structure) -> dict[str, list[Residue]]:
        header, tailer = structure.get_residues(), structure.get_residues()
        next(tailer)

        structure_data = defaultdict(set)
        structure_data['len'] = 0
        for structure_len, (head, tail) in enumerate(zip(header, tailer), start=1):
            if self.amides:
                try:
                    stereo = amide_bond.assign_stereo(head, tail)
                except amide_bond.ProlineException:
                    stereo = 'cis_proline'
                structure_data[stereo].add(head)
            if self.chiralities:
                if structure_len == 1:
                    try:
                        chiral = chirality.assign_chirality_amino_acid(head)
                    except chirality.GlycineException:
                        chiral = 'L'
                    structure_data[chiral].add(head)

                try:
                    chiral = chirality.assign_chirality_amino_acid(tail)
                except chirality.GlycineException:
                    chiral = 'L'
                structure_data[chiral].add(tail)
            structure_data['len'] = structure_len

        return structure_data

    def output_to_terminal(self, structure_data: dict[str, list[Residue]]) -> None:
        click.echo("")
        if self.chiralities:
            self.display_information(
                    structure_data['D'],
                    "D amino acids",
                    structure_data['len'],
                    self.colors['D'],
                    )
        # amide
        if self.amides:
            self.display_information(
                    structure_data['cis'],
                    "Cis amide bonds",
                    structure_data['len'],
                    self.colors['cis'],
                    )
            self.display_information(
                    structure_data['cis_proline'],
                    "Bonds to cis prolines",
                    structure_data['len'],
                    self.colors['cis_proline'],
                    )
            self.display_information(
                    structure_data['strange'],
                    "Strange torsion angles for the amide bond",
                    structure_data['len'],
                    self.colors['strange'],
                    )

        if not structure_data['D'] and not structure_data['cis'] and not structure_data['strange']:
            click.echo(click.style("No irregularities were found.", bold=True, fg='green'))
            click.echo("")

    def display_information(self,
                            residue_list: list[Residue],
                            msg: str,
                            structure_len: int,
                            color: str = 'white') -> None:

        if not residue_list:
            return None

        click.echo(click.style(msg, bold=True, fg=color) + " have been detected at:")
        raw = ", ".join(
                [f"{seq1(residue.resname)}{residue.get_id()[1]:0{len(str(structure_len))}}" \
                        for residue in residue_list]
                )
        wrapped_text = textwrap.wrap(raw, width=self.width)
        for line in wrapped_text:
            click.echo(line)
        click.echo("")
        return None


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
def main(file: str, amides: bool, chiralities: bool, pymol: bool) -> None:
    """CLI application"""
    app = App(amides, chiralities)
    while True:
        try:
            struc = app.load_structure(file)
        except PDBError as e:
            click.echo(click.style(e, fg='white', bg='red', bold=True))
            sys.exit(1)
        data = app.process_structure(struc)
        app.output_to_terminal(data)

        if pymol or click.confirm('Do you want to open the structure in PyMOL?', default=False):
            click.echo("-"*app.width)
            clean = {k:[v.get_id()[1] for v in values] for k, values in data.items() \
                    if k in {'D', 'cis', 'cis_proline', 'strange'}}

            with open('.pdbear_temp', 'wb') as f:
                pickle.dump(file, f, protocol=2)
                pickle.dump(app.colors, f, protocol=2)
                pickle.dump(clean, f, protocol=2)

            path = Path(__file__)
            try:
                subprocess.run(['pymol', '-qm', path.parent / 'script_pymol.py'], check=False)
            except FileNotFoundError as e:
                click.echo(click.style(e + '\nIs `pymol` in PATH?', fg='white', bg='red', bold=True))
            click.echo("-"*app.width)

        click.confirm('Do you want to continue?', abort=True, default=True)
        file = click.prompt("Next Structure", type=click.Path(exists=True))


if __name__ == '__main__':
    main() # pylint: disable=no-value-for-parameter
