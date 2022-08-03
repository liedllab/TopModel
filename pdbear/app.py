"""Entry point for CLI"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import os
import pickle
import subprocess
import sys
import textwrap
import warnings

from Bio.SeqUtils import seq1
from Bio.PDB import PDBParser, Residue, Structure, PDBExceptions
import click

from pdbear.src import amide_bond, chirality
from pdbear.src.pdb_errors import PDBError, GlycineException, ProlineException

class App:
    """App container"""
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

    def load_structure(self, path: str) -> Structure.Structure:
        """Load PDB structure from path"""
        parser = PDBParser()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=PDBExceptions.PDBConstructionWarning)
            structure = parser.get_structure("", path)

        return structure

    def process_structure(
            self,
            structure: Structure.Structure
            ) -> dict[str, list[Residue.Residue]]:
        """Calculate chosen measures from structure"""

        header, tailer = structure.get_residues(), structure.get_residues()
        next(tailer)

        structure_data = defaultdict(set)
        structure_data['len'] = 0
        for structure_len, (head, tail) in enumerate(zip(header, tailer), start=1):
            if self.amides:
                try:
                    stereo = amide_bond.assign_stereo(head, tail)
                except ProlineException:
                    stereo = 'cis_proline'
                structure_data[stereo].add(head)
            if self.chiralities:
                if structure_len == 1:
                    try:
                        chiral = chirality.assign_chirality_amino_acid(head)
                    except GlycineException:
                        chiral = 'L'
                    structure_data[chiral].add(head)

                try:
                    chiral = chirality.assign_chirality_amino_acid(tail)
                except GlycineException:
                    chiral = 'L'
                structure_data[chiral].add(tail)
            structure_data['len'] = structure_len

        return structure_data

    def output_to_terminal(self, structure_data: dict[str, list[Residue.Residue]]) -> None:
        """Organise the output to terminal"""
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
                            residue_list: list[Residue.Residue],
                            msg: str,
                            structure_len: int,
                            color: str = 'white') -> None:
        """Print single information"""

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

def run_pymol(script: str | Path,
              structure: str | Path,
              data: dict[str, list[int]],
              app: App,
              temp: str | Path) -> None:
    """Run PyMOL script and handle terminal output"""

    click.echo("-" * app.width)
    with open(temp, 'wb') as tempfile:
        pickle.dump(structure, tempfile, protocol=2)
        pickle.dump(app.colors, tempfile, protocol=2)
        pickle.dump(data, tempfile, protocol=2)

    try:
        subprocess.run(['pymol', '-qm', script], check=False)
    except FileNotFoundError as error:
        click.echo(click.style(error + '\nIs `pymol` in PATH?', fg='white', bg='red', bold=True))
    finally:
        os.remove(temp)
    click.echo("-" * app.width)


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
    """Check PDB Structures for errors"""
    app = App(amides, chiralities)
    while True:
        struc = app.load_structure(file)
        try:
            data = app.process_structure(struc)
        except PDBError as error:
            click.echo(click.style(error, fg='white', bg='red', bold=True))
            sys.exit(1)
        app.output_to_terminal(data)

        if pymol or click.confirm('Do you want to open the structure in PyMOL?', default=False):
            run_pymol(
                    script=Path(__file__).parent / 'script_pymol.py',
                    structure=file,
                    data={k: [v.get_id()[1] for v in values] for k, values in data.items()
                            if k in {'D', 'cis', 'cis_proline', 'strange'}},
                    app=app,
                    temp='.pdbear.temp',
                    )

        click.confirm('Do you want to continue?', abort=True, default=True)
        file = click.prompt("Next Structure", type=click.Path(exists=True))


if __name__ == '__main__':
    main() # pylint: disable=no-value-for-parameter
