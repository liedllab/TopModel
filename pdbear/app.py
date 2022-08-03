"""Entry point for CLI"""
from __future__ import annotations

from collections import defaultdict
from enum import Enum
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
from pdbear.src.utils import PDBError, GlycineException, ProlineException
from pdbear.src.utils import ChiralCenter, StereoIsomer

class App:
    """App container"""
    def __init__(self, amides, chiralities):
        self.width = os.get_terminal_size().columns
        self.amides = amides
        self.chiralities = chiralities
        self.colors = {
                ChiralCenter.D: 'magenta',
                StereoIsomer.CIS: 'red',
                StereoIsomer.CIS_PROLINE: 'red',
                StereoIsomer.INBETWEEN: 'yellow',
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
            ) -> tuple[dict[str, list[Residue.Residue]], int]:
        """Calculate chosen measures from structure"""

        header, tailer = structure.get_residues(), structure.get_residues()
        next(tailer)

        structure_data = defaultdict(set)
        structure_len = 0
        for structure_len, (head, tail) in enumerate(zip(header, tailer), start=1):
            if self.amides:
                try:
                    stereo = amide_bond.assign_stereo(head, tail)
                except ProlineException:
                    stereo = StereoIsomer.CIS_PROLINE
                structure_data[stereo].add(head)
            if self.chiralities:
                if structure_len == 1:
                    try:
                        chiral = chirality.assign_chirality_amino_acid(head)
                    except GlycineException:
                        chiral = ChiralCenter.L
                    structure_data[chiral].add(head)

                try:
                    chiral = chirality.assign_chirality_amino_acid(tail)
                except GlycineException:
                    chiral = ChiralCenter.L
                structure_data[chiral].add(tail)

        return structure_data, structure_len

    def output_to_terminal(self,
                           structure_data: dict[Enum, list[Residue.Residue]],
                           length: int) -> None:
        """Organise the output to terminal"""
        click.echo("")
        counter = 0
        if self.chiralities:
            counter += self.display_information(
                    structure_data[ChiralCenter.D],
                    "D amino acids",
                    length,
                    self.colors[ChiralCenter.D],
                    )
        # amide
        if self.amides:
            counter += self.display_information(
                    structure_data[StereoIsomer.CIS],
                    "Cis amide bonds",
                    length,
                    self.colors[StereoIsomer.CIS],
                    )
            counter += self.display_information(
                    structure_data[StereoIsomer.CIS_PROLINE],
                    "Bonds to cis prolines",
                    length,
                    self.colors[StereoIsomer.CIS_PROLINE],
                    )
            counter += self.display_information(
                    structure_data[StereoIsomer.INBETWEEN],
                    "Strange torsion angles for the amide bond",
                    length,
                    self.colors[StereoIsomer.INBETWEEN],
                    )

        if counter == 0:
            click.echo(click.style("No irregularities were found.", bold=True, fg='green'))
            click.echo("")

    def display_information(self,
                            residue_list: list[Residue.Residue],
                            msg: str,
                            structure_len: int,
                            color: str = 'white') -> int:
        """Print single information and returns number of residues"""

        if not residue_list:
            return 0

        click.echo(click.style(msg, bold=True, fg=color) + " have been detected at:")
        raw = ", ".join(
                [f"{seq1(residue.resname)}{residue.get_id()[1]:0{len(str(structure_len))}}" \
                        for residue in residue_list]
                )
        wrapped_text = textwrap.wrap(raw, width=self.width)
        for line in wrapped_text:
            click.echo(line)
        click.echo("")
        return len(residue_list)

def run_pymol(script: str | Path,
              structure: str | Path,
              data: dict[Enum, list[int]],
              app: App,
              temp: str | Path) -> None:
    """Run PyMOL script and handle terminal output"""

    click.echo("-" * app.width)
    with open(temp, 'wb') as tempfile:
        pickle.dump(structure, tempfile, protocol=2)
        pickle.dump(
                {k.name : values for k, values in app.colors.items()},
                tempfile, protocol=2)
        pickle.dump(
                {k.name : [v.get_id()[1] for v in values] \
                        for k, values in data.items()},
                tempfile, protocol=2)

    try:
        subprocess.run(['pymol', '-qm', script], check=False)
    except FileNotFoundError as error:
        click.echo(click.style(error + '\nIs `pymol` in PATH?', fg='white', bg='red', bold=True))
    finally:
        os.remove(temp)
    click.echo("-" * app.width)


@click.command()
@click.argument(
        "file",
        required=True,
        type=click.Path(exists=True),
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
            data, n_res = app.process_structure(struc)
            print(data)
        except PDBError as error:
            click.echo(click.style(error, fg='white', bg='red', bold=True))
            sys.exit(1)
        app.output_to_terminal(data, n_res)

        if pymol or click.confirm('Do you want to open the structure in PyMOL?', default=False):
            run_pymol(
                    script=Path(__file__).parent / 'script_pymol.py',
                    structure=file,
                    data=data,
                    app=app,
                    temp='.pdbear.temp',
                    )

        click.confirm('Do you want to continue?', abort=True, default=True)
        file = click.prompt("Next Structure", type=click.Path(exists=True))


if __name__ == '__main__':
    main() # pylint: disable=no-value-for-parameter
