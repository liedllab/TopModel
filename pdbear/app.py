#! /usr/bin/env python3
"""Entry point for CLI"""
from __future__ import annotations

from collections import defaultdict
import itertools
from enum import Enum
from pathlib import Path
import os
import pickle
import subprocess
import sys
import textwrap
import warnings
from typing import NamedTuple, Callable

from Bio.SeqUtils import seq1
from Bio.PDB import PDBParser, Residue, Structure, PDBExceptions
import click

from pdbear import src 

from pdbear.src.utils import PDBError, GlycineException, ProlineException
from pdbear.src.utils import ChiralCenters, StereoIsomers, Clashes, Color



def load_structure(path: str) -> Structure.Structure:
    """Load PDB structure from path"""
    parser = PDBParser()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=PDBExceptions.PDBConstructionWarning)
        structure = parser.get_structure("", path)

    return structure



# handle colors, display and formating

class App:
    """App container"""
    def __init__(self, amides: bool, chiralities: bool, clashes: bool):
        self.funcs: list[Callable] = []
        self.display: dict[Enum, Color] = {}
        if amides:
            self.funcs.append(src.amide_bond.get_stereo)
            self.display.update({
                StereoIsomers.CIS: Color.RED,
                StereoIsomers.CIS_PROLINE: Color.RED,
                StereoIsomers.INBETWEEN: Color.YELLOW,
                })
        if chiralities:
            self.funcs.append(src.chirality.get_chirality)
            self.display.update(
                        {ChiralCenters.D:Color.MAGENTA},
                        )
        if clashes:
            self.funcs.append(src.clashes.get_clashes)
            self.display.update(
                        {Clashes.VDW: Color.BLUE},
                        )
        self.width = os.get_terminal_size().columns
        self.n_res = None
        self.data = None

    def process_structure(self, structure: Structure.Structure) -> None:
        """Calculate chosen measures from structure"""
        structure_data = {}
        for func in self.funcs:
            structure_data.update(func(structure))

        self.n_res = len(str(len(list(structure.get_residues()))))
        self.data = structure_data

    def output_to_terminal(self) -> dict[Enum, list]:
        """Organise the output to terminal"""
        
        res_frmt = lambda residue: f'{seq1(residue.get_resname())}{residue.get_id()[1]:0{self.n_res}}'
        output = defaultdict(list)
        for key, val in self.data.items():
            if key in self.display:
                for item in val:
                    match item:
                        case (_, _):
                            output[key].append('-'.join((res_frmt(res) for res in item)))
                        case Residue.Residue():
                            output[key].append(res_frmt(item))
                        case _:
                            msg = f"For {item} in {key}"
                            raise ValueError(msg)
        
        click.echo("")
        if not output:
            click.echo(click.style("No irregularities were found.", bold=True, fg='green'))
            click.echo("")

        for key, val in output.items():
            self.display_information(key, val)

    def display_information(self, key: Enum, names: list[str]) -> None:
        """Print single information and returns number of residues"""
        
        click.echo(
                click.style(f"{key.name} {key.__class__.__name__}", 
                            bold=True, 
                            fg=self.display[key].value)
                + " have been detected at:")
        raw = ", ".join(names)
        wrapped_text = textwrap.wrap(raw, width=self.width)
        for line in wrapped_text:
            click.echo(line)
        click.echo("")

    def to_pml(self, structure_path: str):
        res_frmt = lambda res: f'resid {res.get_id()[1]}'
        
        commands = [
                f'load {structure_path}', 
                'set cartoon_transparency, 0.3',
                'hide lines',
                'show cartoon',
                'color white',
                ]
        for key, val in self.data.items():
            residues = []
            for item in val:
                match item:
                    case (a, b):
                        residues.append(f'{res_frmt(a)} or {res_frmt(b)}')
                    case c:
                        residues.append(res_frmt(c))
            sel_string = ' or '.join(residues)
            try:
                color = self.display[key].value
            except KeyError:
                pass
            else:
                commands.append(f'color {color}, {sel_string}')
                commands.append(f'show sticks, {sel_string}')
            commands.append(f'select {key.name}, {sel_string}')

        commands.append('color atomic, not elem C')
        commands.append('deselect')

        pml = '\n'.join(commands)
        with open(Path(structure_path).stem + '.pml', 'w') as f:
            f.write(pml)


    def __repr__(self):
        return (f"{self.__class__.name}(width={self.width}px, amides={self.amides},"
                "chiralities={self.chiralities}, clashes={self.clashes})")


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
@click.option("--clashes", is_flag=True, default=True, show_default=True)
@click.option("--pymol", "-p", is_flag=True, default=False, show_default=True)
def main(file: str, amides: bool, chiralities: bool, clashes: bool, pymol: bool) -> None:
    """Check PDB Structures for errors"""
    app = App(amides, chiralities, clashes)
    filepath = Path(file)
    struc = load_structure(filepath)
    try:
        app.process_structure(struc)
    except PDBError as error:
        click.echo(click.style(error, fg='white', bg='red', bold=True))
        sys.exit(1)
    
    app.output_to_terminal()
    app.to_pml(filepath)
    if pymol:
        subprocess.run(['pymol', '-qm', filepath.stem + '.pml'], check=False)


if __name__ == '__main__':
    main() # pylint: disable=no-value-for-parameter
