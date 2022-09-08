#! /usr/bin/env python3
"""Entry point for CLI"""
from __future__ import annotations

from enum import Enum
import itertools
from pathlib import Path
import os
import subprocess
import tempfile
import textwrap
import warnings
from typing import Callable

from Bio.PDB import PDBParser, Structure, PDBExceptions
import click

from pdbear import src
from pdbear.src.utils import ChiralCenters, AmideBonds, Clashes
from pdbear.src.errors import PDBError
from pdbear.src.parser import get_structure


class Color(Enum):
    RED = 'red'
    YELLOW = 'yellow'
    MAGENTA = 'magenta'
    CYAN = 'cyan'


class App:
    """App container"""
    def __init__(self, amides: bool, chiralities: bool, clashes: bool):
        self.funcs: list[Callable] = []
        self.display: dict[Enum, Color] = {}
        if amides:
            self.funcs.append(src.get_amid_stereo)
            self.display.update({
                AmideBonds.CIS: Color.RED,
                AmideBonds.CIS_PROLINE: Color.YELLOW,
                AmideBonds.NON_PLANAR: Color.YELLOW,
                })
        if chiralities:
            self.funcs.append(src.get_chirality)
            self.display.update(
                        {ChiralCenters.D:Color.MAGENTA},
                        )
        if clashes:
            self.funcs.append(src.get_clashes)
            self.display.update(
                        {Clashes.VDW: Color.CYAN},
                        )
        try:
            self.width = os.get_terminal_size().columns
        except OSError:
            self.width = 80

        self.n_res = None
        self.data = None
        self.score = None

    def process_structure(self, structure: Structure.Structure) -> None:
        """Calculate chosen measures from structure"""
        structure_data = {}
        for func in self.funcs:
            structure_data.update(func(structure))

        self.n_res = len(list(structure.get_residues()))
        self.data = structure_data

    def output_to_terminal(self) -> dict[Enum, list]:
        """Organise the output to terminal"""
        if not any(self.data.values()):
            click.echo(click.style("No irregularities were found.", bold=True, fg='green'))

        for key, values in self.data.items():
            if not values:
                continue
            key_string = f'\n{key.name} {key.__class__.__name__}'
            value_string = ', '.join((x.to_cli() for x in values))
            wrapped = textwrap.wrap(value_string, width=self.width)
            try:
                click.echo((click.style(key_string, bold=True, fg=self.display[key].value)
                            + " detected at:"))
            except KeyError:
                continue
            else:
                for line in wrapped:
                    click.echo(line)

    def to_pml(self, structure_path: str) -> str:
        """Transforms output of processed structure into a string which can be used by pymol to open
the structure."""
        commands = [
                f'load {structure_path}',
                'set cartoon_transparency, 0.3',
                'hide lines',
                'show cartoon',
                'color white',
                ]
        for key, values in self.data.items():
            sel_string = ' or '.join((entry.to_pymol() for entry in values))
            try:
                color = self.display[key].value
            except KeyError:
                continue
            else:
                commands.append(f'color {color}, {sel_string}')
                commands.append(f'show sticks, {sel_string}')
            commands.append(f'select {key.name}, {sel_string}')

        commands.append('color atomic, not elem C')
        commands.append('deselect')

        pml = '\n'.join(commands)
        return pml

    def compute_score(self) -> int:
        raw_score = sum(
            (entry.score for entry in itertools.chain(*self.data.values()))
            )
        # wenn bei jeder 3 aminosäure ein fehler ist dann ca 0
        self.score = 1-(3*raw_score / (self.n_res))
        return self.score

    def __repr__(self):
        return f"{self.__class__.__name__}(width={self.width}px)"


@click.command()
@click.argument(
        "files",
        required=True,
        #type=click.Path(exists=True),
        nargs=-1, # unlimited number of arguments
        )
@click.option("--amides/--no-amides", is_flag=True, default=True)
@click.option("--chiralities/--no-chiralities", is_flag=True, default=True)
@click.option("--clashes/--no-clashes", is_flag=True, default=True)
@click.option("--score/--no-score", is_flag=True, default=True)
@click.option("--quiet/--verbose", "-q/-v", is_flag=True, default=False)
@click.option("--pymol", "-p", is_flag=True, default=False, show_default=True,
              help=('Open the structure in PyMOL with the irregularities annotated.'
                  'Requires pymol to be in path')
        )
def main(files: list[str], amides: bool, chiralities: bool, clashes: bool, pymol: bool,
        score: bool, quiet: bool) -> None:
    """Check PDB Structures for errors"""

    app = App(amides, chiralities, clashes)
    paths = (Path(file) for file in files)
    for path in paths:
        click.echo(click.style(f'{path.name.upper()}', bold=True))
        struc = get_structure(str(path))
        try:
            app.process_structure(struc)
        except PDBError as error:
            click.echo(click.style(error, fg='white', bg='red', bold=True))
            raise SystemExit() from error

        if not quiet:
#            click.echo('-'*app.width)
            app.output_to_terminal()
        if score:
            click.echo(click.style('\nSCORE:\t', bold=True) + f'{app.compute_score():#.2%}\n')
        if pymol:
            pml = app.to_pml(path)
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pml') as tmp:
                tmp.write(pml)
                tmp.flush()
                subprocess.call(('pymol', '-Q', tmp.name))


if __name__ == '__main__':
    main() # pylint: disable=no-value-for-parameter
