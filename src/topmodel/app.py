#! /usr/bin/env python3
"""Entry point for CLI"""
from __future__ import annotations
from enum import Enum
from pathlib import Path
import os
import subprocess
import tempfile
import textwrap
from typing import Callable

from Bio.PDB.Structure import Structure
import click

from topmodel import check
from topmodel.util.utils import ChiralCenters, AmideBonds
from topmodel.util.utils import Clashes, StructuralIrregularity
from topmodel.util.errors import MissingInformationError, PDBCodeError
from topmodel.util.parser import get_structure


class Color(Enum):
    """Usable colors by colorama/click."""
    RED = 'red'
    YELLOW = 'yellow'
    MAGENTA = 'magenta'
    CYAN = 'cyan'
    GREEN = 'green'


class App:
    """App container"""
    def __init__(self, amides: bool, chiralities: bool, clashes: bool):
        self.funcs: list[Callable[[Structure], dict[Enum, list[StructuralIrregularity]]]] = []
        self.display: dict[Enum, Color] = {}
        if chiralities:
            self.funcs.append(check.get_chirality)  # type: ignore[attr-defined]
            self.display.update(
                        {ChiralCenters.D: Color.MAGENTA},
                        )
        if amides:
            self.funcs.append(check.get_amide_stereo)  # type: ignore[attr-defined]
            self.display.update({
                AmideBonds.CIS: Color.RED,
                AmideBonds.CIS_PROLINE: Color.GREEN,
                AmideBonds.NON_PLANAR: Color.CYAN,
                })
        if clashes:
            self.funcs.append(check.get_clashes)  # type: ignore[attr-defined]
            self.display.update(
                        {Clashes.VDW: Color.YELLOW},
                        )
        try:
            self.width = min(80, os.get_terminal_size().columns)
        except OSError:
            self.width = 80

        self._n_res: int
        self._data: dict[Enum, list[StructuralIrregularity]]
        self._score: int

    @property
    def n_res(self):
        try:
            return self._n_res
        except AttributeError as error:
            raise ValueError("Structure has not been processed yet.") from error

    @property
    def data(self):
        try:
            return self._data
        except AttributeError as error:
            raise ValueError("Structure has not been processed yet.") from error

    @property
    def score(self):
        try:
            return self._score
        except AttributeError as error:
            raise ValueError("Score has not been computed yet.") from error

    def process_structure(self, structure: Structure) -> None:
        """Calculate chosen measures from structure"""
        structure_data: dict[Enum, list[StructuralIrregularity]] = {}
        for func in self.funcs:
            structure_data.update(func(structure))

        self._n_res = len(list(structure.get_residues()))
        self._data = structure_data
        return None

    def output_to_terminal(self) -> None:
        """Organise the output to terminal"""
        if not any(self.data[x] for x in (AmideBonds.CIS, AmideBonds.NON_PLANAR,
                                          AmideBonds.CIS_PROLINE, ChiralCenters.D, Clashes.VDW)):
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

    def to_pml(self, structure_path: Path) -> str:
        """Transforms output of processed structure into a string which can be used by pymol to open
the structure."""
        commands = []
        if structure_path.exists():
            commands.append(f'load {structure_path}')
        else:
            commands.append(f'fetch {structure_path}')
        commands.append('set cartoon_transparency, 0.3')
        commands.append('hide lines')
        commands.append('show cartoon')
        commands.append('color white')

        for key, values in reversed(self.data.items()):
            sel_string = ' or '.join((entry.to_pymol() for entry in values))
            if sel_string:
                commands.append(f'select {key.name}, {sel_string}')
                try:
                    color = self.display[key].value
                except KeyError:
                    continue
                else:
                    commands.append(f'color {color}, {key.name}')
                    commands.append(f'show sticks, {key.name}')

        commands.append('color atomic, not elem C')
        commands.append('deselect')

        pml = '\n'.join(commands)
        return pml

    def compute_score(self) -> None:
        """Assign a score based on the errors in a Structure."""
        raw_score = len(self.data[Clashes.VDW])
        score = 30*raw_score // self.n_res
        self._score = int(score)
        return None

    def __repr__(self):
        return f"{self.__class__.__name__}(width={self.width}px)"


@click.command()
@click.argument(
        "files",
        required=True,
        nargs=-1,
        )
@click.option("--amides/--no-amides", is_flag=True, default=True)
@click.option("--chiralities/--no-chiralities", is_flag=True, default=True)
@click.option("--clashes/--no-clashes", is_flag=True, default=True)
@click.option("--score/--no-score", is_flag=True, default=False)
@click.option("--quiet/--verbose", is_flag=True, default=False, help='Start PyMOL quietly')
@click.option("--pymol", is_flag=True, default=False, show_default=True,
              help=('Open the structure in PyMOL with the irregularities annotated. '
                    'Requires pymol to be in path.')
              )
def main(files: list[str],
         amides: bool,
         chiralities: bool,
         clashes: bool,
         score: bool,
         quiet: bool,
         pymol: bool,
         ) -> None:
    """Check structure models for errors"""

    paths = (Path(file) for file in files)
    for path in paths:
        app = App(amides, chiralities, clashes)
        click.echo('-'*app.width)
        click.echo(click.style(f'{path.name.upper()}', bold=True))
        try:
            struc = get_structure(path)
        except PDBCodeError as error:
            click.echo(click.style(error, fg='white', bg='red', bold=True))
            raise click.Abort() from error
        try:
            app.process_structure(struc)
        except MissingInformationError as error:
            click.echo(click.style(error, fg='white', bg='red', bold=True))
            raise click.Abort() from error

        app.output_to_terminal()
        if score:
            app.compute_score()
            if app.score == 0:
                score_color = 'green'
            elif app.score < 10:
                score_color = 'yellow'
            else:
                score_color = 'red'
            click.echo(click.style('\nCLASH SCORE: ', bold=True) +
                       click.style(f'{app.score:>3d}', bold=True, fg=score_color))
        if pymol:
            pml = app.to_pml(path)
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pml') as tmp:
                tmp.write(pml)
                tmp.flush()
                command = 'pymol -Q' if quiet else 'pymol'
                subprocess.call((*command.split(), tmp.name))
                if not quiet:
                    click.echo(f'{"SUMMARY":=^{app.width}}')
                    app.output_to_terminal()
                    if score:
                        click.echo(click.style('\nCLASH SCORE: ', bold=True) +
                                   click.style(f'{app.score:>3d}', bold=True, fg=score_color))

    click.echo('-'*app.width)


if __name__ == '__main__':
    main()
