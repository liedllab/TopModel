import warnings

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import click

from pdbear.chirality import get_chiralities
from pdbear.amide_bond import get_stereo

@click.command()
@click.option("--path", help="Path to the PDB file")
def check_structure(path: str) -> None:
    parser = PDBParser()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=PDBConstructionWarning)
        pdb = parser.get_structure("", path)

    chiralities = ((n, chiral_info) for n, chiral_info in enumerate(get_chiralities(pdb), 1) \
                    if chiral_info == 'D')
    stereo = ((n, db_info) for n, db_info in enumerate(get_stereo(pdb), 1) \
                    if db_info != "trans")

    print("Errors in chirality: ", end='')
    for n, _ in chiralities:
        print(f"res {n}", end=', ')


    print("\nCis amides: ", end='')
    for n, _ in chiralities:
        print(f"res {n}-{n+1}", end=', ')

    print("\nFinished.")

if __name__ == '__main__':
    check_structure()
