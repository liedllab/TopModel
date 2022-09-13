from typing import Any
from Bio.PDB.mmtf import MMTFParser as BioMMTFParser


class TopModelMMTFParser(BioMMTFParser):
    """Wrap MMTFParser from BioPython to fullfill the Parser Protocol."""
    def get_structure(self, struc_id: Any, filename: str):
        """Get a structure from a file - given a file path."""
        return super().get_structure(filename)
