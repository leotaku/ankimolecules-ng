import pandas as pd
from anki_model import DeckSet, MoleculeNote, Package
from chembl_webresource_client.new_client import new_client
from pathlib import Path
from pymol2 import PyMOL
from tqdm import tqdm
from typing import NamedTuple
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMolecule
from rdkit.Chem import MolFromMolBlock, MolToMolBlock, AddHs  # type:ignore
from rdkit.rdBase import BlockLogs


class Row(NamedTuple):
    document_location: str
    classification: str
    name: str
    chemblid: str


def write_required_files(dir: Path, chemblid: str):
    molecule = new_client.molecule.get(chemblid)  # type:ignore
    if not molecule["molecule_structures"]:
        tqdm.write(f"skipping {chemblid}")
        return

    with dir.joinpath(f"{chemblid}_3D.mol").open("+w") as f:
        mol = MolFromMolBlock(
            molecule["molecule_structures"]["molfile"], removeHs=False
        )

        # with BlockLogs():
        #     try:
        #         mol_with_hs = AddHs(mol)
        #         MMFFOptimizeMolecule(mol_with_hs)
        #     except:
        #         try:
        #             MMFFOptimizeMolecule(mol)
        #         except:
        #             pass
        #     else:
        #         mol = mol_with_hs

        f.write(MolToMolBlock(mol))

    with PyMOL() as p:
        p.cmd.load(dir.joinpath(f"{chemblid}_3D.mol"))
        p.cmd.color("gray", "(symbol C)")
        p.cmd.show_as("sticks")
        p.cmd.orient()
        p.cmd.zoom(complete=1)
        p.cmd.png(dir.joinpath(f"{chemblid}.png").as_posix(), 1000, 800, dpi=150, ray=1)


upstream_url = "https://docs.google.com/spreadsheets/d/e/2PACX-1vQgbPmrlTty5Q2luk79OigcbyWyQXAQR4xMpxNWJYHwMPpZvGjhBN7wd88vgAyWGMyzwIedvpR4iiNO/pub?output=xlsx"
data_dir = Path("data")
data_dir.mkdir(exist_ok=True, parents=True)

sheets = pd.read_excel(upstream_url, sheet_name=None)
decks = DeckSet("15 Pharmazeutische Chemie")

for sheet_name, sheet in sheets.items():
    for row in tqdm(list(sheet.itertuples()), desc=sheet_name):  # type: ignore
        row: Row

        if pd.isnull(row.chemblid):
            continue

        write_required_files(data_dir, row.chemblid)
        note = MoleculeNote(name=row.name, chemblid=row.chemblid)
        decks.add_note(
            row.classification if pd.notnull(row.classification) else "Other", note
        )

Package(
    deck_or_decks=decks.to_list(), media_files=data_dir.glob("*.png")
).write_to_file(data_dir.joinpath(f"{decks.root_deck_name}.apkg"))
