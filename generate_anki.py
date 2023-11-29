#!/usr/bin/env python

import numpy as np
import pandas as pd
import requests

from anki_model import DeckSet, MoleculeNote, Package
from chembl_webresource_client.new_client import new_client
from rdkit.Chem.rdmolfiles import MolFromMolBlock, MolToMolBlock
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMolecule
from rdkit.Chem.rdDistGeom import EmbedMolecule
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolops import AddHs, RemoveHs
from rdkit.Chem.Draw import MolToImage
from rdkit.Chem.rdDepictor import Compute2DCoords, GenerateDepictionMatching3DStructure
from rdkit.rdBase import BlockLogs
from pathlib import Path
from pymol2 import PyMOL
from tqdm import tqdm
from typing import NamedTuple

class Row(NamedTuple):
    document_location: str
    classification: str
    name: str
    chemblid: str
    pubchemid: str


def fetch_3d_molecule(row: Row) -> Mol:
    if row.pubchemid:
        rsp = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{row.pubchemid}/sdf?record_type=3d"
        )
        if rsp.status_code == 200:
            mol = MolFromMolBlock(rsp.text, removeHs=False)
            mol.SetProp("_Name", row.pubchemid)

            return mol

    if row.chemblid:
        molecule = new_client.molecule.get(row.chemblid)  # type:ignore
        if molecule["molecule_structures"]:
            mol = MolFromMolBlock(
                molecule["molecule_structures"]["molfile"], removeHs=False
            )
            mol.SetProp("_Name", row.chemblid)

            try:
                with BlockLogs():
                    AddHs(mol)
                    EmbedMolecule(mol)
                    MMFFOptimizeMolecule(mol)
            except BaseException:
                pass

            return mol

    raise Exception("could not")


def fetch_2d_molecule(row: Row) -> Mol:
    if row.pubchemid:
        rsp = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{row.pubchemid}/sdf?record_type=2d"
        )
        if rsp.status_code == 200:
            mol = MolFromMolBlock(rsp.text, removeHs=True)
            mol.SetProp("_Name", row.pubchemid)

            return mol

    if row.chemblid:
        molecule = new_client.molecule.get(row.chemblid)  # type:ignore
        if molecule["molecule_structures"]:
            mol = MolFromMolBlock(
                molecule["molecule_structures"]["molfile"], removeHs=True
            )
            mol.SetProp("_Name", row.chemblid)

            return mol

    raise Exception("could not")


def write_files(dir: Path, row: Row):
    if (
        dir.joinpath(f"{row.pubchemid}_3D.png").exists()
        or dir.joinpath(f"{row.chemblid}_3D.png").exists()
    ):
        return

    try:
        mol3d = fetch_3d_molecule(row)
    except Exception:
        return

    id = mol3d.GetProp("_Name")

    with dir.joinpath(f"{id}_3D.mol").open("+w") as f:
        f.write(MolToMolBlock(mol3d))

    with PyMOL() as p:
        p.cmd.load(dir.joinpath(f"{id}_3D.mol"))
        p.cmd.color("gray", "(symbol C)")
        p.cmd.show_as("sticks")
        p.cmd.hide("(symbol H)")
        p.cmd.orient()
        p.cmd.zoom(complete=1)
        p.cmd.png(dir.joinpath(f"{id}_3D.png").as_posix(), 1000, 800, dpi=150, ray=1)

    try:
        mol2d = fetch_2d_molecule(row)
    except Exception:
        return

    try:
        GenerateDepictionMatching3DStructure(mol2d, mol3d)
    except BaseException:
        Compute2DCoords(mol2d)

    MolToImage(mol2d).save(dir.joinpath(f"{id}_2D.png"))


upstream_url = "https://docs.google.com/spreadsheets/d/e/2PACX-1vQgbPmrlTty5Q2luk79OigcbyWyQXAQR4xMpxNWJYHwMPpZvGjhBN7wd88vgAyWGMyzwIedvpR4iiNO/pub?output=xlsx"
data_dir = Path("data")
data_dir.mkdir(exist_ok=True, parents=True)

sheets = pd.read_excel(upstream_url, sheet_name=None)
decks = DeckSet("15 Pharmazeutische Chemie")

for sheet_name, sheet in sheets.items():
    if len(sheet_name) < 1 or sheet_name[0] != "D":
        continue

    sheet = sheet.replace({np.nan: None, "#NAME?": None})

    for row in tqdm(list(sheet.itertuples()), desc=sheet_name):  # type: ignore
        row: Row

        write_files(data_dir, row)

        note = MoleculeNote(
            name=row.name, chemblid=row.chemblid, pubchemid=row.pubchemid
        )
        decks.add_note(
            row.classification if pd.notnull(row.classification) else "Other", note
        )

Package(
    deck_or_decks=decks.to_list(), media_files=data_dir.glob("*.png")
).write_to_file(data_dir.joinpath(f"{decks.root_deck_name}.apkg"))
