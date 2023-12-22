#!/usr/bin/env python

from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import requests
from chembl_webresource_client.new_client import new_client
from pymol2 import PyMOL
from rdkit.Chem.Draw import MolToImage
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdDistGeom import EmbedMolecule
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMolecule
from rdkit.Chem.rdmolfiles import MolFromMolBlock, MolToMolBlock
from rdkit.Chem.rdmolops import AddHs
from rdkit.rdBase import BlockLogs
from tqdm import tqdm

from anki_model import DeckNumberer, DeckSet, MoleculeNote, Package
from spreadsheet_model import Row, SpreadsheetModel


def optimize_mol(mol: Mol) -> Mol:
    try:
        with BlockLogs():
            AddHs(mol)
            EmbedMolecule(mol)
            MMFFOptimizeMolecule(mol)
    except BaseException:
        pass

    return mol


def fetch_chembl(id, kind: Literal["2d", "3d"]) -> Mol:
    rsp = new_client.molecule.get(row.chemblid)  # type:ignore
    mol = MolFromMolBlock(rsp["molecule_structures"]["molfile"], removeHs=False)
    mol.SetProp("_Name", row.chemblid)

    if kind == "3d":
        return optimize_mol(mol)

    return mol


def fetch_pubchem(id, kind: Literal["2d", "3d"]) -> Mol:
    rsp = requests.get(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{row.pubchemid}/sdf?record_type={kind}"
    )
    if rsp.status_code != 200:
        if kind == "3d":
            mol = fetch_pubchem(id, kind="2d")
            return optimize_mol(mol)
        else:
            raise Exception(f"could not fetch from PubChem: {row.pubchemid}")

    mol = MolFromMolBlock(rsp.text, removeHs=False)
    mol.SetProp("_Name", row.pubchemid)

    return mol


def fetch_mol(row: Row, kind: Literal["2d", "3d"]) -> Mol:
    if row.pubchemid:
        return fetch_pubchem(row.pubchemid, kind)

    if row.chemblid:
        return fetch_chembl(row.chemblid, kind)

    raise Exception(f"no valid ID present for: {row.name}")


def files_exist(dir: Path, row: Row):
    result = {}

    for id in [row.chemblid, row.pubchemid]:
        if not id:
            continue

        path3d = dir.joinpath(id + "_3D.png")
        if path3d.exists():
            result["3d"] = path3d.name

        path2d = dir.joinpath(id + "_2D.png")
        if path2d.exists():
            result["2d"] = path2d.name

    return result if len(result) == 2 else None


def write_files(dir: Path, row: Row):
    mol3d = fetch_mol(row, kind="3d")
    pathMol3d = dir.joinpath(f"{mol3d.GetProp('_Name')}_3D.mol")
    path3d = dir.joinpath(f"{mol3d.GetProp('_Name')}_3D.png")

    with pathMol3d.open("+w") as f:
        f.write(MolToMolBlock(mol3d))

    with PyMOL() as p:
        p.cmd.load(pathMol3d)
        p.cmd.color("gray", "(symbol C)")
        p.cmd.show_as("sticks")
        p.cmd.hide("(symbol H)")
        p.cmd.orient()
        p.cmd.zoom(complete=1)
        p.cmd.png(path3d.as_posix(), 1000, 800, dpi=150, ray=1)

    mol2d = fetch_mol(row, kind="2d")
    path2d = dir.joinpath(f"{mol2d.GetProp('_Name')}_2D.png")

    MolToImage(mol2d).save(path2d)

    return {"2d": path2d.name, "3d": path3d.name}


upstream_url = "https://docs.google.com/spreadsheets/d/e/2PACX-1vQgbPmrlTty5Q2luk79OigcbyWyQXAQR4xMpxNWJYHwMPpZvGjhBN7wd88vgAyWGMyzwIedvpR4iiNO/pub?output=xlsx"
data_dir = Path("data")
data_dir.mkdir(exist_ok=True, parents=True)

sheets = pd.read_excel(upstream_url, sheet_name=None, dtype=str)
decks = DeckSet("B15 Pharmazeutische Chemie")
numberer = DeckNumberer()

for sheet_name, sheet in sheets.items():
    if len(sheet_name) < 1 or sheet_name[0] != "#":
        continue

    try:
        sheet = sheet.replace({np.nan: None, "#NAME?": None})
        sheet.document_location = sheet.document_location.replace({None: "N/A"})
        sheet = SpreadsheetModel.validate(sheet)
    except BaseException as e:
        raise Exception(f"validating {sheet_name}") from e

    for row in tqdm(list(sheet.itertuples()), desc=sheet_name):  # type: ignore
        row: Row

        if row.skip:
            continue

        if not (files := files_exist(data_dir, row)):
            files = write_files(data_dir, row)

        note = MoleculeNote(
            name=row.name,
            taxonomy=row.classification.split("::"),
            chemblid=f"{row.chemblid}",
            pubchemid=f"{row.pubchemid}",
            file_2d=files["2d"],
            file_3d=files["3d"],
        )

        decks.add_note(numberer.number(row.classification), note)

print(f"Total molecules: {decks.count()}")

Package(
    deck_or_decks=decks.to_list(), media_files=data_dir.glob("*.png")
).write_to_file(data_dir.joinpath(f"{decks.root_deck_name}.apkg"))
