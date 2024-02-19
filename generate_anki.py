#!/usr/bin/env python

from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

from anki_model import DeckNumberer, DeckSet, MoleculeNote, Package
from molecule_depict import MoleculeDepictor
from spreadsheet_model import Row, SpreadsheetModel

upstream_url = "https://docs.google.com/spreadsheets/d/e/2PACX-1vQgbPmrlTty5Q2luk79OigcbyWyQXAQR4xMpxNWJYHwMPpZvGjhBN7wd88vgAyWGMyzwIedvpR4iiNO/pub?output=xlsx"
data_dir = Path("data")
data_dir.mkdir(exist_ok=True, parents=True)

sheets = pd.read_excel(upstream_url, sheet_name=None, dtype=str, keep_default_na=False)
decks = DeckSet("B15 Pharmazeutische Chemie")
numberer = DeckNumberer()
depictor = MoleculeDepictor(data_dir)

for sheet_name, sheet in sheets.items():
    if len(sheet_name) < 1 or sheet_name[0] != "#":
        continue

    try:
        sheet = SpreadsheetModel.validate(sheet).replace({np.nan: None})
    except BaseException as e:
        raise Exception(f"validating {sheet_name}") from e

    for row in tqdm(list(sheet.itertuples()), desc=sheet_name):  # type: ignore
        row: Row

        if row.skip:
            continue

        output = depictor.depict(row.pubchemid or row.chemblid)

        note = MoleculeNote(
            name=row.name,
            taxonomy=row.classification.split("::"),
            chemblid=f"{row.chemblid}",
            pubchemid=f"{row.pubchemid}",
            file_2d=output.image2d,
            file_3d=output.image3d,
        )

        decks.add_note(numberer.number(row.classification), note)

print(f"Total molecules: {decks.count()}")

Package(
    deck_or_decks=decks.to_list(), media_files=data_dir.glob("*.png")
).write_to_file(data_dir.joinpath(f"{decks.root_deck_name}.apkg"))
