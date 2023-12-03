from os import PathLike
from typing import List

from genanki import Deck, Model, Note, Package, guid_for


class DeckNumberer:
    def __init__(self):
        self._mapping = {}

    def number(self, name: str) -> str:
        segments = []

        for i, segment in enumerate(name.split("::")):
            parent = "::".join(segments)
            self._mapping.setdefault(parent, {})
            self._mapping[parent].setdefault(segment, len(self._mapping[parent]) + 1)

            segments.append(f"{self._mapping[parent][segment]:02d}: {segment}")

        return "::".join(segments)


class DeckSet:
    def __init__(self, root_deck_name: str):
        self.root_deck_name = root_deck_name
        self._decks = {}

    def add_note(self, subdeck_name: str, note: Note):
        full_name = f"{self.root_deck_name}::{subdeck_name}"
        if subdeck_name not in self._decks:
            self._decks[subdeck_name] = Deck(abs(hash(full_name)), full_name)

        self._decks[subdeck_name].add_note(note)

    def count(self) -> int:
        count = 0
        for deck in self._decks.values():
            count += len(deck.notes)
        return count

    def to_list(self) -> List[Deck]:
        return list(self._decks.values())


model = Model(
    1720913573,
    "Molecule",
    sort_field_index=0,
    fields=[
        {"name": "Name"},
        {"name": "ChEMBL ID"},
        {"name": "PubChem ID"},
        {"name": "2D Structure"},
        {"name": "3D Structure"},
    ],
    templates=[
        {
            "name": "3D Structure to Name",
            "qfmt": "{{3D Structure}}",
            "afmt": '{{FrontSide}}<hr id="answer">{{Name}}',
        }
    ],
    css="""
.card {
  font-family: arial;
  font-size: 20px;
  text-align: center;
  color: black;
  background-color: white;
}
    """,
)


class MoleculeNote(Note):
    def __init__(
        self,
        name: str,
        chemblid: str,
        pubchemid: str,
        file_2d: PathLike,
        file_3d: PathLike,
    ):
        super().__init__(
            model=model,
            fields=[
                name,
                chemblid,
                pubchemid,
                f"<img src='{file_2d}' />",
                f"<img src='{file_3d}' />",
            ],
        )

    @property
    def guid(self):
        return guid_for(self.fields[0] + ":::" + self.fields[1])  # type:ignore
