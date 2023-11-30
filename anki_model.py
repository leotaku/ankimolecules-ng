from typing import List
from genanki import Note, Deck, Model, Package, guid_for
from os import PathLike


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

    def to_list(self) -> List[Deck]:
        return list(self._decks.values())


model = Model(
    1720913576,
    "Molecule",
    sort_field_index=1,
    fields=[
        {"name": "ChEMBL ID"},
        {"name": "PubChem ID"},
        {"name": "Name"},
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
    def __init__(self, name: str, chemblid: str, pubchemid: str):
        super().__init__(
            model=model,
            fields=[chemblid, pubchemid, name, f"<img src='{chemblid}.png' />"],
        )

    @property
    def guid(self):
        return guid_for(self.fields[0])  # type:ignore
