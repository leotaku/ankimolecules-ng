from typing import List
from genanki import Note, Deck, Model, Package, guid_for


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
