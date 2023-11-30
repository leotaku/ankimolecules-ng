from typing import NamedTuple

import pandera as pa
from pandera.typing import Series


class Row(NamedTuple):
    document_location: str
    classification: str
    name: str
    chemblid: str | None
    pubchemid: str | None
    skip: str | None


class SpreadsheetModel(pa.DataFrameModel):
    document_location: Series[str] = pa.Field()
    classification: Series[str] = pa.Field()
    name: Series[str] = pa.Field()
    chemblid: Series[str] = pa.Field(nullable=True)
    pubchemid: Series[str] = pa.Field(nullable=True)
    skip: Series[str] = pa.Field(nullable=True)

    class Config:
        add_missing_columns = True
