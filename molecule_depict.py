from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Literal

import requests
from pymol2 import PyMOL
from rdkit.Chem.Draw import MolToImage
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdDepictor import (
    ConstrainedDepictionParams,
    GenerateDepictionMatching2DStructure,
)
from rdkit.Chem.rdDistGeom import EmbedMolecule
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMolecule
from rdkit.Chem.rdmolfiles import MolFromMolBlock, MolToMolBlock
from rdkit.Chem.rdmolops import AddHs
from rdkit.rdBase import BlockLogs
from psico.mcsalign import mcsalign


@dataclass
class MoleculeExport:
    image2d: str
    image3d: str


class MoleculeDepictor:
    def __init__(
        self,
        output_directory: PathLike,
        cache_mol: bool = True,
        cache_img: bool = True,
    ):
        self.output_directory = Path(output_directory)
        self.cache_mol = cache_mol
        self.cache_img = cache_img

        self.previous_mol_2d = None
        self.previous_mol_3d_path = None

    def depict(self, id: str) -> MoleculeExport:
        return MoleculeExport(
            image2d=self._depict_2d(id).name,
            image3d=self._depict_3d(id).name,
        )

    def _fetch_mol_cached(self, id, kind: Literal["2d", "3d"]) -> (Mol, Path):
        molfile = self.output_directory.joinpath(f"{id}_{kind.upper()}.mol")

        if self.cache_mol and molfile.exists():
            mol = MolFromMolBlock(molfile.read_bytes())
        else:
            mol = fetch_mol(id, kind)
            molfile.write_text(MolToMolBlock(mol))

        return mol, molfile

    def _depict_2d(self, id: str) -> Path:
        imagepath = self.output_directory.joinpath(f"{id}_2D.png")
        if self.cache_img and imagepath.exists():
            return imagepath

        mol, _ = self._fetch_mol_cached(id, kind="2d")

        from rdkit.Chem.rdMolAlign import AlignMol

        if self.previous_mol_2d:
            try:
                AlignMol(mol, self.previous_mol_2d)
            except:
                pass

        MolToImage(mol).save(imagepath)

        self.previous_mol_2d = mol
        return imagepath

    def _depict_3d(self, id: str) -> Path:
        imagepath = self.output_directory.joinpath(f"{id}_3D.png")
        if self.cache_img and imagepath.exists():
            return imagepath

        mol, molpath = self._fetch_mol_cached(id, kind="3d")

        with PyMOL() as p:
            p.cmd.load(molpath.as_posix(), "this")
            if self.previous_mol_3d_path:
                p.cmd.load(self.previous_mol_3d_path.as_posix(), "other")
                try:
                    mcsalign("%this", "%other", _self=p.cmd)
                except:
                    pass
                finally:
                    self.previous_mol_3d_path = molpath
                p.cmd.hide("%other")

            p.cmd.color("gray", "(symbol C)")
            p.cmd.show_as("sticks", "bonded")
            p.cmd.show_as("nb_spheres", "not bonded")
            p.cmd.hide("(symbol H)")
            p.cmd.orient()
            p.cmd.zoom(complete=1)
            p.cmd.png(imagepath.as_posix(), 1000, 800, dpi=150, ray=1)

        return imagepath


def extract_mcs_as_mol(mol: Mol, other: Mol) -> Mol | None:
    from rdkit.Chem.rdFMCS import FindMCS
    from rdkit.Chem.rdmolfiles import MolFromSmarts
    from rdkit.Chem.rdchem import RWMol

    mcs = FindMCS([mol, other])
    print(mcs)
    if not mcs:
        return None

    substructure = MolFromSmarts(mcs.smartsString)
    match = mol.GetSubstructMatch(substructure)
    if not match:
        return None

    new_mol = RWMol(mol)
    new_mol.BeginBatchEdit()
    for atom in mol.GetAtoms():
        if (id := atom.GetIdx()) not in match:
            new_mol.RemoveAtom(id)
    new_mol.CommitBatchEdit()

    return new_mol.GetMol()


def optimize_mol(mol: Mol) -> Mol:
    try:
        with BlockLogs():
            AddHs(mol)
            EmbedMolecule(mol)
            MMFFOptimizeMolecule(mol)
    except BaseException:
        pass

    return mol


def fetch_mol(id, kind: Literal["2d", "3d"]) -> Mol:
    if id.startswith("CHEMBL"):
        return fetch_chembl(id, kind)
    else:
        return fetch_pubchem(id, kind)


def fetch_chembl(id, kind: Literal["2d", "3d"]) -> Mol:
    rsp = requests.get(
        f"https://www.ebi.ac.uk/chembl/api/data/molecule/{id}?format=sdf",
    )
    if rsp.status_code != 200:
        raise Exception(f"could not fetch from ChEMBL: '{id}'")

    mol = MolFromMolBlock(rsp.text, removeHs=True)
    mol.SetProp("_Name", id)

    if kind == "3d":
        return optimize_mol(mol)

    return mol


def fetch_pubchem(id, kind: Literal["2d", "3d"]) -> Mol:
    rsp = requests.get(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{id}/sdf?record_type={kind}"
    )
    if rsp.status_code != 200:
        if kind == "3d":
            mol = fetch_pubchem(id, kind="2d")
            return optimize_mol(mol)
        else:
            raise Exception(f"could not fetch from PubChem: '{id}'")

    mol = MolFromMolBlock(rsp.text, removeHs=True)
    mol.SetProp("_Name", id)

    return mol
