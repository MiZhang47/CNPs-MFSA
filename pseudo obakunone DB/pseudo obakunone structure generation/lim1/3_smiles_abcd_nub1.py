import pandas as pd
from rdkit import Chem
from rdkit.Chem import GetSymmSSSR, Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
import os

standard_numbers_dict = {
    "A": [1, 2, 3, 4, 5, 6],
    "B": [3, 2, 7, 8, 9, 10],
    "C": [9, 10, 11, 12, 13, 14],
    "D": [13, 14, 15, 16, 17]
}

base_path = r"C:\Users\xxx\Desktop\D_SMILES_x215\Unique_Rings"
input_files = {
    "A": os.path.join(base_path, "A_Ring_Unique.csv"),
    "B": os.path.join(base_path, "B_Ring_Unique.csv"),
    "C": os.path.join(base_path, "C_Ring_Unique.csv"),
    "D": os.path.join(base_path, "D_Ring_Unique.csv")
}
output_img_folder = os.path.join(base_path, "Ring_Atom_Images")
os.makedirs(output_img_folder, exist_ok=True)

def get_ring_atoms(smiles, ring_type):
    mol = Chem.MolFromSmiles(smiles)
    rings = [list(r) for r in GetSymmSSSR(mol)]

    if ring_type == "D":
        five_membered = [r for r in rings if len(r) == 5]
        for r in five_membered:
            if all(mol.GetAtomWithIdx(i).GetSymbol() == "C" for i in r):
                return r
        six_membered = [r for r in rings if len(r) == 6]
        for r in six_membered:
            atom_syms = [mol.GetAtomWithIdx(i).GetSymbol() for i in r]
            if atom_syms.count("O") == 1 and atom_syms.count("C") == 5:
                return r
        return []
    else:
        carbon_rings = [r for r in rings if all(mol.GetAtomWithIdx(i).GetSymbol() == "C" for i in r)]
        candidates = [r for r in carbon_rings if 5 <= len(r) <= 6]
        return max(candidates, key=len) if candidates else []

def assign_standard_numbers(ring_type, atom_indices, mol):
    std_nums = standard_numbers_dict[ring_type]
    if ring_type == "D":
        atom_indices = [idx for idx in atom_indices if mol.GetAtomWithIdx(idx).GetSymbol() == "C"][:5]
        std_nums = std_nums[:len(atom_indices)]
    for i in range(len(atom_indices)):
        idx = atom_indices[i]
        atom = mol.GetAtomWithIdx(idx)
        branches = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in atom_indices]
        if any(nbr.GetSymbol() == "C" and nbr.GetDegree() == 1 for nbr in branches):
            reordered = atom_indices[i:] + atom_indices[:i]
            return reordered[:len(std_nums)], std_nums
    return atom_indices[:len(std_nums)], std_nums

def set_atommap_numbers(mol, index_map):
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in index_map:
            atom.SetAtomMapNum(index_map[idx])
    return mol

def mol_to_mapped_smiles(mol):
    return Chem.MolToSmiles(mol, isomericSmiles=False)

def get_canonical_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol, canonical=True)
    except:
        return None
    return None

def draw_numbered_molecule(smiles, atom_indices, output_path, label_text=None):
    from PIL import Image, ImageDraw
    import io

    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DCairo(400, 400)
    options = drawer.drawOptions()
    for idx in atom_indices:
        options.atomLabels[idx] = str(idx)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    img_bytes = drawer.GetDrawingText()

    if label_text:
        img = Image.open(io.BytesIO(img_bytes))
        draw = ImageDraw.Draw(img)
        draw.text((10, 10), label_text, fill=(0, 0, 0))
        img.save(output_path)
    else:
        with open(output_path, "wb") as f:
            f.write(img_bytes)

for ring_type, file_path in input_files.items():
    df = pd.read_csv(file_path)
    new_indices_list = []
    std_numbers_list = []
    mapped_smiles_list = []
    canonical_list = []

    canonical_seen = {}

    for i, row in df.iterrows():
        smiles = row["Ring_SMILES"]
        mol = Chem.MolFromSmiles(smiles)
        ring_atoms = get_ring_atoms(smiles, ring_type)

        if not ring_atoms:
            new_indices_list.append("Error")
            std_numbers_list.append("Error")
            mapped_smiles_list.append("Error")
            canonical_list.append("Error")
            continue

        reordered, std = assign_standard_numbers(ring_type, ring_atoms, mol)
        new_indices_list.append(reordered)
        std_numbers_list.append(std)

        index_map = dict(zip(reordered, std))
        mol_mapped = set_atommap_numbers(Chem.MolFromSmiles(smiles), index_map)
        mapped_smiles = mol_to_mapped_smiles(mol_mapped)
        mapped_smiles_list.append(mapped_smiles)

        canonical = get_canonical_smiles(smiles)
        canonical_list.append(canonical)

    df["Ring_Atom_Indices"] = new_indices_list
    df["Ring_Assigned_Numbers"] = std_numbers_list
    df["Mapped_SMILES"] = mapped_smiles_list
    df["Canonical_SMILES"] = canonical_list

    df_unique = df[df["Canonical_SMILES"] != "Error"].drop_duplicates(subset=["Canonical_SMILES"]).copy()
    df_unique["Label"] = [f"{ring_type}{i+1}" for i in range(len(df_unique))]

    for i, row in df_unique.iterrows():
        label = row["Label"]
        smiles = row["Ring_SMILES"]
        indices = row["Ring_Atom_Indices"]
        if isinstance(indices, str) or indices == "Error":
            continue
        output_path = os.path.join(output_img_folder, f"{label}_Ring.png")
        draw_numbered_molecule(smiles, eval(str(indices)), output_path, label_text=label)

    output_path = os.path.join(base_path, f"{ring_type}_Ring_Unique_Renumbered.csv")
    df_unique.to_csv(output_path, index=False, encoding="utf-8-sig")
