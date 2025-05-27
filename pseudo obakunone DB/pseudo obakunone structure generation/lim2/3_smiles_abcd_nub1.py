import pandas as pd
from rdkit import Chem
from rdkit.Chem import GetSymmSSSR, Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
import os

standard_numbers_dict = {
    "A": [1, 5, 6, 4, 3, 2],
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

    if ring_type == "A":
        seven = [r for r in rings if len(r) == 7]
        for r in seven:
            syms = [mol.GetAtomWithIdx(i).GetSymbol() for i in r]
            if syms.count("O") == 1 and syms.count("C") == 6:
                return r
        return []

    elif ring_type != "D":
        carbon_rings = [r for r in rings if all(mol.GetAtomWithIdx(i).GetSymbol() == "C" for i in r)]
        candidates = [r for r in carbon_rings if 5 <= len(r) <= 6]
        return max(candidates, key=len) if candidates else []

    else:
        six = [r for r in rings if len(r) == 6]
        for r in six:
            syms = [mol.GetAtomWithIdx(i).GetSymbol() for i in r]
            if syms.count("O") == 1 and syms.count("C") == 5:
                return r
        return []

def assign_standard_numbers(ring_type, atom_indices, mol):
    std_nums = standard_numbers_dict[ring_type]
    if ring_type in ["A", "D"]:
        atom_indices = [idx for idx in atom_indices if mol.GetAtomWithIdx(idx).GetSymbol() == "C"]
        atom_indices = atom_indices[:len(std_nums)]
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

def draw_numbered_molecule(smiles, atom_indices, output_path, label=None):
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DCairo(400, 400)
    options = drawer.drawOptions()
    for idx in atom_indices:
        options.atomLabels[idx] = str(idx)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    with open(output_path, "wb") as f:
        f.write(drawer.GetDrawingText())

for ring_type, file_path in input_files.items():
    df = pd.read_csv(file_path)
    new_indices_list = []
    std_numbers_list = []
    mapped_smiles_list = []
    canonical_smiles_list = []

    for i, row in df.iterrows():
        smiles = row["Ring_SMILES"]
        mol = Chem.MolFromSmiles(smiles)
        ring_atoms = get_ring_atoms(smiles, ring_type)

        if not ring_atoms:
            new_indices_list.append("Error")
            std_numbers_list.append("Error")
            mapped_smiles_list.append("Error")
            canonical_smiles_list.append("Error")
            continue

        reordered, std = assign_standard_numbers(ring_type, ring_atoms, mol)
        new_indices_list.append(reordered)
        std_numbers_list.append(std)

        index_map = dict(zip(reordered, std))
        mol_mapped = set_atommap_numbers(Chem.MolFromSmiles(smiles), index_map)
        mapped_smiles = mol_to_mapped_smiles(mol_mapped)
        canonical = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), canonical=True)

        mapped_smiles_list.append(mapped_smiles)
        canonical_smiles_list.append(canonical)

    df["Ring_Atom_Indices"] = new_indices_list
    df["Ring_Assigned_Numbers"] = std_numbers_list
    df["Mapped_SMILES"] = mapped_smiles_list
    df["Canonical_SMILES"] = canonical_smiles_list

    df_unique = df[df["Mapped_SMILES"] != "Error"].drop_duplicates(subset=["Mapped_SMILES"]).copy()
    df_unique.reset_index(drop=True, inplace=True)
    df_unique["Label"] = [f"{ring_type}{i+1}" for i in range(len(df_unique))]

    for _, row in df_unique.iterrows():
        try:
            indices = eval(str(row["Ring_Atom_Indices"]))
            label = row["Label"]
            path = os.path.join(output_img_folder, f"{label}_Ring.png")
            draw_numbered_molecule(row["Ring_SMILES"], indices, path, label=label)
        except:
            continue

    output_csv = os.path.join(base_path, f"{ring_type}_Ring_Unique_Renumbered.csv")
    df_unique.to_csv(output_csv, index=False, encoding='utf-8-sig')
