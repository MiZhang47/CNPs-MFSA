import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import GetSymmSSSR, Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

standard_numbers_dict = {
    "A": [1, 2, 3, 4, 5, 6],
    "B": [2, 1, 6, 7, 8, 9, 10, 11],
    "C": [8, 9, 12, 13, 14, 15]
}

base_path = r"C:\Users\xxx\Desktop\D_SMILES_x215\Unique_Rings"
input_files = {
    "A": os.path.join(base_path, "A_Ring_Unique.csv"),
    "B": os.path.join(base_path, "B_Ring_Unique.csv"),
    "C": os.path.join(base_path, "C_Ring_Unique.csv")
}

output_img_folder = os.path.join(base_path, "Ring_Atom_Images")
os.makedirs(output_img_folder, exist_ok=True)

def get_ring_atoms(smiles, ring_type):
    mol = Chem.MolFromSmiles(smiles)
    rings = [list(r) for r in GetSymmSSSR(mol)]

    if ring_type == "A":
        candidates = [r for r in rings if len(r) == 6]
        for r in candidates:
            for idx in r:
                atom = mol.GetAtomWithIdx(idx)
                methyls = [n for n in atom.GetNeighbors() if n.GetSymbol() == "C" and n.GetDegree() == 1]
                if len(methyls) >= 2:
                    return r
        return []
    elif ring_type == "B":
        candidates = [r for r in rings if len(r) == 8 and all(mol.GetAtomWithIdx(i).GetSymbol() == "C" for i in r)]
        return candidates[0] if candidates else []
    elif ring_type == "C":
        candidates = [r for r in rings if len(r) == 6 and all(mol.GetAtomWithIdx(i).GetSymbol() == "C" for i in r)]
        return candidates[0] if candidates else []
    return []

def assign_standard_numbers(ring_type, atom_indices, mol):
    std_nums = standard_numbers_dict[ring_type]
    atom_indices = atom_indices[:len(std_nums)]
    std_nums = std_nums[:len(atom_indices)]
    for i in range(len(atom_indices)):
        idx = atom_indices[i]
        atom = mol.GetAtomWithIdx(idx)
        if ring_type == "A":
            methyls = [n for n in atom.GetNeighbors() if n.GetSymbol() == "C" and n.GetDegree() == 1]
            if len(methyls) >= 2:
                reordered = atom_indices[i:] + atom_indices[:i]
                return reordered[:len(std_nums)], std_nums
    return atom_indices, std_nums

def set_atommap_numbers(mol, index_map):
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in index_map:
            atom.SetAtomMapNum(index_map[idx])
    return mol

def draw_numbered_molecule(smiles, atom_indices, output_path):
    mol = Chem.MolFromSmiles(smiles)
    Chem.RemoveStereochemistry(mol)
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
    cleaned_smiles_list = []

    for i, row in df.iterrows():
        original_smiles = row["Ring_SMILES"]
        mol = Chem.MolFromSmiles(original_smiles)
        if mol is None:
            new_indices_list.append("Error")
            std_numbers_list.append("Error")
            mapped_smiles_list.append("Error")
            cleaned_smiles_list.append("Error")
            continue

        Chem.RemoveStereochemistry(mol)
        clean_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
        cleaned_smiles_list.append(clean_smiles)

        ring_atoms = get_ring_atoms(clean_smiles, ring_type)
        if not ring_atoms:
            new_indices_list.append("Error")
            std_numbers_list.append("Error")
            mapped_smiles_list.append("Error")
            continue

        reordered, std = assign_standard_numbers(ring_type, ring_atoms, mol)
        new_indices_list.append(reordered)
        std_numbers_list.append(std)

        mol_mapped = Chem.MolFromSmiles(clean_smiles)
        mol_mapped = set_atommap_numbers(mol_mapped, dict(zip(reordered, std)))
        Chem.RemoveStereochemistry(mol_mapped)
        mapped_smiles = Chem.MolToSmiles(mol_mapped, isomericSmiles=False)
        mapped_smiles_list.append(mapped_smiles)

    df["Ring_SMILES"] = cleaned_smiles_list
    df["Ring_Atom_Indices"] = new_indices_list
    df["Ring_Assigned_Numbers"] = std_numbers_list
    df["Mapped_SMILES"] = mapped_smiles_list

    df_unique = df[df["Mapped_SMILES"] != "Error"].drop_duplicates(subset=["Mapped_SMILES"]).copy()
    df_unique.reset_index(drop=True, inplace=True)
    df_unique["Label"] = [f"{ring_type}{i+1}" for i in range(len(df_unique))]

    for _, row in df_unique.iterrows():
        try:
            indices = eval(str(row["Ring_Atom_Indices"]))
            label = row["Label"]
            path = os.path.join(output_img_folder, f"{label}_Ring.png")
            draw_numbered_molecule(row["Ring_SMILES"], indices, path)
        except:
            continue

    output_csv = os.path.join(base_path, f"{ring_type}_Ring_Unique_Renumbered.csv")
    df_unique.to_csv(output_csv, index=False, encoding='utf-8-sig')
