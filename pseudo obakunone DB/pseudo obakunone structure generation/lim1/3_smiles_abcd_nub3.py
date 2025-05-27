import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

base_path = r"C:\Users\xxx\Desktop\D_SMILES_x215\Unique_Rings"
file_ab = os.path.join(base_path, "Fused_AB_SMILES.csv")
file_cd = os.path.join(base_path, "Fused_CD_SMILES.csv")
output_csv = os.path.join(base_path, "Fused_ABCD_SMILES.csv")
output_img_dir = os.path.join(base_path, "Fused_ABCD_Images")
os.makedirs(output_img_dir, exist_ok=True)

def get_atommap_mapping(mol):
    return {atom.GetAtomMapNum(): atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum() in [9, 10]}

def fuse_by_common_edge(mol1, mol2, map1, map2, label1, label2):
    idx1_a = map1[label1]
    idx1_b = map1[label2]
    idx2_a = map2[label1]
    idx2_b = map2[label2]
    offset = mol1.GetNumAtoms()

    def get_neighbors(mol, src, exclude):
        nbs, bts = [], []
        for nbr in mol.GetAtomWithIdx(src).GetNeighbors():
            n_idx = nbr.GetIdx()
            if n_idx == exclude:
                continue
            bts.append(mol.GetBondBetweenAtoms(src, n_idx).GetBondType())
            nbs.append(n_idx + offset)
        return nbs, bts

    nbs_a, bts_a = get_neighbors(mol2, idx2_a, idx2_b)
    nbs_b, bts_b = get_neighbors(mol2, idx2_b, idx2_a)

    combo = Chem.CombineMols(mol1, mol2)
    em = Chem.EditableMol(combo)

    remove = sorted([idx2_a + offset, idx2_b + offset], reverse=True)
    for i in remove:
        em.RemoveAtom(i)

    def shift(i): return i - sum(1 for r in remove if i > r)
    nbs_a = [shift(i) for i in nbs_a]
    nbs_b = [shift(i) for i in nbs_b]

    for nb, bt in zip(nbs_a, bts_a):
        em.AddBond(idx1_a, nb, bt)
    for nb, bt in zip(nbs_b, bts_b):
        em.AddBond(idx1_b, nb, bt)

    mol_fused = em.GetMol()
    Chem.SanitizeMol(mol_fused, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
    return mol_fused

df_ab = pd.read_csv(file_ab)
df_cd = pd.read_csv(file_cd)

fused_data = []
count = 0

for i, row_ab in df_ab.iterrows():
    for j, row_cd in df_cd.iterrows():
        try:
            smi_ab = row_ab["Fused_AB_SMILES"]
            smi_cd = row_cd["Fused_CD_SMILES"]

            mol_ab = Chem.MolFromSmiles(smi_ab)
            mol_cd = Chem.MolFromSmiles(smi_cd)

            map_ab = get_atommap_mapping(mol_ab)
            map_cd = get_atommap_mapping(mol_cd)

            mol_fused = fuse_by_common_edge(mol_ab, mol_cd, map_ab, map_cd, 9, 10)
            fused_smi = Chem.MolToSmiles(mol_fused, isomericSmiles=True)

            img = Draw.MolToImage(mol_fused, size=(400, 400))
            img.save(os.path.join(output_img_dir, f"Fused_ABCD_{i}_{j}.png"))

            fused_data.append({
                "A_Label": row_ab.get("A_Label", f"A{i}"),
                "B_Label": row_ab.get("B_Label", f"B{i}"),
                "C_Label": row_cd.get("C_Label", f"C{j}"),
                "D_Label": row_cd.get("D_Label", f"D{j}"),
                "A_SMILES": row_ab.get("A_SMILES", ""),
                "A_Atom_Indices": row_ab.get("A_Atom_Indices", ""),
                "A_Assigned_Numbers": row_ab.get("A_Assigned_Numbers", ""),
                "B_SMILES": row_ab.get("B_SMILES", ""),
                "B_Atom_Indices": row_ab.get("B_Atom_Indices", ""),
                "B_Assigned_Numbers": row_ab.get("B_Assigned_Numbers", ""),
                "C_SMILES": row_cd.get("C_SMILES", ""),
                "C_Atom_Indices": row_cd.get("C_Atom_Indices", ""),
                "C_Assigned_Numbers": row_cd.get("C_Assigned_Numbers", ""),
                "D_SMILES": row_cd.get("D_SMILES", ""),
                "D_Atom_Indices": row_cd.get("D_Atom_Indices", ""),
                "D_Assigned_Numbers": row_cd.get("D_Assigned_Numbers", ""),
                "Fused_AB_SMILES": smi_ab,
                "Fused_CD_SMILES": smi_cd,
                "Fused_ABCD_SMILES": fused_smi
            })
            count += 1
        except:
            continue

df_out = pd.DataFrame(fused_data)
df_out.to_csv(output_csv, index=False)
