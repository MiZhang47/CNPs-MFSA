import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

base_path = r"C:\Users\xxx\Desktop\D_SMILES_x215\Unique_Rings"
file_a = os.path.join(base_path, "A_Ring_Unique_Renumbered.csv")
file_b = os.path.join(base_path, "B_Ring_Unique_Renumbered.csv")
file_c = os.path.join(base_path, "C_Ring_Unique_Renumbered.csv")
file_d = os.path.join(base_path, "D_Ring_Unique_Renumbered.csv")

output_ab_csv = os.path.join(base_path, "Fused_AB_SMILES.csv")
output_cd_csv = os.path.join(base_path, "Fused_CD_SMILES.csv")
output_img_folder = os.path.join(base_path, "Fused_Images")
os.makedirs(output_img_folder, exist_ok=True)

def color_atoms(smiles, atom_indices_dict, color_map):
    mol = Chem.MolFromSmiles(smiles)
    for label, isotope in color_map.items():
        idx = atom_indices_dict.get(label)
        if idx is not None:
            mol.GetAtomWithIdx(idx).SetIsotope(isotope)
    return mol

def mark_colored_atoms(mol, color_map):
    inverse_map = {v: k for k, v in color_map.items()}
    for atom in mol.GetAtoms():
        isotope = atom.GetIsotope()
        if isotope in inverse_map:
            atom.SetAtomMapNum(inverse_map[isotope])
        atom.SetIsotope(0)
    return mol

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

    fused_tmp = em.GetMol()
    for nb, bt in zip(nbs_a, bts_a):
        if not fused_tmp.GetBondBetweenAtoms(idx1_a, nb):
            em.AddBond(idx1_a, nb, bt)
    for nb, bt in zip(nbs_b, bts_b):
        if not fused_tmp.GetBondBetweenAtoms(idx1_b, nb):
            em.AddBond(idx1_b, nb, bt)

    mol_fused = em.GetMol()
    Chem.SanitizeMol(mol_fused)
    return mol_fused

df_a = pd.read_csv(file_a)
df_b = pd.read_csv(file_b)
df_c = pd.read_csv(file_c)
df_d = pd.read_csv(file_d)

fused_ab_data = []
fused_cd_data = []
color_tag = {9: 109, 10: 110}

for i, row_a in df_a.iterrows():
    for j, row_b in df_b.iterrows():
        try:
            smi_a = row_a['Ring_SMILES']
            smi_b = row_b['Ring_SMILES']
            idx_a = eval(row_a['Ring_Atom_Indices'])
            idx_b = eval(row_b['Ring_Atom_Indices'])
            num_a = eval(row_a['Ring_Assigned_Numbers'])
            num_b = eval(row_b['Ring_Assigned_Numbers'])
            map_a = dict(zip(num_a, idx_a))
            map_b = dict(zip(num_b, idx_b))

            mol_a = Chem.MolFromSmiles(smi_a)
            mol_b = color_atoms(smi_b, map_b, color_tag)

            mol_ab = fuse_by_common_edge(mol_a, mol_b, map_a, map_b, 3, 2)
            mol_ab = mark_colored_atoms(mol_ab, color_tag)

            fused_smiles = Chem.MolToSmiles(mol_ab, isomericSmiles=True)
            b_colored = Chem.MolToSmiles(mol_b)

            fused_ab_data.append({
                "A_Label": row_a.get("Label", f"A{i}"),
                "B_Label": row_b.get("Label", f"B{j}"),
                "A_SMILES": smi_a,
                "B_SMILES": b_colored,
                "Fused_AB_SMILES": fused_smiles
            })

            img = Draw.MolToImage(mol_ab, size=(300, 300))
            img.save(os.path.join(output_img_folder, f"Fused_AB_{i}_{j}.png"))

        except:
            continue

for i, row_c in df_c.iterrows():
    for j, row_d in df_d.iterrows():
        try:
            smi_c = row_c['Ring_SMILES']
            smi_d = row_d['Ring_SMILES']
            idx_c = eval(row_c['Ring_Atom_Indices'])
            idx_d = eval(row_d['Ring_Atom_Indices'])
            num_c = eval(row_c['Ring_Assigned_Numbers'])
            num_d = eval(row_d['Ring_Assigned_Numbers'])
            map_c = dict(zip(num_c, idx_c))
            map_d = dict(zip(num_d, idx_d))

            mol_c = color_atoms(smi_c, map_c, color_tag)
            mol_d = Chem.MolFromSmiles(smi_d)

            mol_cd = fuse_by_common_edge(mol_c, mol_d, map_c, map_d, 13, 14)
            mol_cd = mark_colored_atoms(mol_cd, color_tag)

            fused_smiles = Chem.MolToSmiles(mol_cd, isomericSmiles=True)
            c_colored = Chem.MolToSmiles(mol_c)

            fused_cd_data.append({
                "C_Label": row_c.get("Label", f"C{i}"),
                "D_Label": row_d.get("Label", f"D{j}"),
                "C_SMILES": c_colored,
                "D_SMILES": smi_d,
                "Fused_CD_SMILES": fused_smiles
            })

            img = Draw.MolToImage(mol_cd, size=(300, 300))
            img.save(os.path.join(output_img_folder, f"Fused_CD_{i}_{j}.png"))

        except:
            continue

pd.DataFrame(fused_ab_data).to_csv(output_ab_csv, index=False)
pd.DataFrame(fused_cd_data).to_csv(output_cd_csv, index=False)
