import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

base_path = r"C:\Users\xxx\Desktop\D_SMILES_x215\Unique_Rings"
file_ab = os.path.join(base_path, "Fused_AB_SMILES.csv")
file_c = os.path.join(base_path, "C_Ring_Unique_Renumbered.csv")
output_csv = os.path.join(base_path, "Fused_ABC_SMILES.csv")
output_img_dir = os.path.join(base_path, "Fused_ABC_Images")
os.makedirs(output_img_dir, exist_ok=True)

def get_atommap_mapping(mol):
    return {atom.GetAtomMapNum(): atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum() in [7, 8]}

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
df_c = pd.read_csv(file_c)
fused_data = []
count = 0

for i, row_ab in df_ab.iterrows():
    for j, row_c in df_c.iterrows():
        try:
            smi_ab = row_ab["Fused_AB_SMILES"]
            smi_c = row_c["Ring_SMILES"]

            mol_ab = Chem.MolFromSmiles(smi_ab)
            mol_c = Chem.MolFromSmiles(smi_c)

            map_ab = get_atommap_mapping(mol_ab)
            idx_c = eval(row_c["Ring_Atom_Indices"])
            num_c = eval(row_c["Ring_Assigned_Numbers"])
            map_c = dict(zip(num_c, idx_c))

            if 7 not in map_ab or 8 not in map_ab or 7 not in map_c or 8 not in map_c:
                print(f"⚠️ Skipped fusion {i}_{j}: missing atom 7 or 8")
                continue

            mol_fused = fuse_by_common_edge(mol_ab, mol_c, map_ab, map_c, 7, 8)
            fused_smi = Chem.MolToSmiles(mol_fused, isomericSmiles=True)

            img = Draw.MolToImage(mol_fused, size=(400, 400))
            img.save(os.path.join(output_img_dir, f"Fused_ABC_{i}_{j}.png"))

            fused_data.append({
                "A_Label": row_ab.get("A_Label", f"A_{i}"),
                "B_Label": row_ab.get("B_Label", f"B_{i}"),
                "C_Label": row_c.get("Label", f"C_{j}"),
                "A_SMILES": row_ab.get("A_SMILES", ""),
                "B_SMILES": row_ab.get("B_SMILES", ""),
                "C_SMILES": smi_c,
                "Fused_AB_SMILES": smi_ab,
                "Fused_ABC_SMILES": fused_smi,
                "A_Ring_Atom_Indices": row_ab.get("A_Ring_Atom_Indices", ""),
                "A_Ring_Assigned_Numbers": row_ab.get("A_Ring_Assigned_Numbers", ""),
                "B_Ring_Atom_Indices": row_ab.get("B_Ring_Atom_Indices", ""),
                "B_Ring_Assigned_Numbers": row_ab.get("B_Ring_Assigned_Numbers", ""),
                "C_Ring_Atom_Indices": row_c.get("Ring_Atom_Indices", ""),
                "C_Ring_Assigned_Numbers": row_c.get("Ring_Assigned_Numbers", ""),
                "A_FileName": row_ab.get("A_FileName", f"A_{i}"),
                "B_FileName": row_ab.get("B_FileName", f"B_{i}"),
                "C_FileName": row_c.get("FileName", f"C_{j}")
            })
            count += 1

        except Exception as e:
            print(f"❌ Fusion failed {i}_{j}: {e}")
            print(f"   ↳ AB_SMILES: {smi_ab}")
            print(f"   ↳ C_SMILES : {smi_c}")

df_out = pd.DataFrame(fused_data)
df_out.to_csv(output_csv, index=False, encoding="utf-8-sig")
print(f"\n✅ Total successful fusions: {count}, results saved to: {output_csv}")
