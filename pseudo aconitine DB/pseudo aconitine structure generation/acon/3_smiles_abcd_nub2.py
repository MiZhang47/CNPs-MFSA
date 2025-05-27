import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

base_path = r"C:\Users\xxx\Desktop\D_SMILES_x215\Unique_Rings"
file_a = os.path.join(base_path, "A_Ring_Unique_Renumbered.csv")
file_b = os.path.join(base_path, "B_Ring_Unique_Renumbered.csv")
output_ab_csv = os.path.join(base_path, "Fused_AB_SMILES.csv")
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

def fuse_by_common_edge_two(mol_a, mol_b, map_a, map_b, shared_labels):
    idx_a = [map_a[l] for l in shared_labels]
    idx_b = [map_b[l] for l in shared_labels]
    offset = mol_a.GetNumAtoms()

    anchors, nbs, bts = [], [], []
    for l in shared_labels:
        src_b = map_b[l]
        tgt_a = map_a[l]
        atom_b = mol_b.GetAtomWithIdx(src_b)
        for nbr in atom_b.GetNeighbors():
            n_idx = nbr.GetIdx()
            if n_idx in idx_b:
                continue
            anchors.append(tgt_a)
            nbs.append(n_idx + offset)
            bts.append(mol_b.GetBondBetweenAtoms(src_b, n_idx).GetBondType())

    combo = Chem.CombineMols(mol_a, mol_b)
    em = Chem.EditableMol(combo)

    remove = sorted([idx + offset for idx in idx_b], reverse=True)
    for i in remove:
        em.RemoveAtom(i)

    def shift(i): return i - sum(1 for r in remove if i > r)
    nbs = [shift(i) for i in nbs]
    anchors = [shift(i) for i in anchors]

    fused_tmp = em.GetMol()
    for a, nb, bt in zip(anchors, nbs, bts):
        if not fused_tmp.GetBondBetweenAtoms(a, nb):
            em.AddBond(a, nb, bt)

    mol_fused = em.GetMol()
    Chem.SanitizeMol(mol_fused)
    return mol_fused

df_a = pd.read_csv(file_a)
df_b = pd.read_csv(file_b)
fused_ab_data = []
color_tag = {8: 8, 7: 7}
label_counter = 1

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

            mol_ab = fuse_by_common_edge_two(mol_a, mol_b, map_a, map_b, [1, 2])
            mol_ab = mark_colored_atoms(mol_ab, color_tag)

            fused_smiles = Chem.MolToSmiles(mol_ab, isomericSmiles=True)
            b_colored = Chem.MolToSmiles(mol_b)

            fused_ab_data.append({
                "Label": f"AB{label_counter}",
                "A_Label": row_a.get("Label", f"A_{i}"),
                "B_Label": row_b.get("Label", f"B_{j}"),
                "A_SMILES": smi_a,
                "B_SMILES": b_colored,
                "Fused_AB_SMILES": fused_smiles,
                "A_Ring_Atom_Indices": row_a['Ring_Atom_Indices'],
                "A_Ring_Assigned_Numbers": row_a['Ring_Assigned_Numbers'],
                "A_FileName": row_a.get('FileName', f"A_{i}"),
                "B_Ring_Atom_Indices": row_b['Ring_Atom_Indices'],
                "B_Ring_Assigned_Numbers": row_b['Ring_Assigned_Numbers'],
                "B_FileName": row_b.get('FileName', f"B_{j}"),
                "Compound_ID_A": row_a.get('Compound_ID', ''),
                "Compound_ID_B": row_b.get('Compound_ID', ''),
                "Ring_Type_A": row_a.get('Ring_Type', ''),
                "Ring_Type_B": row_b.get('Ring_Type', '')
            })

            img = Draw.MolToImage(mol_ab, size=(300, 300))
            img.save(os.path.join(output_img_folder, f"Fused_AB_{i}_{j}.png"))

            label_counter += 1

        except Exception as e:
            print(f"❌ AB fusion failed {i}_{j}: {e}")

pd.DataFrame(fused_ab_data).to_csv(output_ab_csv, index=False, encoding="utf-8-sig")
print("✅ AB fusion completed. Results saved to:", output_ab_csv)
