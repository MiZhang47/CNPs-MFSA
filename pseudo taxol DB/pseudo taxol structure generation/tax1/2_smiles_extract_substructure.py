import os
import pandas as pd
import traceback
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import SmilesWriteParams

input_csv_path = r"C:/Users/xxx/Desktop/D_SMILES_x215/D_SMILES_x215_Output.csv"
output_csv_path = r"C:/Users/xxx/Desktop/D_SMILES_x215/D_SMILES_x215_RingsExtracted_NoFusion.csv"
output_img_dir = r"C:/Users/xxx/Desktop/D_SMILES_x215/standardized_mol_images"
ring_unique_dir = r"C:/Users/xxx/Desktop/D_SMILES_x215/Unique_Rings"
os.makedirs(output_img_dir, exist_ok=True)
os.makedirs(ring_unique_dir, exist_ok=True)

df = pd.read_csv(input_csv_path)
scaffold_rows = df[df["FileName"].str.contains("aco1_0")]
target_rows = df[~df["FileName"].str.contains("aco1_0")]

scaffold_mols = []
for _, row in scaffold_rows.iterrows():
    raw_smiles = row["SMILES"]
    mol_with_dummy = Chem.MolFromSmiles(raw_smiles)
    if mol_with_dummy is None:
        continue
    editable = Chem.EditableMol(mol_with_dummy)
    to_remove = [atom.GetIdx() for atom in mol_with_dummy.GetAtoms() if atom.GetAtomicNum() == 0]
    for idx in sorted(to_remove, reverse=True):
        editable.RemoveAtom(idx)
    scaffold = editable.GetMol()
    Chem.SanitizeMol(scaffold)
    rdDepictor.Compute2DCoords(scaffold)
    scaffold_mols.append((row["FileName"], scaffold))

def extract_ring_submol(mol, ring_atoms, all_ring_atoms, owned_atoms):
    ring_atom_set = set(ring_atoms)
    all_rings_set = set(all_ring_atoms)
    selected_atoms = set(ring_atoms)
    ri = mol.GetRingInfo()
    all_rings = [set(r) for r in ri.AtomRings()]

    def expand_substituent(atom_idx, visited):
        if atom_idx in visited:
            return
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        for r in all_rings:
            if atom_idx in r and len(r) <= 6:
                for ridx in r:
                    if ridx not in visited and ridx not in ring_atom_set and ridx not in all_rings_set:
                        expand_substituent(ridx, visited)
                break
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx not in ring_atom_set and n_idx not in all_rings_set and n_idx not in visited:
                expand_substituent(n_idx, visited)

    visited = set()
    for idx in owned_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx not in ring_atom_set and n_idx not in all_rings_set:
                expand_substituent(n_idx, visited)

    selected_atoms.update(visited)

    aromatic_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    aromatic_rings = [set(r) for r in ri.AtomRings() if all(a in aromatic_atoms for a in r)]
    for aro_ring in aromatic_rings:
        if selected_atoms & aro_ring:
            selected_atoms.update(aro_ring)

    selected_atoms = sorted(selected_atoms)
    atom_map = {idx: i for i, idx in enumerate(selected_atoms)}
    bonds = [bond.GetIdx() for bond in mol.GetBonds()
             if bond.GetBeginAtomIdx() in selected_atoms and bond.GetEndAtomIdx() in selected_atoms]

    submol = Chem.PathToSubmol(mol, bonds, atomMap=atom_map)
    return submol, atom_map, selected_atoms

results = []
failed = []
ring_records = {'A': [], 'B': [], 'C': []}

for _, row in target_rows.iterrows():
    file_name = row['FileName']
    raw_smiles = row['SMILES']
    mol_with_dummy = Chem.MolFromSmiles(raw_smiles)

    if mol_with_dummy is None:
        print(f"❌ Invalid SMILES skipped: {file_name}")
        failed.append(file_name)
        continue

    editable = Chem.EditableMol(mol_with_dummy)
    to_remove = [atom.GetIdx() for atom in mol_with_dummy.GetAtoms() if atom.GetAtomicNum() == 0]
    for idx in sorted(to_remove, reverse=True):
        editable.RemoveAtom(idx)
    mol = editable.GetMol()
    Chem.SanitizeMol(mol)

    matched = False
    for _, scaffold in scaffold_mols:
        match = mol.GetSubstructMatch(scaffold)
        if match:
            matched = True
            core_atoms = list(match)
            non_core_atoms = [i for i in range(mol.GetNumAtoms()) if i not in core_atoms]
            new_order = core_atoms + non_core_atoms
            mol = Chem.RenumberAtoms(mol, new_order)
            break

    if not matched:
        print(f"⚠️ Scaffold not matched: {file_name}")
        failed.append(file_name)
        continue

    try:
        rdDepictor.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)
        opts = drawer.drawOptions()
        for i in range(mol.GetNumAtoms()):
            opts.atomLabels[i] = str(i)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        img_path = os.path.join(output_img_dir, file_name.replace(".mol", "_aligned.png"))
        with open(img_path, "wb") as f:
            f.write(drawer.GetDrawingText())

        rings = mol.GetRingInfo().AtomRings()
        six_membered = [set(r) for r in rings if len(r) == 6]
        seven_membered = [set(r) for r in rings if len(r) == 7]
        five_membered = [set(r) for r in rings if len(r) == 5]

        A_ring = B_ring = C_ring = None
        for a in six_membered:
            for b in seven_membered:
                if len(a & b) >= 2:
                    for c in five_membered:
                        if len(b & c) >= 2:
                            A_ring = list(a)
                            B_ring = list(b)
                            C_ring = list(c)
                            break
                    if C_ring:
                        break
            if C_ring:
                break

        if not (A_ring and B_ring and C_ring):
            raise ValueError("A-B-C rings not found with common edge")

        A_ring_set = set(A_ring)
        B_ring_set = set(B_ring)
        C_ring_set = set(C_ring)

        AB_shared = A_ring_set & B_ring_set
        BC_shared = B_ring_set & C_ring_set

        A_owned = A_ring_set - B_ring_set
        B_owned = B_ring_set - A_ring_set - C_ring_set
        C_owned = C_ring_set - B_ring_set

        B_owned |= AB_shared
        C_owned |= BC_shared

        all_ring_atoms = A_ring_set | B_ring_set | C_ring_set
        ring_defs = [
            ('A', A_ring, A_owned),
            ('B', B_ring, B_owned),
            ('C', C_ring, C_owned),
        ]

        submols_smiles = []
        atom_lists = []
        for label, ring, owned_atoms in ring_defs:
            submol, atom_map, selected_atoms = extract_ring_submol(mol, ring, all_ring_atoms, owned_atoms)
            params = SmilesWriteParams()
            params.canonical = False
            params.atomOrdering = list(range(submol.GetNumAtoms()))
            submol_smiles = Chem.MolToSmiles(submol, params)
            submols_smiles.append(submol_smiles)
            atom_str = "-".join(str(i) for i in selected_atoms)
            atom_lists.append(atom_str)
            ring_records[label].append({
                "FileName": file_name,
                "Ring_SMILES": submol_smiles,
                "Ring_Atoms": atom_str
            })

        results.append({
            "FileName": file_name,
            "Standardized_SMILES": Chem.MolToSmiles(mol, canonical=False),
            "Standard_Atom_Indexes": "-".join(str(i) for i in range(mol.GetNumAtoms())),
            "A_Ring_SMILES": submols_smiles[0],
            "A_Ring_Atoms": atom_lists[0],
            "B_Ring_SMILES": submols_smiles[1],
            "B_Ring_Atoms": atom_lists[1],
            "C_Ring_SMILES": submols_smiles[2],
            "C_Ring_Atoms": atom_lists[2],
        })

        print(f"{file_name} ✅ Done")

    except Exception as e:
        print(f"{file_name} ❌ Error: {e}")
        traceback.print_exc()
        failed.append(file_name)

if results:
    pd.DataFrame(results).to_csv(output_csv_path, index=False, encoding='utf-8-sig')
    print(f"\n✅ Saved to: {output_csv_path}")

for ring_label in ['A', 'B', 'C']:
    ring_df = pd.DataFrame(ring_records[ring_label])
    ring_unique = ring_df.drop_duplicates(subset=["Ring_SMILES"])
    unique_path = os.path.join(ring_unique_dir, f"{ring_label}_Ring_Unique.csv")
    ring_unique.to_csv(unique_path, index=False, encoding='utf-8-sig')
    print(f"✅ {ring_label} ring unique count: {len(ring_unique)}, saved to {unique_path}")

if failed:
    fail_log = output_csv_path.replace('.csv', '_failed.csv')
    pd.DataFrame({'FailedFile': failed}).to_csv(fail_log, index=False, encoding='utf-8-sig')
    print(f"⚠️ Failed structures saved to: {fail_log}")
