import os
import ast
import pandas as pd
from rdkit import Chem
from rdkit.Chem import RWMol
from rdkit.Chem import SmilesWriteParams
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

input_csv = r"C:\Users\xxx\Desktop\Fused_ABC_SMILES_With_Ring_Info.csv"

output_stage1_csv = r"C:\Users\xxx\Desktop\Stage1_Bridged.csv"
output_stage2_csv = r"C:\Users\xxx\Desktop\Stage2_ExtraBond.csv"
output_stage3_nh_csv = r"C:\Users\xxx\Desktop\Stage3_NH.csv"
output_stage3_nc_csv = r"C:\Users\xxx\Desktop\Stage3_NC.csv"
output_stage3_ncc_csv = r"C:\Users\xxx\Desktop\Stage3_NCC.csv"

output_stage1_img = r"C:\Users\xxx\Desktop\Stage1_Bridge_Images"
output_stage2_img = r"C:\Users\xxx\Desktop\Stage2_ExtraBond_Images"
output_stage3_nh_img = r"C:\Users\xxx\Desktop\Stage3_NH_Images"
output_stage3_nc_img = r"C:\Users\xxx\Desktop\Stage3_NC_Images"
output_stage3_ncc_img = r"C:\Users\xxx\Desktop\Stage3_NCC_Images"
os.makedirs(output_stage1_img, exist_ok=True)
os.makedirs(output_stage2_img, exist_ok=True)
os.makedirs(output_stage3_nh_img, exist_ok=True)
os.makedirs(output_stage3_nc_img, exist_ok=True)
os.makedirs(output_stage3_ncc_img, exist_ok=True)

bridge_smiles = "CO[CH:1][CH:2]O"
bridge_mol = Chem.MolFromSmiles(bridge_smiles)

def get_mapped_atoms(mol):
    result = {}
    for atom in mol.GetAtoms():
        if atom.HasProp("molAtomMapNumber"):
            result[atom.GetAtomMapNum()] = atom.GetIdx()
    return result

def safe_eval(val):
    try:
        return ast.literal_eval(val) if isinstance(val, str) else val
    except:
        return []

params = SmilesWriteParams()
params.kekuleSmiles = False
params.canonical = False
params.isomericSmiles = True

df = pd.read_csv(input_csv)
stage1_smiles, stage2_smiles = [], []
stage3_nh_smiles, stage3_nc_smiles, stage3_ncc_smiles = [], [], []

for i, row in df.iterrows():
    try:
        base_smiles = row["Fused_ABC_SMILES"]
        mol = Chem.MolFromSmiles(base_smiles, sanitize=False)
        if mol is None:
            for lst in [stage1_smiles, stage2_smiles, stage3_nh_smiles, stage3_nc_smiles, stage3_ncc_smiles]:
                lst.append("MolError")
            continue

        five_rings = safe_eval(row.get("Five_Rings", "[]"))
        six_rings = safe_eval(row.get("Six_Rings", "[]"))
        seven_rings = safe_eval(row.get("Seven_Rings", "[]"))
        five_std = safe_eval(row.get("Five_Ring_Standard_Numbers", "[]"))
        six_std = safe_eval(row.get("Six_Ring_Standard_Numbers", "[]"))
        seven_std = safe_eval(row.get("Seven_Ring_Standard_Numbers", "[]"))

        if len(five_rings) != 1 or 7 not in five_std or len(seven_rings) != 1 or 2 not in seven_std:
            for lst in [stage1_smiles, stage2_smiles, stage3_nh_smiles, stage3_nc_smiles, stage3_ncc_smiles]:
                lst.append("NoBridgeTarget")
            continue

        atom1 = five_rings[0][five_std.index(7)]
        atom2 = seven_rings[0][seven_std.index(2)]

        rw_mol = RWMol(mol)
        bridge_rw = RWMol(bridge_mol)
        bridge_map = get_mapped_atoms(bridge_rw)
        for atom in bridge_rw.GetAtoms():
            atom.ClearProp("molAtomMapNumber")

        combined = Chem.CombineMols(rw_mol, bridge_rw)
        combined_rw = RWMol(combined)
        offset = rw_mol.GetNumAtoms()
        combined_rw.AddBond(atom1, bridge_map[1] + offset, Chem.rdchem.BondType.SINGLE)
        combined_rw.AddBond(atom2, bridge_map[2] + offset, Chem.rdchem.BondType.SINGLE)

        stage1 = Chem.MolToSmiles(combined_rw, params)
        stage1_smiles.append(stage1)

        Chem.rdDepictor.Compute2DCoords(combined_rw)
        drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)
        drawer.drawOptions().addAtomIndices = True
        drawer.DrawMolecule(combined_rw, highlightAtoms=[bridge_map[1]+offset, bridge_map[2]+offset])
        drawer.FinishDrawing()
        with open(os.path.join(output_stage1_img, f"stage1_bridge_{i+1}.png"), "wb") as f:
            f.write(drawer.GetDrawingText())

        extra1 = seven_rings[0][seven_std.index(7)] if 7 in seven_std else None
        extra2 = seven_rings[0][seven_std.index(11)] + 1 if 11 in seven_std else None
        if extra1 is None or extra2 is None:
            for lst in [stage2_smiles, stage3_nh_smiles, stage3_nc_smiles, stage3_ncc_smiles]:
                lst.append("NoExtraBond")
            continue

        if not combined_rw.GetBondBetweenAtoms(extra1, extra2):
            combined_rw.AddBond(extra1, extra2, Chem.rdchem.BondType.SINGLE)

        stage2 = Chem.MolToSmiles(combined_rw, params)
        stage2_smiles.append(stage2)

        Chem.rdDepictor.Compute2DCoords(combined_rw)
        drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)
        drawer.drawOptions().addAtomIndices = True
        drawer.DrawMolecule(combined_rw, highlightAtoms=[extra1, extra2])
        drawer.FinishDrawing()
        with open(os.path.join(output_stage2_img, f"stage2_extrabond_{i+1}.png"), "wb") as f:
            f.write(drawer.GetDrawingText())

        pos1 = six_rings[0][six_std.index(1)] + 1 if 1 in six_std else None
        pos2 = seven_rings[0][seven_std.index(11)] + 1 if 11 in seven_std else None
        if pos1 is None or pos2 is None:
            for lst in [stage3_nh_smiles, stage3_nc_smiles, stage3_ncc_smiles]:
                lst.append("NoNLinkTarget")
            continue

        def add_n_group(base_rw, group_type, img_folder, filename):
            mol = RWMol(base_rw)
            n = mol.AddAtom(Chem.Atom(7))
            mol.GetAtomWithIdx(n).SetNumExplicitHs(0)
            highlight = [pos1, pos2, n]

            if group_type == "NH":
                h = mol.AddAtom(Chem.Atom(1))
                mol.AddBond(n, h, Chem.rdchem.BondType.SINGLE)
                highlight.append(h)
            elif group_type == "NC":
                c = mol.AddAtom(Chem.Atom(6))
                mol.AddBond(n, c, Chem.rdchem.BondType.SINGLE)
                highlight.append(c)
            elif group_type == "NCC":
                c1 = mol.AddAtom(Chem.Atom(6))
                c2 = mol.AddAtom(Chem.Atom(6))
                mol.AddBond(n, c1, Chem.rdchem.BondType.SINGLE)
                mol.AddBond(c1, c2, Chem.rdchem.BondType.SINGLE)
                highlight.extend([c1, c2])

            mol.AddBond(n, pos1, Chem.rdchem.BondType.SINGLE)
            mol.AddBond(n, pos2, Chem.rdchem.BondType.SINGLE)
            smiles = Chem.MolToSmiles(mol, params)

            Chem.rdDepictor.Compute2DCoords(mol)
            drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)
            drawer.drawOptions().addAtomIndices = True
            drawer.DrawMolecule(mol, highlightAtoms=highlight)
            drawer.FinishDrawing()
            with open(os.path.join(img_folder, filename), "wb") as f:
                f.write(drawer.GetDrawingText())
            return smiles

        stage3_nh = add_n_group(combined_rw, "NH", output_stage3_nh_img, f"stage3_nh_{i+1}.png")
        stage3_nc = add_n_group(combined_rw, "NC", output_stage3_nc_img, f"stage3_nc_{i+1}.png")
        stage3_ncc = add_n_group(combined_rw, "NCC", output_stage3_ncc_img, f"stage3_ncc_{i+1}.png")

        stage3_nh_smiles.append(stage3_nh)
        stage3_nc_smiles.append(stage3_nc)
        stage3_ncc_smiles.append(stage3_ncc)

    except Exception as e:
        err = f"Error: {str(e)}"
        for lst in [stage1_smiles, stage2_smiles, stage3_nh_smiles, stage3_nc_smiles, stage3_ncc_smiles]:
            lst.append(err)

df["Stage1_Bridged_SMILES"] = stage1_smiles
df.to_csv(output_stage1_csv, index=False)

df["Stage2_ExtraBond_SMILES"] = stage2_smiles
df.to_csv(output_stage2_csv, index=False)

df["Stage3_NH_SMILES"] = stage3_nh_smiles
df.to_csv(output_stage3_nh_csv, index=False)

df["Stage3_NC_SMILES"] = stage3_nc_smiles
df.to_csv(output_stage3_nc_csv, index=False)

df["Stage3_NCC_SMILES"] = stage3_ncc_smiles
df.to_csv(output_stage3_ncc_csv, index=False)
