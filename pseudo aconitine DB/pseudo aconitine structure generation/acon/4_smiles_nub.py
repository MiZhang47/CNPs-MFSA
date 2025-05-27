import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdmolops, Draw

input_csv = r"C:\Users\xxx\Desktop\D_SMILES_x215\Unique_Rings\Fused_ABC_SMILES.csv"
output_csv = r"C:\Users\xxx\Desktop\Fused_ABC_SMILES_With_Ring_Info.csv"
output_img_dir = r"C:\Users\xxx\Desktop\Fused_ABC_Ring_Images"
os.makedirs(output_img_dir, exist_ok=True)

standard_numbers_dict = {
    "A": [1, 2, 3, 4, 5, 6],
    "B": [8, 10, 11, 9, 1, 2, 7],
    "C": [7, 12, 13, 14, 8]
}

def filter_six_membered_rings_with_shared_atoms(six_rings, seven_rings):
    filtered = []
    for six in six_rings:
        six_set = set(six)
        for seven in seven_rings:
            shared_atoms = six_set & set(seven)
            if len(shared_atoms) >= 2:
                filtered.append(six)
                break
    return filtered

def get_rings_by_size(smi):
    try:
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if mol is None:
            return "Parse_Error", "Parse_Error", "Parse_Error", None

        five = []
        six = []
        seven = []

        for ring in rdmolops.GetSymmSSSR(mol):
            ring = list(ring)
            if len(ring) == 5:
                five.append(ring)
            elif len(ring) == 6:
                six.append(ring)
            elif len(ring) == 7:
                seven.append(ring)

        filtered_six = filter_six_membered_rings_with_shared_atoms(six, seven)

        return five, filtered_six, seven, mol

    except Exception as e:
        return f"Error: {str(e)}", "Error", "Error", None

df = pd.read_csv(input_csv)
smiles_list = df["Fused_ABC_SMILES"].tolist()

compound_ids = []
five_rings = []
six_rings = []
seven_rings = []
five_numbering = []
six_numbering = []
seven_numbering = []

for i, smi in enumerate(smiles_list):
    compound_id = f"ABC{i+1}"
    compound_ids.append(compound_id)

    five, six, seven, mol = get_rings_by_size(smi)
    five_rings.append(five)
    six_rings.append(six)
    seven_rings.append(seven)

    five_numbering.append(standard_numbers_dict["C"] if isinstance(five, list) and len(five) == 1 else [])
    six_numbering.append(standard_numbers_dict["A"] if isinstance(six, list) and len(six) == 1 else [])
    seven_numbering.append(standard_numbers_dict["B"] if isinstance(seven, list) and len(seven) == 1 else [])

    try:
        if isinstance(mol, Chem.Mol):
            for atom in mol.GetAtoms():
                atom.SetProp("atomLabel", str(atom.GetIdx()))

            drawer = Draw.MolDraw2DCairo(500, 500)
            drawer.drawOptions().addAtomIndices = True
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()

            img_path = os.path.join(output_img_dir, f"{compound_id}.png")
            with open(img_path, "wb") as f:
                f.write(drawer.GetDrawingText())
            print(f"🖼️ Image generated: {compound_id}")
    except Exception as e:
        print(f"❌ Image generation failed: {compound_id}, error: {e}")

df["Compound_ID"] = compound_ids
df["Five_Rings"] = five_rings
df["Six_Rings"] = six_rings
df["Seven_Rings"] = seven_rings
df["Five_Ring_Standard_Numbers"] = five_numbering
df["Six_Ring_Standard_Numbers"] = six_numbering
df["Seven_Ring_Standard_Numbers"] = seven_numbering

df.to_csv(output_csv, index=False)
print(f"\n✅ CSV saved: {output_csv}")
print(f"✅ Images saved to: {output_img_dir}")
