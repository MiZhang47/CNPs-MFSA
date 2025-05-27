import os
import pandas as pd
from rdkit import Chem

mol_folder = 'C:/Users/xxx/Desktop/D_SMILES_x215/'
output_file = 'C:/Users/xxx/Desktop/D_SMILES_x215/D_SMILES_x215_Output.csv'

if not os.path.exists(mol_folder):
    raise FileNotFoundError(f"Folder {mol_folder} does not exist!")

smiles_list = []

for file_name in os.listdir(mol_folder):
    if file_name.endswith('.mol'):
        mol_path = os.path.join(mol_folder, file_name)

        try:
            with open(mol_path, 'r', encoding='utf-8') as f:
                mol_block = f.read()

            mol = Chem.MolFromMolBlock(mol_block, sanitize=False)
            if mol is not None:
                Chem.SanitizeMol(mol)
                smiles = Chem.MolToSmiles(mol)
                smiles_list.append({'FileName': file_name, 'SMILES': smiles})
            else:
                print(f"File {file_name} could not be parsed into a molecule (mol is None).")
        except Exception as e:
            print(f"File {file_name} failed to read. Skipped. Reason: {e}")

if smiles_list:
    df = pd.DataFrame(smiles_list)
    df.to_csv(output_file, index=False)
    print(f"SMILES saved to {output_file}")
else:
    print("No valid .mol files found or SMILES generation failed.")
