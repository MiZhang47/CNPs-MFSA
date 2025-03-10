"""
Last updated: 2024/12
"""

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

path_csv = os.path.join('C:/Users/zhang/Desktop/StrucSimulation6/', 'MD_67_epo_diol_ene(1296).csv')
df = pd.read_csv(path_csv)

# Define folder for save image
img_folder = 'C:/Users/zhang/Desktop/StrucSimulation6/MD_67_epo_diol_ene(1296)/'
if not os.path.exists(img_folder):
    os.makedirs(img_folder)

# Read rows 'Number' and 'New SMILES'
for index, row in df.iterrows():
    number = row['Number']
    smiles = row['New SMILES']

    # Create mol file from 'New SMILES'
    mol = Chem.MolFromSmiles(smiles)

    # Draw structures
    img = Draw.MolToImage(mol)

    # Save image
    img_path = os.path.join(img_folder, f'Structure_{number}.png')
    Draw.MolToFile(mol, img_path, size=(300, 300), imageType='png', dpi=300)

    # Add new row structures
    df.loc[index, 'Structure'] = img_path

# save csv file with structure
df.to_csv(os.path.join('C:/Users/zhang/Desktop/StrucSimulation6/', 'MD_67_epo_diol_ene_struc(1296).csv'), index=False)