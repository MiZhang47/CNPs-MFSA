"""
Last updated: 2024/12
MF calculation in positive & negative ion modes (M, M+NH3, M+HCOOH)
"""

import re
import os
import csv
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# Redefine the new fragments
frag1 = 'C(=O)O'  # carbonyl group
frag2 = 'C=C'  # double bond
frag3 = 'C(O)'  # hydroxyl group
frag4 = 'C(OC)'  # methoxy group
frag5 = 'C1OC1'  # epoxy group
frag6 = 'C'  # carbon

# Function to generate the SMILES strings according to the new rules
def generate_smiles_new():
    results = []
    no = 1

    # frag6 + frag1, f=1~18-a, a=1
    for f in range(1, 18):
        smiles = frag6 * f + frag1
        results.append((f"R{no}", smiles))
        no += 1

    # frag6 + frag2 + frag1, f=1~18-a-2b, a=1, b=1~3
    for b in range(1, 4):
        for f in range(1, 18 - 2 * b):
            smiles = frag6 * f + frag2 * b + frag1
            results.append((f"R{no}", smiles))
            no += 1

    # frag6 + frag5 + frag2 + frag1, f=1~10-a-2b-2e, a=1, b=2~3, e=1
    for b in range(2, 4):
        for f in range(1, 9 - 2 * b):
            smiles = frag6 * f + frag5 + frag2 * b + frag1
            results.append((f"R{no}", smiles))
            no += 1

    # frag6 + frag4 + frag2 + frag1, f=1~10-a-2b-2d, a=1, b=2~3, d=1~2
    for d in range(1, 3):
        for b in range(2, 4):
            for f in range(1, 10 - 2 * b - 2 * d):
                smiles = frag6 * f + frag4 * d + frag2 * b + frag1
                results.append((f"R{no}", smiles))
                no += 1

    # frag6 + frag3 + frag2 + frag1, f=1~10-a-2b-2c, a=1, b=2~3, c=1~2
    for c in range(1, 3):
        for b in range(2, 4):
            for f in range(1, 10 - 2 * b - 2 * c):
                smiles = frag6 * f + frag3 * c + frag2 * b + frag1
                results.append((f"R{no}", smiles))
                no += 1

    # frag6 + frag3 + frag4 + frag2 + frag1, f=1~10-a-2b-c-d, a=1, b=2~3, d=1, c=1
    for c in range(1, 2):
        for b in range(2, 4):
            for f in range(1, 10 - 2 * b - c):
                smiles = frag6 * f + frag3 * c + frag4 + frag2 * b + frag1
                results.append((f"R{no}", smiles))
                no += 1

    return results

# Generation Molecular formula in positive ion mode
def calculate_molecular_formula(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        formula = rdMolDescriptors.CalcMolFormula(mol)
        # Subtract H2O from the molecular formula
        elements = {}
        # Parse the molecular formula to count elements
        for element in re.findall(r'([A-Z][a-z]?)(\d*)', formula):
            el, count = element
            count = int(count) if count else 1
            elements[el] = elements.get(el, 0) + count
        # Subtract H2 and O
        elements['H'] = elements.get('H', 0) - 2
        elements['O'] = elements.get('O', 0) - 1
        # Reconstruct the formula
        new_formula = ''
        for el, count in elements.items():
            if count > 0:  # Ensure no negative counts
                new_formula += f"{el}{count if count > 1 else ''}"
        return new_formula
    else:
        return "Invalid SMILES"

# Generate the new SMILES strings
smiles_results_new = generate_smiles_new()

# Define the output path for the CSV file
output_path = os.path.join('C:/Users/zhang/Desktop/', 'Acyl.csv')

# Write to CSV with the specified format
with open(output_path, mode='w', newline='', encoding='utf-8') as file:
    writer = csv.writer(file)
    writer.writerow(['no.', 'SMILES', 'Molecular Formula'])
    for r, smiles in smiles_results_new:
        mf = calculate_molecular_formula(smiles)
        writer.writerow([r, smiles, mf])