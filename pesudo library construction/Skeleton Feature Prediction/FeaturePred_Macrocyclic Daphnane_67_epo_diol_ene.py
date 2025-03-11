"""
Last updated: 2024/04
Structure-based MS feature prediction for normal macrocyclic daphnane
"""

import os
import pandas as pd

# Define shift value of A parts
A_shift_set = {
    'A8': {'C': 0, 'H': 0, 'O': 0},  # 1-alkyl-3-one, reference A ring structure
    'A9': {'C': 0, 'H': 0, 'O': 1},  # 1-alkyl-2-ol-3-one
    'A10': {'C': 0, 'H': 2, 'O': 0},  # 1-alkyl-3-ol
    'A11': {'C': 0, 'H': 2, 'O': 1},  # 1-alkyl-2-ol-3-ol
    'A12': {'C': 0, 'H': -2, 'O': 1},  # 3,4-seco ring
    'A13': {'C': 0, 'H': 0, 'O': 1},  # bicyclo[2.2.1]heptane ring
}

# Define shift value of B parts
B_shift_set = {
    'B1': {'C': 0, 'H': 0, 'O': 0},  # 5-ol-6,7-epoxy, reference B ring structure
    'B2': {'C': 0, 'H': 0, 'O': -1},  # 5-dehydro-6,7-epoxy
    'B3': {'C': 0, 'H': 2, 'O': 1},  # 5-ol-6,7-diol
    'B4': {'C': 0, 'H': 2, 'O': 0},  # 5-dehydro-6,7-diol
    'B5': {'C': 0, 'H': 0, 'O': -1},  # 5-ol-6,7-ene
    'B6': {'C': 0, 'H': 2, 'O': -2},  # 5-dehydro-6,7-ene
}

# Define shift value of C parts
C_shift_set = {
    'C1': {'C': 0, 'H': 0, 'O': 0},  # 9,13,14-orthoester, reference C ring structure
    'C2': {'C': 0, 'H': 0, 'O': 1},  # 9,13,14-orthoester, 12-ol
    'C5': {'C': 0, 'H': 0, 'O': 1},  # 9,13,14-orthoester, 18-ol
    'C6': {'C': 0, 'H': 0, 'O': 2},  # 9,13,14-orthoester, 12-ol-18-ol
}

# Define shift value of M parts
M_shift_set = {
    # C9 macrocyclic ring (M9)
    'M1': {'C': -1, 'H': -2, 'O': 0},  # without acylation
    # C10 macrocyclic ring (M10)
    'M2': {'C': 0, 'H': 0, 'O': 0},  # without acylation, reference C ring structure
    'M3': {'C': 0, 'H': 0, 'O': 1},  # 2'-ol
    'M4': {'C': 0, 'H': 0, 'O': 2},  # 2'-ol-7'-ol
    'M5': {'C': 0, 'H': 0, 'O': 3},  # 2'-ol-6'-ol-7'-ol
    'M6': {'C': 0, 'H': 2, 'O': 4},  # 2'-ol-5'-ol-6'-ol-7'-ol
    # C14 macrocyclic ring (M14)
    'M7': {'C': 4, 'H': 6, 'O': 0},  # 5',6'-ene
    # C16 macrocyclic ring (M16)
    'M8': {'C': 6, 'H': 10, 'O': 0},  # 2',3'-ene
    'M9': {'C': 6, 'H': 12, 'O': 1},  # 15'-ol
}

# Define reference feature molecular formula for macrocyclic daphnane
ref_feature_setMD = [  # reference MS features from TDPE-5
    {'C': 30, 'H': 44, 'O': 8},  # charge = +1, Theo.mass 533.31089
    {'C': 30, 'H': 42, 'O': 7},  # charge = +1, Theo.mass 515.30033
    {'C': 30, 'H': 40, 'O': 6},  # charge = +1, Theo.mass 497.28977
    {'C': 30, 'H': 38, 'O': 5},  # charge = +1, Theo.mass 479.27920
    {'C': 30, 'H': 36, 'O': 4},  # charge = +1, Theo.mass 461.26864
    {'C': 29, 'H': 38, 'O': 4},  # charge = +1, Theo.mass 451.28429
    {'C': 29, 'H': 36, 'O': 4},  # charge = +1, Theo.mass 449.26864
    {'C': 30, 'H': 34, 'O': 3},  # charge = +1, Theo.mass 443.25807
    {'C': 29, 'H': 36, 'O': 3},  # charge = +1, Theo.mass 433.27372
    {'C': 29, 'H': 34, 'O': 2},  # charge = +1, Theo.mass 415.26316
    {'C': 28, 'H': 34, 'O': 2},  # charge = +1, Theo.mass 403.26316
]

path_csv = os.path.join('C:/Users/zhang/Desktop/StrucSimulation5/', 'MD_67_epo_diol_ene(1296).csv')
df = pd.read_csv(path_csv)  # Read the CSV file using the pandas and store the data in a DataFrame object named df.

for index, row in df.iterrows():  # Iterate over each row in the DataFrame.
    new_feature_set = []  # Create an empty list to store the calculated new feature molecular formulas.
    for feature in ref_feature_setMD:  # Traverse the reference feature set to calculate the new molecular formula.
        # Create a dictionary to store the calculated new feature molecular formulas.
        modified_feature = {
            'C': feature['C'] + A_shift_set[row['A parts']]['C'] + B_shift_set[row['B parts']]['C']
                + C_shift_set[row['C parts']]['C'] + M_shift_set[row['Macro parts']]['C'],
            'H': feature['H'] + A_shift_set[row['A parts']]['H'] + B_shift_set[row['B parts']]['H']
                + C_shift_set[row['C parts']]['H'] + M_shift_set[row['Macro parts']]['H'],
            'O': feature['O'] + A_shift_set[row['A parts']]['O'] + B_shift_set[row['B parts']]['O']
                + C_shift_set[row['C parts']]['O'] + M_shift_set[row['Macro parts']]['O']
        }
        new_feature_set.append(modified_feature)  # Add the calculated new feature formula to the new_feature_set list.

    for i, feature in enumerate(new_feature_set, 1):
        # Iterate over each element in the new_feature_set list.
        # Return an enumeration object that yields a pair of values:
        # the index (starting from 1) and the corresponding element (dictionary of feature molecular formulas)
        if feature['O'] == 0:  # Check whether the number of O atoms in the current feature molecular formula is 0.
            df.loc[index, f'MS Feature {i}'] = f"C{feature['C']}H{feature['H']}"
            # Store the new feature formula (excluding the O atom) as a string in the current row (index)
            # of the DataFrame with column name "MS Feature i".
        else:
            df.loc[index, f'MS Feature {i}'] = f"C{feature['C']}H{feature['H']}O{feature['O']}"
            # Store the new feature formula (including the O atom) as a string in the current row (index)
            # of the DataFrame with column name "MS Feature i".

output_path = os.path.join('C:/Users/zhang/Desktop/FeaturePredict5/', 'MD_67_epo_diol_ene(1296)_feature.csv')
df.to_csv(output_path, index=False)