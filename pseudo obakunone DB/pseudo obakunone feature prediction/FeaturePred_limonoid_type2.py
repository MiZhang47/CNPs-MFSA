"""
Last updated: 2025/05
Structure-based MS feature prediction for pDB-limonoid type II
"""

import os
import pandas as pd

# Define shift value of A parts
A_shift_set = {
    'A1': {'C': 0, 'H': 0, 'O': 0},  # 1,2-en-3-one
    'A2': {'C': 0, 'H': 2, 'O': 0},  # 1,2-en-3-ol
}

# Define shift value of B parts
B_shift_set = {
    'B1': {'C': 2, 'H': 2, 'O': 1},  # 7-OAc
    'B2': {'C': 0, 'H': 0, 'O': 0},  # 7-ol
    'B3': {'C': 0, 'H': -2, 'O': 0}  # 7-one
}

# Define shift value of C parts
C_shift_set = {
    'C1': {'C': 0, 'H': 0, 'O': 0},  # no sub
}

# Define shift value of D parts
D_shift_set = {
    'D6': {'C': 0, 'H': 2, 'O': 1},  # 14,15-epoxy, 16-one, furan
}

ref_feature_TypeII = [
    {'C': 5,  'H': 4,  'O': 1},  # charge = +1, Theo.mass 81.034
    {'C': 8,  'H': 10, 'O': 0},  # charge = +1, Theo.mass 107.0858
    {'C': 6,  'H': 4,  'O': 2},  # charge = +1, Theo.mass 109.0287
    {'C': 9,  'H': 12, 'O': 0},  # charge = +1, Theo.mass 121.1013
    {'C': 9,  'H': 10, 'O': 1},  # charge = +1, Theo.mass 135.0805
    {'C': 9,  'H': 12, 'O': 1},  # charge = +1, Theo.mass 137.0961
    {'C': 10, 'H': 12, 'O': 1},  # charge = +1, Theo.mass 149.0961
    {'C': 11, 'H': 12, 'O': 1},  # charge = +1, Theo.mass 161.096
    {'C': 11, 'H': 14, 'O': 1},  # charge = +1, Theo.mass 163.1116
    {'C': 12, 'H': 14, 'O': 1},  # charge = +1, Theo.mass 175.1116
    {'C': 13, 'H': 14, 'O': 1},  # charge = +1, Theo.mass 187.1116
    {'C': 14, 'H': 14, 'O': 1},  # charge = +1, Theo.mass 199.1116
    {'C': 15, 'H': 16, 'O': 1},  # charge = +1, Theo.mass 213.1272
    {'C': 17, 'H': 20, 'O': 1},  # charge = +1, Theo.mass 241.1584
    {'C': 16, 'H': 18, 'O': 2},  # charge = +1, Theo.mass 243.1376
    {'C': 16, 'H': 16, 'O': 3},  # charge = +1, Theo.mass 257.1169
    {'C': 16, 'H': 18, 'O': 3},  # charge = +1, Theo.mass 259.1326
    {'C': 19, 'H': 24, 'O': 1},  # charge = +1, Theo.mass 269.1895
    {'C': 17, 'H': 18, 'O': 3},  # charge = +1, Theo.mass 271.1325
    {'C': 20, 'H': 26, 'O': 2},  # charge = +1, Theo.mass 299.2
    {'C': 26, 'H': 28, 'O': 3},  # charge = +1, Theo.mass 389.2105
    {'C': 26, 'H': 30, 'O': 4},  # charge = +1, Theo.mass 407.221
]


path_csv = os.path.join('C:/Users/zhang/Desktop/limonoid/', 'lim_epo_ABCD_SMILES.csv')
df = pd.read_csv(path_csv)  # Read the CSV file using the pandas and store the data in a DataFrame object named df.

for index, row in df.iterrows():  # Iterate over each row in the DataFrame.
    new_feature_set = []  # Create an empty list to store the calculated new feature molecular formulas.
    for feature in ref_feature_TypeII:  # Traverse the reference feature set to calculate the new molecular formula.
        # Create a dictionary to store the calculated new feature molecular formulas.
        modified_feature = {
            'C': feature['C'] + A_shift_set[row['A parts']]['C'] + B_shift_set[row['B parts']]['C']
                + C_shift_set[row['C parts']]['C'] + D_shift_set[row['D parts']]['C'],
            'H': feature['H'] + A_shift_set[row['A parts']]['H'] + B_shift_set[row['B parts']]['H']
                + C_shift_set[row['C parts']]['H'] + D_shift_set[row['D parts']]['H'],
            'O': feature['O'] + A_shift_set[row['A parts']]['O'] + B_shift_set[row['B parts']]['O']
                + C_shift_set[row['C parts']]['O'] + D_shift_set[row['D parts']]['O'],
        }
        new_feature_set.append(modified_feature)  # Add the calculated new feature formula to the new_feature_set list.

    for i, feature in enumerate(new_feature_set, 1):
        # Iterate over each element in the new_feature_set list.
        # Return an enumeration object that yields a pair of values:
        # the index (starting from 1) and the corresponding element (dictionary of feature molecular formulas)
        if feature['O'] == 0:  # Check whether the number of O atoms in the current feature molecular formula is 0.
            df.loc[index, f'ion formula {i}'] = f"C{feature['C']}H{feature['H']}"
            # Store the new feature formula (excluding the O atom) as a string in the current row (index)
            # of the DataFrame with column name "ion formula i".
        else:
            df.loc[index, f'ion formula {i}'] = f"C{feature['C']}H{feature['H']}O{feature['O']}"
            # Store the new feature formula (including the O atom) as a string in the current row (index)
            # of the DataFrame with column name "ion formula i".

output_path = os.path.join('C:/Users/zhang/Desktop/limonoid/', 'lim_epo_ABCD_SMILES_feature'
                                                                                   'ion.csv')
df.to_csv(output_path, index=False)