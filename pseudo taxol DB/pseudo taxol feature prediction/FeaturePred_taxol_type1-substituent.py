"""
Last updated: 2025/05
Structure-based MS feature prediction for pDB-limonoid type I
"""

import os
import pandas as pd

# Define shift value of A parts
A_shift_set = {
    'A1': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'A2': {'C': 2, 'H': 2, 'O': 1, 'N': 0},
    'A3': {'C': 2, 'H': 2, 'O': 0, 'N': 0},
    'A4': {'C': 1, 'H': 2, 'O': 0, 'N': 0},
    'A5': {'C': -1, 'H': 2, 'O': 0, 'N': 0},
}

# Define shift value of B parts
B_shift_set = {
    'B1': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'B2': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'B3': {'C': 0, 'H': 0, 'O': 0, 'N': 0},  # initial
    'B4': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'B5': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'B6': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'B7': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
}

# Define shift value of C parts
C_shift_set = {
    'C1': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'C2': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
}

ref_feature_TypeI = [
    {'C': 16, 'H': 15, 'O': 4, 'N': 1},  # charge = +1, Theo.mass 286
]

path_csv = os.path.join('C:/Users/zhang/Desktop/taxane/', 'tax_ABC_SMILES.csv')
df = pd.read_csv(path_csv)  # Read the CSV file using the pandas and store the data in a DataFrame object named df.

for index, row in df.iterrows():  # Iterate over each row in the DataFrame.
    new_feature_set = []  # Create an empty list to store the calculated new feature molecular formulas.
    for feature in ref_feature_TypeI:  # Traverse the reference feature set to calculate the new molecular formula.
        # Create a dictionary to store the calculated new feature molecular formulas.
        modified_feature = {
            'C': feature['C'] + A_shift_set[row['A parts']]['C'] + B_shift_set[row['B parts']]['C']
                + C_shift_set[row['C parts']]['C'],
            'H': feature['H'] + A_shift_set[row['A parts']]['H'] + B_shift_set[row['B parts']]['H']
                + C_shift_set[row['C parts']]['H'],
            'O': feature['O'] + A_shift_set[row['A parts']]['O'] + B_shift_set[row['B parts']]['O']
                + C_shift_set[row['C parts']]['O'],
            'N': feature['N'] + A_shift_set[row['A parts']]['N'] + B_shift_set[row['B parts']]['N']
                 + C_shift_set[row['C parts']]['N']
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
            df.loc[index, f'ion formula {i}'] = f"C{feature['C']}H{feature['H']}O{feature['O']}N{feature['N']}"
            # Store the new feature formula (including the O atom) as a string in the current row (index)
            # of the DataFrame with column name "ion formula i".

output_path = os.path.join('C:/Users/zhang/Desktop/taxane/', 'tax_ABC_SMILES_feature'
                                                                                   'ion-2.csv')
df.to_csv(output_path, index=False)