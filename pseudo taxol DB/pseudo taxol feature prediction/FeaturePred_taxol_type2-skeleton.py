"""
Last updated: 2025/05
Structure-based MS feature prediction for pDB-limonoid type I
"""

import os
import pandas as pd

# Define shift value of A parts
A_shift_set = {
    'A6': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'A7': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'A8': {'C': 0, 'H': -2, 'O': 0, 'N': 0},
    'A9': {'C': 0, 'H': 0, 'O': -1, 'N': 0},
    'A10': {'C': 0, 'H': 0, 'O': 0, 'N': 0}
}

# Define shift value of B parts
B_shift_set = {
    'B8': {'C': 0, 'H': 0, 'O': -1, 'N': 0},
    'B9': {'C': 0, 'H': -2, 'O': -2, 'N': 0},
    'B10': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'B11': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'B12': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'B13': {'C': 0, 'H': 0, 'O': 1, 'N': 0},
    'B14': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'B15': {'C': 0, 'H': 0, 'O': -1, 'N': 0},
    'B16': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'B17': {'C': 0, 'H': 0, 'O': -1, 'N': 0},
    'B18': {'C': 0, 'H': 0, 'O': 0, 'N': 0}
}

# Define shift value of C parts
C_shift_set = {
    'C3': {'C': 0, 'H': 0, 'O': -1, 'N': 0},
    'C4': {'C': 0, 'H': 0, 'O': -1, 'N': 0},
    'C5': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'C6': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'C7': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'C8': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'C9': {'C': 0, 'H': 0, 'O': -1, 'N': 0},
    'C10': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'C11': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
    'C12': {'C': 0, 'H': 0, 'O': 0, 'N': 0},
}

ref_feature_TypeII = [
    {'C': 20, 'H': 32, 'O': 6, 'N': 0},  # charge = +1, Theo.mass 369
    {'C': 20, 'H': 30, 'O': 5, 'N': 0},  # charge = +1, Theo.mass 350
    {'C': 20, 'H': 28, 'O': 4, 'N': 0},  # charge = +1, Theo.mass 333
    {'C': 20, 'H': 26, 'O': 3, 'N': 0},  # charge = +1, Theo.mass 314
    {'C': 20, 'H': 24, 'O': 2, 'N': 0},  # charge = +1, Theo.mass 296
    {'C': 20, 'H': 22, 'O': 1, 'N': 0},  # charge = +1, Theo.mass 278
    {'C': 20, 'H': 20, 'O': 0, 'N': 0},  # charge = +1, Theo.mass 260
]

path_csv = os.path.join('C:/Users/zhang/Desktop/taxane/', 'tax2_ABC_SMILES.csv')
df = pd.read_csv(path_csv)  # Read the CSV file using the pandas and store the data in a DataFrame object named df.

for index, row in df.iterrows():  # Iterate over each row in the DataFrame.
    new_feature_set = []  # Create an empty list to store the calculated new feature molecular formulas.
    for feature in ref_feature_TypeII:  # Traverse the reference feature set to calculate the new molecular formula.
        # Create a dictionary to store the calculated new feature molecular formulas.
        modified_feature = {
            'C': feature['C'] + A_shift_set[row['A parts']]['C'] + B_shift_set[row['B parts']]['C']
                + C_shift_set[row['C parts']]['C'],
            'H': feature['H'] + A_shift_set[row['A parts']]['H'] + B_shift_set[row['B parts']]['H']
                + C_shift_set[row['C parts']]['H'],
            'O': feature['O'] + A_shift_set[row['A parts']]['O'] + B_shift_set[row['B parts']]['O']
                + C_shift_set[row['C parts']]['O']
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

output_path = os.path.join('C:/Users/zhang/Desktop/taxane/', 'tax2_ABC_SMILES_feature'
                                                                                   'ion.csv')
df.to_csv(output_path, index=False)