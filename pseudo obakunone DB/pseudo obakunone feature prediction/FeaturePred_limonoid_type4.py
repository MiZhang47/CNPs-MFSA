"""
Last updated: 2025/05
Structure-based MS feature prediction for pDB-limonoid type I
"""

import os
import pandas as pd

# Define shift value of A parts
A_shift_set = {
    'A6': {'C': 0, 'H': 0, 'O': 0},
    'A7': {'C': 0, 'H': 2, 'O': 0},
    'A8': {'C': 2, 'H': 4, 'O': 2},
    'A9': {'C': 0, 'H': 2, 'O': 1},
}

# Define shift value of B parts
B_shift_set = {
    'B1': {'C': 2, 'H': 4, 'O': 1},  # 7-OAc
    'B2': {'C': 0, 'H': 2, 'O': 0},  # 7-ol
    'B3': {'C': 0, 'H': 0, 'O': 0}  # 7-one
}

# Define shift value of C parts
C_shift_set = {
    'C1': {'C': 0, 'H': 0, 'O': 0},  # no sub
}

# Define shift value of D parts
D_shift_set = {
    'D7': {'C': 0, 'H': 0, 'O': 0},
    'D8': {'C': 0, 'H': 0, 'O': 2}
}

ref_feature_TypeI = [
    {'C': 5,  'H': 4,  'O': 1},  # charge = +1, Theo.mass 81.03406525
    {'C': 5,  'H': 2,  'O': 2},  # charge = +1, Theo.mass 95.01320648
    {'C': 6,  'H': 6,  'O': 1},  # charge = +1, Theo.mass 95.04955292
    {'C': 9,  'H': 10, 'O': 1},  # charge = +1, Theo.mass 135.0805206
    {'C': 9,  'H': 12, 'O': 2},  # charge = +1, Theo.mass 153.0908813
    {'C': 10, 'H': 8,  'O': 2},  # charge = +1, Theo.mass 161.0595398
    {'C': 10, 'H': 10, 'O': 2},  # charge = +1, Theo.mass 163.0753021
    {'C': 11, 'H': 10, 'O': 2},  # charge = +1, Theo.mass 175.0752869
    {'C': 12, 'H': 10, 'O': 2},  # charge = +1, Theo.mass 187.0752258
    {'C': 14, 'H': 12, 'O': 1},  # charge = +1, Theo.mass 197.0958862
    {'C': 14, 'H': 14, 'O': 2},  # charge = +1, Theo.mass 215.1062469
    {'C': 15, 'H': 14, 'O': 3},  # charge = +1, Theo.mass 243.1014252
    {'C': 21, 'H': 18, 'O': 2},  # charge = +1, Theo.mass 303.1378174
    {'C': 22, 'H': 18, 'O': 2},  # charge = +1, Theo.mass 315.1376343
    {'C': 21, 'H': 20, 'O': 3},  # charge = +1, Theo.mass 321.1479187
    {'C': 22, 'H': 18, 'O': 3},  # charge = +1, Theo.mass 331.1323242
    {'C': 22, 'H': 20, 'O': 3},  # charge = +1, Theo.mass 333.1480103
    {'C': 23, 'H': 16, 'O': 3},  # charge = +1, Theo.mass 341.1168823
    {'C': 22, 'H': 20, 'O': 4},  # charge = +1, Theo.mass 349.1431885
    {'C': 23, 'H': 18, 'O': 4},  # charge = +1, Theo.mass 359.1274414
    {'C': 23, 'H': 20, 'O': 5},  # charge = +1, Theo.mass 377.1371765
    {'C': 25, 'H': 26, 'O': 4},  # charge = +1, Theo.mass 391.1894531
    {'C': 25, 'H': 28, 'O': 4},  # charge = +1, Theo.mass 393.2051086
    {'C': 25, 'H': 28, 'O': 5},  # charge = +1, Theo.mass 409.2002869
    {'C': 25, 'H': 30, 'O': 5},  # charge = +1, Theo.mass 411.2156067
    {'C': 26, 'H': 26, 'O': 5},  # charge = +1, Theo.mass 419.1845398
    {'C': 26, 'H': 28, 'O': 6},  # charge = +1, Theo.mass 437.1952515
    {'C': 26, 'H': 30, 'O': 7},  # charge = +1, Theo.mass 455.2057495
]


path_csv = os.path.join('C:/Users/zhang/Desktop/limonoid/', 'lim_oba_ABCD_SMILES.csv')
df = pd.read_csv(path_csv)  # Read the CSV file using the pandas and store the data in a DataFrame object named df.

for index, row in df.iterrows():  # Iterate over each row in the DataFrame.
    new_feature_set = []  # Create an empty list to store the calculated new feature molecular formulas.
    for feature in ref_feature_TypeI:  # Traverse the reference feature set to calculate the new molecular formula.
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

output_path = os.path.join('C:/Users/zhang/Desktop/limonoid/', 'lim_oba_ABCD_SMILES_feature'
                                                                                   'ion.csv')
df.to_csv(output_path, index=False)