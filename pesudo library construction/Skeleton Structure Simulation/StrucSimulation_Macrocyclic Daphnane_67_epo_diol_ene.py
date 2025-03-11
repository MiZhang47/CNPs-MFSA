"""
Last updated: 2024/12
Structural simulation for macrocyclic daphnane
"""

import csv
import os

''' Define the SMILES fragments for the A rings '''
A6 = "C(C)C(=O)"  # 1-alkyl-3-one A8
A7 = "C(O)(C)C(=O)"  # 1-alkyl-2-ol-3-one A9
A8 = "C(C)C(O*)"  # 1-alkyl-3-ol A10
A9 = "C(O)(C)C(O*)"  # 1-alkyl-2-ol-3-ol A11

''' Define the SMILES fragments for the A rings (seco ring type)'''
A10 = "C(C)C(=O)"  # 3,4-seco ring A12

''' Define the SMILES fragments for the A rings (bicyclo ring type)'''
A11 = "C7(OC6(C)OC7=O"  # bicyclo[2.2.1]heptane ring A13

''' Define the B ring SMILES fragments '''
# five-member ring type
# "xxx": macrocyclic ring; "yyy": A ring; "zzz": C ring.
B1_1 = "*OCC12OC1C1C3OxxxC6yyyC(O)(C2O)C6C1(O4)zzz"  # B: 5-ol-6,7-epoxy, format 1
B2_1 = "*OCC1(O)C(O)C2C3OxxxC6yyyC(O)(C1O)C6C2(O4)zzz"  # B: 5-ol-6,7-diol, format 1
B3_1 = "*OCC1=CC2C3OxxxC6yyyC(O)(C1O)C6C2(O4)zzz"  # B: 5-ol-6,7-ene, format 1
B4_1 = "*OCC12OC1C1C3OxxxC6yyyC(O)(C2)C6C1(O4)zzz"  # B: 5-dehydro-6,7-epoxy, format 1
B5_1 = "*OCC1(O)C(O)C2C3OxxxC6yyyC(O)(C1)C6C2(O4)zzz"  # B: 5-dehydro-6,7-diol, format 1
B6_1 = "*OCC1=CC2C3OxxxC6yyyC(O)(C1)C6C2(O4)zzz"  # B: 5-dehydro-6,7-ene, format 1
# seco ring type
# "xxx": macrocyclic ring; "mmm": A ring; "zzz": C ring.
B1_2 = "*OCC12OC1C1C3OxxxC(mmmO)C(C(=O)C2O)C1(O4)zzz"  # B: 5-ol-6,7-epoxy, format 2
B2_2 = "*OCC1(O)C(O)C2C3OxxxC(mmmO)C(C(=O)C1O)C2(O4)zzz"  # B: 5-ol-6,7-diol, format 2
B3_2 = "*OCC1=CC2C3OxxxC(mmmO)C(C(=O)C1O)C2(O4)zzz"  # B: 5-ol-6,7-ene, format 2
B4_2 = "*OCC12OC1C1C3OxxxC(mmmO)C(C(=O)C2)C1(O4)zzz"  # B: 5-dehydro-6,7-epoxy, format 2
B5_2 = "*OCC1(O)C(O)C2C3OxxxC(mmmO)C(C(=O)C1)C2(O4)zzz"  # B: 5-dehydro-6,7-diol, format 2
B6_2 = "*OCC1=CC2C3OxxxC(mmmO)C(C(=O)C1)C2(O4)zzz"  # B: 5-dehydro-6,7-ene, format 2
# Bicyclo ring type
# "xxx": macrocyclic ring; "nnn": A ring; "zzz": C ring.
B1_3 = "*OCC12OC1C1C3OxxxC6C(nnn)C2O)C1(O4)zzz"  # B: 5-ol-6,7-epoxy, format 3
B2_3 = "*OCC1(O)C(O)C2C3OxxxC6C(nnn)C1O)C2(O4)zzz"  # B: 5-ol-6,7-diol, format 3
B3_3 = "*OCC1=CC2C3OxxxC6C(nnn)C1O)C2(O4)zzz"  # B: 5-ol-6,7-ene, format 3
B4_3 = "*OCC12OC1C1C3OxxxC6C(nnn)C2)C1(O4)zzz"  # B: 5-dehydro-6,7-epoxy, format 3
B5_3 = "*OCC1(O)C(O)C2C3OxxxC6C(nnn)C1)C2(O4)zzz"  # B: 5-dehydro-6,7-diol, format 3
B6_3 = "*OCC1=CC2C3OxxxC6C(nnn)C1)C2(O4)zzz"  # B: 5-dehydro-6,7-ene, format 3

''' Define the SMILES fragments for the C rings '''
C1 = "C(C)CC3(C(=C)C)O5"  # 9,13,14-orthoester
C2 = "C(C)C(O*)C3(C(=C)C)O5"  # 9,13,14-orthoester, 12-ol
C5 = "C(CO*)CC3(C(=C)C)O5"  # 9,13,14-orthoester, 18-ol
C6 = "C(CO*)C(O*)C3(C(=C)C)O5"  # 9,13,14-orthoester, 12-ol-18-ol

''' Define the SMILES fragments for the Macrocyclic rings '''
# C9 macrocyclic ring (M9)
M1 = "C45CCCCCCCC"  # without acylation (C-1', 2', 3', 4', 5', 6', 7', 8', 9')

# C10 macrocyclic ring (M10)
M2 = "C45CCCCCCCC(C)"  # without acylation (C-1', 2', 3', 4', 5', 6', 7', 8', 9', (10'))
M3 = "C45C(O)CCCCCCC(C)"  # 2'-ol
M4 = "C45C(O)CCCCC(O*)CC(C)"  # 2'-ol-7'-ol
M5 = "C45C(O)CCCC(O*)C(O*)CC(C)"  # 2'-ol-6'-ol-7'-ol
M6 = "C45C(O)CCC(O*)C(O*)C(O*)CC(C)"  # 2'-ol-5'-ol-6'-ol-7'-ol

# C14 macrocyclic ring (M14)
M7 = "C45CCCC=CCCCCC(CCC)"  # 5',6'-ene (C-1', 2', 3', 4', 5', 6', 7', 8', 9', 10', 11', (12', 13', 14'))

# C16 macrocyclic ring (M16)
M8 = "C45C=CCCCCCCCCCCCC(C)"  # 2',3'-ene (C-1', 2', 3', 4', 5', 6', 7', 8', 9', 10', 11', 12', 13', 14', 15', (16'))
M9 = "C45CCCCCCCCCCCCCC(C)(O)"  # 15'-ol

''' Generate the final SMILES by combining the MDS, A, C, and M fragments '''
new_smiles = []
combinations = []

A_fragments = [A6, A7, A8, A9]
A10_fragments = [A10]
A11_fragments = [A11]
B_fragments = [B1_1, B2_1, B3_1, B4_1, B5_1, B6_1,
                B1_2, B2_2, B3_2, B4_2, B5_2, B6_2,
                B1_3, B2_3, B3_3, B4_3, B5_3, B6_3]
C_fragments = [C1, C2, C5, C6]
M_fragments = [M1, M2, M3, M4, M5, M6, M7, M8, M9]

''' Use enumerate to iterate over the B_fragments, C_fragments, and M_fragments lists. '''
for idx, B in enumerate(B_fragments, 1):
    for c_idx, C in enumerate(C_fragments, 1):
        for m_idx, M in enumerate(M_fragments, 1):
            # Checks whether the substring "yyy" is contained in the current element B.
            if "yyy" in B:
                for a_idx, A in enumerate(A_fragments, 1):
                    # If it contains "yyy", traverse the A_fragments list
                    # and replace "xxx", "yyy" and "zzz" in B with M, A and C to generate a new smiles string.
                    smiles = B.replace("xxx", M).replace("yyy", A).replace("zzz", C)
                    new_smiles.append(smiles)  # Add it to the new_smiles list.
                    # The corresponding combination information is added to the combinations list at the same time.
                    combinations.append((f"B{idx}", f"A{a_idx}", f"C{c_idx}", f"M{m_idx}"))
            # Checks whether the substring "mmm" is contained in the current element B.
            elif "mmm" in B:
                for a_idx, A10 in enumerate(A10_fragments, 1):
                    # If it contains "yyy", traverse the A_fragments list
                    # and replace "xxx", "nnn" and "zzz" in B with M, A5 and C to generate a new smiles string.
                    smiles = B.replace("xxx", M).replace("mmm", A10).replace("zzz", C)
                    new_smiles.append(smiles)  # Add it to the new_smiles list.
                    # The corresponding combination information is added to the combinations list at the same time.
                    combinations.append((f"B{idx}", f"A5{a_idx}", f"C{c_idx}", f"M{m_idx}"))
            # Checks whether the substring "nnn" is contained in the current element B.
            elif "nnn" in B:
                for a_idx, A11 in enumerate(A11_fragments, 1):
                    # If it contains "yyy", traverse the A_fragments list
                    # and replace "xxx", "nnn" and "zzz" in B with M, A6 and C to generate a new smiles string.
                    smiles = B.replace("xxx", M).replace("nnn", A11).replace("zzz", C)
                    new_smiles.append(smiles)  # Add it to the new_smiles list.
                    # The corresponding combination information is added to the combinations list at the same time.
                    combinations.append((f"B{idx}", f"A6{a_idx}", f"C{c_idx}", f"M{m_idx}"))

# Define a dictionary to store the parts that need to be replaced and their corresponding replacement values.
replacement_dict = {
    'A1': 'A6', 'A2': 'A7', 'A3': 'A8', 'A4': 'A9', 'A51': 'A10', 'A61': 'A11',
    'B1': 'B1', 'B2': 'B2', 'B3': 'B3', 'B4': 'B4', 'B5': 'B5', 'B6': 'B6',
    'B7': 'B1', 'B8': 'B2', 'B9': 'B3', 'B10': 'B4', 'B11': 'B5', 'B12': 'B6',
    'B13': 'B1', 'B14': 'B2', 'B15': 'B3', 'B16': 'B4', 'B17': 'B5', 'B18': 'B6',
    'C1': 'C1', 'C2': 'C2', 'C3': 'C5', 'C4': 'C6'
}

# Define the output path for the CSV file
path_csv = os.path.join('C:/Users/zhang/Desktop/StrucSimulation6/', 'MD_67_epo_diol_ene(1296).csv')

# Write the final SMILES, their numbers and the parts to the CSV file
# Open a file named path_csv using the with statement.
# 'w' means to open the file in write mode, and newline='' means not to add an extra newline when writing.
# csvfile is a file object used to refer to the file in subsequent operations.
with open(path_csv, 'w', newline='') as csvfile:
    # Create a CSV file writer (csv_writer) that will be used to write data to the CSV file.
    csv_writer = csv.writer(csvfile)
    # Write a list (header row) to a CSV file using the writerow method of csv_writer.
    csv_writer.writerow(["Number", "A parts", "B parts", "C parts", "Macro parts", "New SMILES"])
    # Use the enumerate function to iterate over the elements of the new_smiles and combinations lists.
    # (use the zip function to combine them).
    # The enumerate function returns a tuple containing two elements: one is the index of the element (starting from 1),
    # and another one is the element itself (a tuple containing smile and parts).
    # idx is the index of the element, smile and parts are the element itself.
    for idx, (smile, combination) in enumerate(zip(new_smiles, combinations), 1):
        # Use list comprehension and the replacement_dict dictionary to replace the parts information.
        replaced_parts = [replacement_dict.get(part, part) for part in combination]
        # Adjust the order of A parts and B parts
        replaced_parts = [replaced_parts[1], replaced_parts[0], replaced_parts[2], replaced_parts[3]]
        # Use the writerow method of csv_writer to write a list (rows of data) to a CSV file.
        csv_writer.writerow([idx, *replaced_parts, smile])