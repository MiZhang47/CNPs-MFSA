"""
Last updated: 2024/12
Structural simulation for normal non-macrocyclic daphnane
"""
import csv
import os

''' Define the SMILES fragments for the A rings '''
A4_3 = "(C(O*)C(C)CC3"  # 1,2-dihydro-3-ol, format 3

''' Define the SMILES fragments for the B rings '''
# xxx: A ring + C ring
B7 = "*OCC1(O)C2OC3xxxC1O"  # 5-ol-4,7-epoxy
B8 = "*OCC12OC3xxxC2O"  # 5-ol-4,6-epoxy

''' Define the SMILES fragments for the C rings '''
C3_3 = "C3(O)C(C)CC(O)(C(=C)C)C(O*)C23)"  # 9,13,14-triol, format 3
C3_4 = "C3(O)C(C)CC(O)(C(=C)C)C(O*)C3C1O)"  # 9,13,14-triol, format 4

C4_3 = "C3(O)C(C)C(O*)C(O)(C(=C)C)C(O*)C23)"  # 12-ol and 9,13,14-triol, format 3
C4_4 = "C3(O)C(C)C(O*)C(O)(C(=C)C)C(O*)C3C1O)"  # 12-ol and 9,13,14-triol, format 4

''' Create lists for A and C rings with the correct combinations '''
A = [[A4_3],
    [A4_3]
]

C = [
    [C3_3, C4_3],
    [C3_4, C4_4]
]

''' Generate the AC combinations '''
AC_combinations = []  # Create an empty list AC_combinations to store the resulting triples.
AC_parts = []  # Create an empty list AC_parts.

for i in range(2):  # Define a loop where i will iterate from 0 to 1.
    # Use the enumerate function to iterate over the ith element of list A.
    # The enumerate function returns a tuple containing two elements: one is the index of the element (starting from 1),
    # and another one is the element itself. j is the index of the element and a is the element itself.
    for j, a in enumerate(A[i], 1):
        # Use the enumerate function to iterate over the ith element of list C.
        # The enumerate function returns a tuple containing two elements:
        # one is the index of the element (starting from 1), and another one is the element itself.
        # k is the index of the element and c is the element itself.
        for k, c in enumerate(C[i], 1):
            # Add the triplet to the AC_combinations list. The triplet includes: the concatenation of a and c,
            # c string formatted as "A" followed by the index k of the element,
            # and a string formatted as "A" followed by the index j of the element.
            AC_combinations.append((a + c, "A{}".format(j), "C{}".format(k)))

''' Define the B ring fragments '''
B_fragments = [B7, B8]

''' Generate the final SMILES by combining the B and CA fragments '''
new_smiles = []  # Create an empty list new_smiles.
new_parts = []  # Create an empty list new_parts.

# Use the enumerate function to traverse the B_fragments list.
# The enumerate function returns a tuple containing two elements: one is the index of the element (starting from 1),
# and another one is the element itself. idx is the index of the element and B is the element itself.
for idx, B in enumerate(B_fragments, 1):
    # Checks whether the B string contains two '1'.
    # If included, execute the next block of code, otherwise skip to the else part.
    if B.count('1') == 2:
        # If B contains two '1', iterate through the first 128 elements of the AC_combinations list.
        # These elements should be a triplet consisting of AC, a_part and c_part.
        for AC, a_part, c_part in AC_combinations[:2]:  # AxC=1x2=2
            # Replace the "xxx" string in B with AC, and add the replaced string to the new_smiles list.
            new_smiles.append(B.replace("xxx", AC))
            # Add a triple ("B{}".format(idx), a_part, c_part) to the new_parts list.
            # "B{}".format(idx) will return a string formatted as "B" followed by the index of the element.
            new_parts.append(("B{}".format(idx), a_part, c_part))
    else:  # If B does not contain two '1', then the next block of code is executed.
        # If B does not contain two '1', traverse all elements after the 128th element of the AC_combinations list.
        # These elements should be a triplet consisting of AC, a_part and c_part.
        for AC, a_part, c_part in AC_combinations[2:]:  # AxC=1x2=2
            # Replace the "xxx" string in B with AC, and add the replaced string to the new_smiles list.
            new_smiles.append(B.replace("xxx", AC))
            # Add a triple ("B{}".format(idx), a_part, c_part) to the new_parts list.
            # "B{}".format(idx) will return a string formatted as "B" followed by the index of the element.
            new_parts.append(("B{}".format(idx), a_part, c_part))

# Define a dictionary to store the parts that need to be replaced and their corresponding replacement values.
replacement_dict = {'A1': 'A4', 'B1': 'B7', 'B2': 'B8', 'C1': 'C3', 'C2': 'C4'}

# Define the output path for the CSV file
path_csv = os.path.join('C:/Users/zhang/Desktop/StrucSimulation6/', 'D_47_46_epo(4).csv')

# Write the new SMILES and their parts to the CSV file
# Open a file named path_csv using the with statement.
# 'w' means to open the file in write mode, and newline='' means not to add an extra newline when writing.
# csvfile is a file object used to refer to the file in subsequent operations.
with open(path_csv, 'w', newline='') as csvfile:
    # Create a CSV file writer (csv_writer) that will be used to write data to the CSV file.
    csv_writer = csv.writer(csvfile)
    # Write a list (header row) to a CSV file using the writerow method of csv_writer.
    csv_writer.writerow(["Number", "A parts", "B parts", "C parts", "New SMILES"])
    # Iterate through the new_smiles and new_parts lists, obtaining the SMILES strings and parts information.
    for idx, (smile, parts) in enumerate(zip(new_smiles, new_parts), 1):
        # Use list comprehension and the replacement_dict dictionary to replace the parts information.
        replaced_parts = [replacement_dict.get(part, part) for part in parts]
        # Adjust the order of A parts and B parts
        replaced_parts = [replaced_parts[1], replaced_parts[0], replaced_parts[2]]
        # Write the replaced parts information and the corresponding SMILES strings to the CSV file.
        csv_writer.writerow([idx, *replaced_parts, smile])