"""
Last updated: 2024/12
Structural simulation for normal non-macrocyclic daphnane
"""
import csv
import os

''' Define the SMILES fragments for the A rings '''
A1_1 = "C1C=C(C)C(=O)C1(O)"  # 1,2-en-3-one, format 1
A1_2 = "C2C=C(C)C(=O)C2(O)"  # 1,2-en-3-one, format 2

A2_1 = "C1CC(C)C(=O)C1(O)"  # 1,2-dihydro-3-one, format 1
A2_2 = "C2CC(C)C(=O)C2(O)"  # 1,2-dihydro-3-one, format 2

A3_1 = "C1C=C(C)C(O*)C1(O)"  # 1,2-en-3-ol, format 1
A3_2 = "C2C=C(C)C(O*)C2(O)"  # 1,2-en-3-ol, format 2

A4_1 = "C1CC(C)C(O*)C1(O)"  # 1,2-dihydro-3-ol, format 1
A4_2 = "C2CC(C)C(O*)C2(O)"  # 1,2-dihydro-3-ol, format 2

A5_1 = "C1=CC(C)C(=O)C1(O)"  # 1,10-en-3-one, format 1
A5_2 = "C2=CC(C)C(=O)C2(O)"  # 1,10-en-3-one, format 2

''' Define the SMILES fragments for the B rings '''
# xxx: C ring + A ring
B2 = "*OCC1(O)C(O)xxxC1O"  # 5-ol-6,7-diol
B5 = "*OCC1(O)C(O)xxxC1"  # 5-dehydro-6,7-diol

''' Define the SMILES fragments for the C rings '''
C1_1 = "C1C3OC4(*)OC3(C(=C)C)CC(C)C1(O4)"  # 9,13,14-orthoester, format 1
C1_2 = "C2C3OC4(*)OC3(C(=C)C)CC(C)C2(O4)"  # 9,13,14-orthoester, format 2

C2_1 = "C1C3OC4(*)OC3(C(=C)C)C(O*)C(C)C1(O4)"  # 12-ol and 9,13,14-orthoester, format 1
C2_2 = "C2C3OC4(*)OC3(C(=C)C)C(O*)C(C)C2(O4)"  # 12-ol and 9,13,14-orthoester, format 2

C3_1 = "C1C(O*)C(O)(C(=C)C)CC(C)C1(O)"  # 9,13,14-triol, format 1
C3_2 = "C2C(O*)C(O)(C(=C)C)CC(C)C2(O)"  # 9,13,14-triol, format 2

C4_1 = "C1C(O*)C(O)(C(=C)C)C(O*)C(C)C1(O)"  # 12-ol and 9,13,14-triol, format 1
C4_2 = "C2C(O*)C(O)(C(=C)C)C(O*)C(C)C2(O)"  # 12-ol and 9,13,14-triol, format 2

''' Create lists for A and C rings with the correct combinations '''
A = [
    [A1_1, A2_1, A3_1, A4_1, A5_1],
    [A1_2, A2_2, A3_2, A4_2, A5_2]
]

C = [
    [C1_1, C2_1, C3_1, C4_1],
    [C1_2, C2_2, C3_2, C4_2]
]

''' Generate the CA combinations '''
CA_combinations = []  # Create an empty list CA_combinations to store the resulting triples.
CA_parts = []  # Create an empty list CA_parts.

for i in range(2):  # Define a loop where i will iterate from 0 to 1.
    # Use the enumerate function to iterate over the ith element of list C.
    # The enumerate function returns a tuple containing two elements: one is the index of the element (starting from 1),
    # and another one is the element itself. j is the index of the element and c is the element itself.
    for j, c in enumerate(C[i], 1):
        # Use the enumerate function to iterate over the ith element of list A.
        # The enumerate function returns a tuple containing two elements:
        # one is the index of the element (starting from 1), and another one is the element itself.
        # k is the index of the element and a is the element itself.
        for k, a in enumerate(A[i], 1):
            # Add the triplet to the CA_combinations list. The triplet includes: the concatenation of c and a,
            # a string formatted as "C" followed by the index j of the element,
            # and a string formatted as "A" followed by the index k of the element.
            CA_combinations.append((c + a, "C{}".format(j), "A{}".format(k)))

''' Define the B ring fragments '''
B_fragments = [B2, B5]

''' Generate the final SMILES by combining the B and CA fragments '''
new_smiles = []  # Create an empty list new_smiles.
new_parts = []  # Create an empty list new_parts.

# Use the enumerate function to traverse the B_fragments list.
# The enumerate function returns a tuple containing two elements: one is the index of the element (starting from 1),
# and another one is the element itself. idx is the index of the element and B is the element itself.
for idx, B in enumerate(B_fragments, 1):
    # Checks whether the current element B contains both the characters '1' and '2'.
    # If included, execute the next block of code, otherwise skip to the else part.
    if '1' in B and '2' in B:
        # If B contains '1' and '2', iterate through the first 128 elements of the CA_combinations list.
        # These elements should be a triplet consisting of CA, c_part and a_part.
        for CA, c_part, a_part in CA_combinations[:20]:  # CxA=4x5=20
            # Replace the "xxx" string in B with CA, and add the replaced string to the new_smiles list.
            new_smiles.append(B.replace("xxx", CA))
            # Add a triple ("B{}".format(idx), a_part, c_part) to the new_parts list.
            # "B{}".format(idx) will return a string formatted as "B" followed by the index of the element.
            new_parts.append(("B{}".format(idx), a_part, c_part))
    else:  # If B does not contain both '1' and '2', then the next block of code is executed.
        # If B does not contain '1' and '2', traverse all elements after the 40th element of the CA_combinations list.
        # These elements should be a triplet consisting of CA, c_part and a_part.
        for CA, c_part, a_part in CA_combinations[20:]:  # CxA=4x5=20
            # Replace the "xxx" string in B with CA, and add the replaced string to the new_smiles list.
            new_smiles.append(B.replace("xxx", CA))
            # Add a triple ("B{}".format(idx), a_part, c_part) to the new_parts list.
            # "B{}".format(idx) will return a string formatted as "B" followed by the index of the element.
            new_parts.append(("B{}".format(idx), a_part, c_part))

# Define a dictionary to store the parts that need to be replaced and their corresponding replacement values.
replacement_dict = {'B1': 'B2', 'B2': 'B5'}

# Define the output path for the CSV file
path_csv = os.path.join('C:/Users/zhang/Desktop/StrucSimulation6/', 'D_67_diol(40).csv')

# Write the new SMILES and their parts to the CSV file.
with open(path_csv, 'w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    # Write the header of the CSV file.
    csv_writer.writerow(["Number", "A parts", "B parts", "C parts", "New SMILES"])
    # Iterate through the new_smiles and new_parts lists, obtaining the SMILES strings and parts information.
    for idx, (smile, parts) in enumerate(zip(new_smiles, new_parts), 1):
        # Use list comprehension and the replacement_dict dictionary to replace the parts information.
        replaced_parts = [replacement_dict.get(part, part) for part in parts]
        # Adjust the order of A parts and B parts
        replaced_parts = [replaced_parts[1], replaced_parts[0], replaced_parts[2]]
        # Write the replaced parts information and the corresponding SMILES strings to the CSV file.
        csv_writer.writerow([idx, *replaced_parts, smile])