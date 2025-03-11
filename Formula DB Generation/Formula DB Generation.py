"""
Last updated: 2025/01
Generation molecular formula and save to SQLite in positive & negative ion modes (M, M+NH3, M+HCOOH)
# Number of ion formula in Positive ion database: 63228
# Number of ion formula in Negative ion database: 63228
"""

import sqlite3

# Calculate the Index of Hydrogen Deficiency (IHD) for an organic compound
# based on the number of carbons (nC), hydrogens (nH), and oxygens (nO)
def calculate_IHD(nC, nH, nO):
    # Formula for calculating IHD, which measures the degree of unsaturation in a molecule
    return 1 + (nC + nO - nH / 2) / 2

# Generate a list of organic compounds based on possible combinations of C, H, and O
def generate_organic_compounds():
    compounds = []  # List to store the generated compounds
    # Iterate through possible numbers of carbon atoms (from 1 to 60)
    for nC in range(1, 61):
        # Iterate through possible numbers of hydrogen atoms (from 2 to 80, in steps of 2)
        for nH in range(2, 81, 2):
            # Iterate through possible numbers of oxygen atoms (from 1 to 20)
            for nO in range(1, 21):
                IHD = calculate_IHD(nC, nH, nO)  # Calculate the IHD for the compound
                # Only store compounds with an IHD between 3 and 25 (inclusive)
                if 3 <= IHD <= 25:
                    compound = (nC, nH, nO)
                    compounds.append(compound)
    return compounds

# Simulate the addition of functional groups to a list of compounds
def add_additional_groups(compounds):
    # Add an NH3 group to each compound
    nh3_added = [(nC, nH + 3, nO, 1) for nC, nH, nO in compounds]
    # Add a CH2O2 group to each compound
    ch2o2_added = [(nC + 1, nH + 2, nO + 2) for nC, nH, nO in compounds]
    # Add a CH3NH2 group to each compound
    ch5n_added = [(nC + 1, nH + 5, nO, 1) for nC, nH, nO in compounds]
    return nh3_added, ch2o2_added, ch5n_added

# Create a SQLite database and set up a table for storing compounds
def create_database(db_name):
    # Connect to the specified SQLite database file or create it if it doesn't exist
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    # Create a table named 'compounds' with a single column 'formula' if it does not already exist
    c.execute('''CREATE TABLE IF NOT EXISTS compounds (formula TEXT)''')
    conn.commit()  # Commit the changes to the database
    return conn  # Return the connection object for further operations

# Save a list of compounds and their NH3-modified versions into a database
def save_to_positive_ion_db(conn, compounds, nh3_added, ch5n_added):
    c = conn.cursor()  # Create a cursor object to interact with the database

    # Loop through each compound in the list and save it to the database
    for nC, nH, nO in compounds:
        formula = f"C{nC}H{nH}O{nO}"  # Format the formula as a string
        c.execute("INSERT INTO compounds (formula) VALUES (?)", (formula,))  # Execute SQL command to insert the formula

    # Loop through each NH3-modified compound in the list and save it to the database
    for nC, nH, nO, nN in nh3_added:
        formula = f"C{nC}H{nH}O{nO}N{nN}"  # Format the modified formula with NH3 group
        c.execute("INSERT INTO compounds (formula) VALUES (?)",
                (formula,))  # Execute SQL command to insert the modified formula

    # Loop through each CH3NH2-modified compound in the list and save it to the database
    for nC, nH, nO, nN in ch5n_added:
        formula = f"C{nC}H{nH}O{nO}N{nN}"  # Format the modified formula with NH3 group
        c.execute("INSERT INTO compounds (formula) VALUES (?)",
                (formula,))  # Execute SQL command to insert the modified formula

    conn.commit()  # Commit all changes to the database to make sure they are saved


# Save a list of compounds and their CH2O2-modified versions into a database
def save_to_negative_ion_db(conn, compounds, ch2o2_added):
    c = conn.cursor()  # Create a cursor object to interact with the database

    # Loop through each compound in the list and save it to the database
    for nC, nH, nO in compounds:
        formula = f"C{nC}H{nH}O{nO}"  # Format the formula as a string
        c.execute("INSERT INTO compounds (formula) VALUES (?)", (formula,))  # Execute SQL command to insert the formula

    # Loop through each CH2O2-modified compound in the list and save it to the database
    for nC, nH, nO in ch2o2_added:
        formula = f"C{nC}H{nH}O{nO}"  # Format the modified formula with CH2O2 group
        c.execute("INSERT INTO compounds (formula) VALUES (?)",
                (formula,))  # Execute SQL command to insert the modified formula

    conn.commit()  # Commit all changes to the database to make sure they are saved
# Generate molecular formulas
compounds = generate_organic_compounds()

# Add additional NH3 and HCOOH
nh3_added, ch2o2_added, ch5n_added = add_additional_groups(compounds)

# Create Positive Ion database and table
positive_ion_conn = create_database('Positive Ion Formula.db')
# Save to Positive Ion SQLite database
save_to_positive_ion_db(positive_ion_conn, compounds, nh3_added, ch5n_added)
# Close the Positive Ion database connection
positive_ion_conn.close()

# Create Negative Ion database and table
negative_ion_conn = create_database('Negative Ion Formula.db')
# Save to Negative Ion SQLite database
save_to_negative_ion_db(negative_ion_conn, compounds, ch2o2_added)
# Close the Negative Ion database connection
negative_ion_conn.close()