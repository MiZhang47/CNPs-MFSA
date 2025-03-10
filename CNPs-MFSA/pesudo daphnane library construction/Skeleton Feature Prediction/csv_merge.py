"""
Last updated: 2023/05/07 csv file merge
"""

import os
import pandas as pd

# Specify the directory path
path = 'C:/Users/zhang/Desktop/FeaturePredict2/'

# Get the file names of all csv files in the directory
csv_files = [f for f in os.listdir(path) if f.endswith('.csv')]

# Create an empty Excel file
writer = pd.ExcelWriter(path + 'StrucMSFeature2.xlsx', engine='xlsxwriter')

# Iterate through each csv file, and write its data to a separate sheet in the Excel file
for csv_file in csv_files:
    # Read the data from the csv file
    df = pd.read_csv(path + csv_file)

    # Use the file name as the sheet name
    sheet_name = csv_file[:-4]  # Remove the file extension '.csv'

    # Write the data to a sheet in the Excel file
    df.to_excel(writer, sheet_name=sheet_name, index=False)

# Save the Excel file
writer.save()