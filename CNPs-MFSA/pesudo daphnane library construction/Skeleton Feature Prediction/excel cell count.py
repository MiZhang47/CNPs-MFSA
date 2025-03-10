import pandas as pd

# Read the Excel file
file_path = "C:/Users/zhang/Desktop/FeaturePredict2/StrucMSFeature2.xlsx"
xl = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")

# Initialize an empty list to store the data from "MS Feature" columns
ms_feature_data = []

# Iterate through all sheets and extract the columns containing "MS Feature"
for sheet_name, sheet_data in xl.items():
    for col_name in sheet_data.columns:
        if "MS Feature" in col_name:
            ms_feature_data.extend(sheet_data[col_name].dropna().tolist())

# Perform statistics on the data
total_cells = len(ms_feature_data)
unique_cells = len(set(ms_feature_data))

# Remove duplicates while preserving the order
seen = set()
ms_feature_data_no_duplicates = [x for x in ms_feature_data if not (x in seen or seen.add(x))]
cells_after_removing_duplicates = len(ms_feature_data_no_duplicates)

# Assign the resulting list to the variable feature_formulas
feature_formulas = ms_feature_data_no_duplicates

C30_feature = []
if "C29" and "C29" "C29" and "C30" in feature_formulas:
    C30_feature.append(feature_formulas)
elif "C34" and "C30" in feature_formulas:

    print(f"Total cells: {total_cells}")
    print(f"Unique cells: {unique_cells}")
    print(f"Cells after removing duplicates (preserving one): {cells_after_removing_duplicates}")

# Print the feature_formulas list
    print("feature_formulas = ", feature_formulas)