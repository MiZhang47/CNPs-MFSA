import pandas as pd
import os

input_csv_path = os.path.join('C:/Users/zhang/Desktop/', 'Thyme56_Annotation Result.csv')

output_csv_path = os.path.join('C:/Users/zhang/Desktop/', 'Plant_Skeleton_Percentage.csv')

plant_samples = [
    'Dac', 'Dat', 'Dau', 'Dax', 'Dbh', 'Dde', 'Dfe', 'Dgi', 'Dgr', 'Dho',
    'Djeft', 'Djef', 'Dki', 'Dob', 'Dodf', 'Dols', 'Dor', 'Doss', 'Dpa', 'Dpe',
    'Dpof', 'Dpol', 'Dpos', 'Dpsf', 'Dpsft', 'Dre', 'Dta', 'Dyu', 'Dgef', 'Dge',
    'Wal', 'Wca', 'Wde', 'Wdo', 'Wil', 'Wis', 'Wla', 'Wle', 'Wli', 'Wlr', 'Wlu',
    'Wmi', 'Wnu', 'Wpa', 'Wpi', 'Wsc', 'Wst', 'Wtr', 'Eal', 'Ecb', 'Ecf', 'Ecft',
    'Ecl', 'Ecs', 'Ega', 'Sch'
]

df = pd.read_csv(input_csv_path)

skeleton_col = 'skeleton type'

all_skeleton_types = df[skeleton_col].unique()

result_df = pd.DataFrame(index=all_skeleton_types)

for plant in plant_samples:
    if plant in df.columns:
        counts = df[df[plant] > 0].groupby(skeleton_col).size()
        result_df[plant] = counts.reindex(all_skeleton_types, fill_value=0)
        total_count = result_df[plant].sum()
        if total_count > 0:
            result_df[plant] = result_df[plant] / total_count

result_df.index.name = skeleton_col
result_df.to_csv(output_csv_path)