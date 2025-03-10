import os
import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

input_csv_path = os.path.join('C:/Users/zhang/Desktop/','Plant_Skeleton_Percentage.csv')

data = pd.read_csv(input_csv_path)

plant_samples = [
    'Dac', 'Dat', 'Dau', 'Dax', 'Dbh', 'Dde', 'Dfe', 'Dgi', 'Dgr', 'Dho',
    'Djeft', 'Djef', 'Dki', 'Dob', 'Dodf', 'Dols', 'Dor', 'Doss', 'Dpa', 'Dpe',
    'Dpof', 'Dpol', 'Dpos', 'Dpsf', 'Dpsft', 'Dre', 'Dta', 'Dyu', 'Dgef', 'Dge',
    'Wal', 'Wca', 'Wde', 'Wdo', 'Wil', 'Wis', 'Wla', 'Wle', 'Wli', 'Wlr', 'Wlu',
    'Wmi', 'Wnu', 'Wpa', 'Wpi', 'Wsc', 'Wst', 'Wtr', 'Eal', 'Ecb', 'Ecf', 'Ecft',
    'Ecl', 'Ecs', 'Ega', 'Sch'
]

skeleton_types = [
    'A1B1C1', 'A4B2C1', 'A1B1C2', 'A4B1C1', 'A1B1C4', 'A4B7C3',
    'A8B1C1M3', 'A8B1C1M5', 'A11B1C1M2', 'A8B1C1M6', 'A8B1C1M2',
    'A11B1C1M3', 'A8B1C1M7', 'A11B1C1M4', 'A8B1C1M4', 'A11B1C1M6'
]

filtered_data = data[data.columns.intersection(plant_samples)]

numeric_data = filtered_data.fillna(0).astype(float)

cosine_sim_matrix = cosine_similarity(numeric_data.T)

similarity_df = pd.DataFrame(
    cosine_sim_matrix,
    index=plant_samples,
    columns=plant_samples
)

output_csv_path = os.path.join(
    'C:/Users/zhang/Desktop/Thyme56-centroid_rerun410/',
    'plant_sample_similarity_Percentage.csv'
)
similarity_df.to_csv(output_csv_path)