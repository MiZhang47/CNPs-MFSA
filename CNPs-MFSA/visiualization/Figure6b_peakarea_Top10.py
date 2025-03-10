import pandas as pd
import os
import matplotlib.pyplot as plt

input_csv_path = 'C:/Users/zhang/Desktop/Thyme56_Annotation Result.csv'

data = pd.read_csv(input_csv_path)

plant_samples = [
    'Dac', 'Dat', 'Dau', 'Dax', 'Dbh', 'Dde', 'Dfe', 'Dgi', 'Dgr', 'Dho',
    'Djeft', 'Djef', 'Dki', 'Dob', 'Dodf', 'Dols', 'Dor', 'Doss', 'Dpa', 'Dpe',
    'Dpof', 'Dpol', 'Dpos', 'Dpsf', 'Dpsft', 'Dre', 'Dta', 'Dyu', 'Dgef', 'Dge',
    'Wal', 'Wca', 'Wde', 'Wdo', 'Wil', 'Wis', 'Wla', 'Wle', 'Wli', 'Wlr', 'Wlu',
    'Wmi', 'Wnu', 'Wpa', 'Wpi', 'Wsc', 'Wst', 'Wtr', 'Eal', 'Ecb', 'Ecf', 'Ecft',
    'Ecl', 'Ecs', 'Ega', 'Sch'
]

skeleton_types = data['skeleton type'].unique()

values = data.groupby(['skeleton type'])[plant_samples].sum()

color_map = {
    "A1B1C2": "#a6bddb",
    "A1B1C1": "#ccebc5",
    "A4B2C1": "#fdb462",
    "A4B1C1": "#fb8072",
    "A1B1C4": "#80b1d3",
    "A4B7C3": "#b3de69",
    "A8B1C1M2": "#fccde5",
    "A8B1C1M5": "#d9d9d9",
    "A8B1C1M3": "#bc80bd",
    "A8B1C1M4": "#ffed6f",
    "A8B1C1M6": "#8dd3c7",
    "A8B1C1M7": "#bebada",
    "A11B1C1M6": "#fb9a99",
    "A11B1C1M2": "#e31a1c",
    "A11B1C1M3": "#ff7f00",
    "A11B1C1M4": "#33a02c"
}

fig, ax = plt.subplots(figsize=(10, 8))

for i, skeleton in enumerate(skeleton_types):
    for j, plant in enumerate(plant_samples):
        value = values.loc[skeleton, plant] if skeleton in values.index else 0
        if value > 0:
            ax.scatter(
                i, j, color=color_map.get(skeleton, "#000000"),
                s=value * 10, alpha=0.8
            )

ax.set_xticks(range(len(skeleton_types)))
ax.set_xticklabels(skeleton_types, rotation=90, fontname='Arial', fontsize=12)
ax.set_yticks(range(len(plant_samples)))
ax.set_yticklabels(plant_samples, fontname='Arial', fontsize=12)

output_path = os.path.join('C:/Users/zhang/Desktop', 'Figure6b_peakarea_Top10.png')
plt.savefig(output_path, dpi=600, bbox_inches='tight')

plt.show()