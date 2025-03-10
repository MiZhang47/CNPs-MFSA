import os
import pandas as pd
from venn import venn
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

input_csv_path = os.path.join('C:/Users/zhang/Desktop/',
                              'Thyme56_Annotation Result.csv')
df = pd.read_csv(input_csv_path)

plt.rcParams['font.family'] = 'Arial'

plant_samples = [
    'Dac', 'Dat', 'Dau', 'Dax', 'Dbh', 'Dde', 'Dfe', 'Dgi', 'Dgr', 'Dho',
    'Djeft', 'Djef', 'Dki', 'Dob', 'Dodf', 'Dols', 'Dor', 'Doss', 'Dpa', 'Dpe',
    'Dpof', 'Dpol', 'Dpos', 'Dpsf', 'Dpsft', 'Dre', 'Dta', 'Dyu', 'Dgef', 'Dge',
    'Wal', 'Wca', 'Wde', 'Wdo', 'Wil', 'Wis', 'Wla', 'Wle', 'Wli', 'Wlr', 'Wlu',
    'Wmi', 'Wnu', 'Wpa', 'Wpi', 'Wsc', 'Wst', 'Wtr', 'Eal', 'Ecb', 'Ecf', 'Ecft',
    'Ecl', 'Ecs', 'Ega', 'Sch'
]

group_D = [col for col in plant_samples if col.startswith('D')]
group_W = [col for col in plant_samples if col.startswith('W')]
group_E = [col for col in plant_samples if col.startswith('E')]
group_S = [col for col in plant_samples if col.startswith('S')]

compounds_D = set(df.loc[df[group_D].sum(axis=1) > 0, 'Compound Number'])
compounds_W = set(df.loc[df[group_W].sum(axis=1) > 0, 'Compound Number'])
compounds_E = set(df.loc[df[group_E].sum(axis=1) > 0, 'Compound Number'])
compounds_S = set(df.loc[df[group_S].sum(axis=1) > 0, 'Compound Number'])

data = {
    "Daphne": compounds_D,
    "Wikstroemia": compounds_W,
    "Edgeworthia": compounds_E,
    "Stellera": compounds_S,
}

outline_colors = [
    "#354D6F",
    "#D9775B",
    "#60A6BD",
    "#F2CF7A",
]

fill_colors = [
    (53/255, 77/255, 111/255, 0.3),
    (217/255, 119/255, 91/255, 0.3),
    (96/255, 166/255, 189/255, 0.3),
    (242/255, 207/255, 122/255, 0.3),
]

plt.figure(figsize=(10, 10))
venn_diagram = venn(data, fontsize=16)

patches = venn_diagram.patches
for i, patch in enumerate(patches):
    if patch:
        patch.set_facecolor(fill_colors[i])
        patch.set_edgecolor(outline_colors[i % len(outline_colors)])
        patch.set_linewidth(1.5)

legend_labels = ["Daphne", "Wikstroemia", "Edgeworthia", "Stellera"]
legend_patches = [Patch(facecolor=fill_colors[i][:3] + (1,), edgecolor=outline_colors[i], label=legend_labels[i]) for i in range(len(legend_labels))]
plt.legend(handles=legend_patches, loc='upper right', title="Groups", fontsize=12, title_fontsize=14, frameon=True)

output_image_path = os.path.join('C:/Users/zhang/Desktop/', 'Figure5d_Distribution across four genera.png')
plt.savefig(output_image_path, dpi=600, bbox_inches='tight')

plt.show()