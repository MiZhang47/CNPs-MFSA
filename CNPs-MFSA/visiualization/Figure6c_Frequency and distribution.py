import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'Arial'

input_csv_path = os.path.join('C:/Users/zhang/Desktop/', 'Thyme56_Annotation Result.csv')

output_image_path = os.path.join('C:/Users/zhang/Desktop/', 'Figure6c_Frequency and distribution.png')

plant_samples = [
    'Dac', 'Dat', 'Dau', 'Dax', 'Dbh', 'Dde', 'Dfe', 'Dgi', 'Dgr', 'Dho',
    'Djeft', 'Djef', 'Dki', 'Dob', 'Dodf', 'Dols', 'Dor', 'Doss', 'Dpa', 'Dpe',
    'Dpof', 'Dpol', 'Dpos', 'Dpsf', 'Dpsft', 'Dre', 'Dta', 'Dyu', 'Dgef', 'Dge',
    'Wal', 'Wca', 'Wde', 'Wdo', 'Wil', 'Wis', 'Wla', 'Wle', 'Wli', 'Wlr', 'Wlu',
    'Wmi', 'Wnu', 'Wpa', 'Wpi', 'Wsc', 'Wst', 'Wtr', 'Eal', 'Ecb', 'Ecf', 'Ecft',
    'Ecl', 'Ecs', 'Ega', 'Sch'
]

try:
    df = pd.read_csv(input_csv_path)

    plant_sample_columns = [col for col in plant_samples if col in df.columns]

    ranges = [
        (1e5, 1e6),
        (1e6, 1e7),
        (1e7, 1e8),
        (1e8, 1e9),
        (1e9, 1e10)
    ]
    counts = []

    for lower, upper in ranges:
        count = df[plant_sample_columns].astype(float).applymap(lambda x: lower <= x < upper).sum().sum()
        counts.append(count)

    bins = [1e5, 1e6, 1e7, 1e8, 1e9, 1e10]

    plt.figure(figsize=(1.8, 4))
    plt.hist(
        bins[:-1], bins=bins, weights=counts, color='lightgray', edgecolor='black', rwidth=1.0
    )

    plt.xlabel("Peak Area", fontsize=12)
    plt.ylabel("Count", fontsize=12)
    plt.xscale('log')
    plt.xticks(bins, ['E5', '6', '7', '8', '9', '10'], fontsize=11)
    plt.yticks(fontsize=12)
    plt.ylim(0, 1100)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()
    plt.savefig(output_image_path, dpi=600)
    plt.show()

    print(f"{output_image_path}")

except Exception as e:
    print(f"{e}")