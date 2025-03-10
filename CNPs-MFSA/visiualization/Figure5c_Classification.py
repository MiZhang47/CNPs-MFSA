import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

input_csv_path = os.path.join('C:/Users/zhang/Desktop/',
                              'Thyme56_Annotation Result.csv')

data = pd.read_csv(input_csv_path)

type_counts = data['type'].value_counts()
skeleton_counts = data.groupby('type')['skeleton type'].value_counts()

inner_labels = type_counts.index
inner_sizes = type_counts.values

outer_labels = [idx[1] for idx in skeleton_counts.index]
outer_sizes = skeleton_counts.values

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12

inner_colors = ['#45505F', '#C5826F']
outer_colors = [
    '#7D95B0', '#8C8C8C', '#789AA5', '#D6D1C9', '#F0EFE9', '#D7B298',
    '#C2C5A8', '#9FC3B5', '#628476', '#F5E8D8', '#E9CB83', '#899CA3',
    '#A8BE85', '#CB783C', '#597964', '#156D88'
]

fig, ax = plt.subplots(figsize=(10, 5))

wedges_inner, texts_inner = ax.pie(
    inner_sizes, radius=0.8, labels=None,
    wedgeprops=dict(width=0.3, edgecolor='w'), colors=inner_colors[:len(inner_sizes)]
)

wedges_outer = ax.pie(
    outer_sizes, radius=1.1, labels=None, autopct=None,
    wedgeprops=dict(width=0.3, edgecolor='w'), colors=outer_colors[:len(outer_sizes)]
)

legend_labels = [f"{label} ({size})" for label, size in zip(outer_labels, outer_sizes)]
plt.legend(
    wedges_outer[0],
    legend_labels,
    loc="center left", bbox_to_anchor=(1.2, 0.5), fontsize=10
)

ax.set(aspect="equal")

plt.title("Skeleton Type Distribution", fontsize=14)

output_path = os.path.join('C:/Users/zhang/Desktop/', 'Figure5c_Classification.png')
plt.savefig(output_path, dpi=600, bbox_inches='tight')

plt.show()