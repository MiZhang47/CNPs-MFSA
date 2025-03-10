import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram

plt.rcParams['font.family'] = 'Arial'

output_csv_path = os.path.join('C:/Users/zhang/Desktop/', 'Plant_Skeleton_Percentage.csv')
result_df = pd.read_csv(output_csv_path, index_col=0)

input_csv_path = os.path.join('C:/Users/zhang/Desktop/', 'plant_sample_similarity_Percentage.csv')
similarity_df = pd.read_csv(input_csv_path, index_col=0)

threshold = 0.7
filtered_similarity = similarity_df.where(similarity_df > threshold, 0)
distance_matrix = 1 - filtered_similarity
distance_matrix.fillna(1, inplace=True)

linkage_matrix = linkage(distance_matrix, method="average")

dendro = dendrogram(linkage_matrix, labels=similarity_df.columns, no_plot=True)
ordered_samples = dendro["ivl"]

result_df = result_df[ordered_samples]

percent_df = result_df.div(result_df.sum(axis=0), axis=1) * 100

color_palette = [
    '#45505F', '#C5826F', '#789AA5', '#D6D1C9', '#F0EFE9', '#D7B298',
    '#C2C5A8', '#9FC3B5', '#628476', '#F5E8D8', '#E9CB83', '#899CA3',
    '#A8BE85', '#CB783C', '#597964', '#156D88'
]
color_map = {skeleton: color_palette[i % len(color_palette)] for i, skeleton in enumerate(percent_df.index)}

fig = plt.figure(figsize=(20, 10))

ax_dendro = plt.subplot2grid((2, 1), (0, 0), rowspan=1, colspan=1)
dendrogram(linkage_matrix, labels=similarity_df.columns, ax=ax_dendro, leaf_rotation=90, color_threshold=0.7)
ax_dendro.set_ylabel("Distance", fontsize=16)
ax_dendro.tick_params(axis="x", labelsize=14)
ax_dendro.tick_params(axis="y", labelsize=14)

ax_bar = plt.subplot2grid((2, 1), (1, 0), rowspan=1, colspan=1)
cumulative_bottom = [0] * len(percent_df.columns)
bar_width = 0.8
for skeleton in percent_df.index:
    ax_bar.bar(
        percent_df.columns,
        percent_df.loc[skeleton],
        bottom=cumulative_bottom,
        label=skeleton,
        color=color_map[skeleton],
        alpha=0.8,
        edgecolor="black",
        linewidth=0.5,
        width=bar_width
    )
    cumulative_bottom = [cumulative_bottom[i] + percent_df.loc[skeleton][i] for i in range(len(cumulative_bottom))]

ax_bar.set_ylabel("Percentage (%)", fontsize=16)
ax_bar.set_ylim(0, 100)
ax_bar.set_yticks(range(0, 101, 20))
ax_bar.tick_params(axis="x", labelsize=14, rotation=90)
ax_bar.tick_params(axis="y", labelsize=14)
ax_bar.legend(bbox_to_anchor=(0.5, -0.25), loc="upper center", ncol=8, fontsize=14)
ax_bar.grid(axis="y", linestyle="--", linewidth=0.7, alpha=0.7)
ax_bar.set_xlim(-0.5, len(percent_df.columns) - 0.5)

output_image_path = os.path.join('C:/Users/zhang/Desktop/', 'Figure5e5f_Hierarchical clustering and Proportional distribution .png')
plt.tight_layout()
plt.savefig(output_image_path)
plt.show()