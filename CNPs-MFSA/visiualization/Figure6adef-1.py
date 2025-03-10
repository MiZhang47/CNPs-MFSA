import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MultipleLocator, FuncFormatter

input_csv_path = 'C:/Users/zhang/Desktop/Thyme56-centroid_rerun410/Thyme56_Annotation Result_395_id_remove_adduct.csv'
df = pd.read_csv(input_csv_path)

plant_samples = [
    'Dac', 'Dat', 'Dau', 'Dax', 'Dbh', 'Dde', 'Dfe', 'Dgi', 'Dgr', 'Dho',
    'Djeft', 'Djef', 'Dki', 'Dob', 'Dodf', 'Dols', 'Dor', 'Doss', 'Dpa', 'Dpe',
    'Dpof', 'Dpol', 'Dpos', 'Dpsf', 'Dpsft', 'Dre', 'Dta', 'Dyu', 'Dgef', 'Dge',
    'Wal', 'Wca', 'Wde', 'Wdo', 'Wil', 'Wis', 'Wla', 'Wle', 'Wli', 'Wlr', 'Wlu',
    'Wmi', 'Wnu', 'Wpa', 'Wpi', 'Wsc', 'Wst', 'Wtr', 'Eal', 'Ecb', 'Ecf', 'Ecft',
    'Ecl', 'Ecs', 'Ega', 'Sch'
]
# Calculate the count of non-zero values for each compound across all plant samples
df['Count'] = df[plant_samples].astype(bool).sum(axis=1)

# Prepare the x and y values for the scatter plot
x_scatter = df['Compound Number']
y_scatter = df['Count']

# Create a custom colormap: light gray to deep blue
colors_list = ["white", "blue"]  # Define color gradient
cmap = LinearSegmentedColormap.from_list("custom_gray_to_blue", colors_list)
norm = mcolors.Normalize(vmin=y_scatter.min(), vmax=y_scatter.max())  # Normalize based on count range
colors = cmap(norm(y_scatter))

# Create a figure and define a GridSpec layout
fig = plt.figure(figsize=(12, 5))
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1.5], width_ratios=[4, 1], hspace=0.2, wspace=0.2)

# Set the global font to Arial
plt.rcParams['font.family'] = 'Arial'

fig = plt.figure(figsize=(16, 16))
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 4], width_ratios=[4, 1], hspace=0.05, wspace=0.1)

# 1. Top scatter plot showing distribution across plant samples
ax_top = fig.add_subplot(gs[0, 0])
scatter = ax_top.scatter(x_scatter, y_scatter, c=colors, cmap=cmap, edgecolor='black')

# Label compounds with occurrences greater than 15
for i, count in enumerate(y_scatter):
    if count > 25:
        ax_top.text(x_scatter[i], y_scatter[i] + 0.5, str(x_scatter[i]), fontsize=16, ha='center', va='bottom')

# Add colorbar for the counts
# cbar = plt.colorbar(scatter, ax=ax_top, label='Count of Non-Zero Entries')
# cbar.set_label('Number of Occurrences', fontsize=14)

# Set plot labels and title
# ax_top.set_xlabel('Compound Number', fontsize=14)
ax_top.set_ylabel('Count', fontsize=16)
# ax_top.set_title('Distribution of Compounds Across Plant Samples', fontsize=16)

# Set x-axis ticks to match the provided sample image
ax_top.set_xticks(range(0, int(df['Compound Number'].max()) + 40, 40))
ax_top.set_ylabel('Counts', fontsize=16)
ax_top.tick_params(axis='x', which='major', labelsize=14)  # Change 'labelsize' to your desired font size
ax_top.tick_params(axis='y', which='major', labelsize=14)  # Change 'labelsize' to your desired font size
# Set axis limits
ax_top.set_xlim(-5, 405)
ax_top.set_ylim(-1, 42)

# 2. Bottom scatter plot visualizing compound distribution across plant samples
# Process data for the scatter plot
long_df = df.melt(
    id_vars=['Compound Number', 'row m/z', 'elemental composition', 'adduct ion', 'row retention time', 'Skeleton SMILES', 'Acyl SMILES',
             'identification'],
    value_vars=plant_samples, var_name='Plant Sample', value_name='Concentration')

long_df = long_df[long_df['Concentration'] > 0]

# Set `Plant Sample` order as reversed to match the expected display order
reversed_order = list(reversed(plant_samples))
long_df['Plant Sample'] = pd.Categorical(long_df['Plant Sample'], categories=reversed_order, ordered=True)

# Create color mapping
color_map = {'possibly undescribed': 'red', 'possibly known': 'grey'}
long_df['Color'] = long_df['identification'].apply(
    lambda x: 'possibly undescribed' if x == 'possibly undescribed' else 'possibly known')

# Draw scatter plot
ax_scatter = fig.add_subplot(gs[1, 0])
sizes = long_df['Concentration'] / 1000000
colors = long_df['Color'].map(color_map)
scatter_plot = ax_scatter.scatter(long_df['Compound Number'], long_df['Plant Sample'], s=sizes, c=colors, alpha=0.6)

# Add legend for scatter plot
# handles = [plt.Line2D([0], [0], marker='o', color='w', label='Possibly Undescribed', markersize=10, markerfacecolor='red'),
           # plt.Line2D([0], [0], marker='o', color='w', label='Possibly Known', markersize=10, markerfacecolor='lightgrey')]
# ax_scatter.legend(handles=handles)

plant_samples_count = len(plant_samples)
padding = 1

ax_scatter.set_ylim(-padding, plant_samples_count)

ax_scatter.set_yticks(range(plant_samples_count))
ax_scatter.set_yticklabels(plant_samples, fontsize=12)

ax_scatter.tick_params(axis='x', which='major', labelsize=14)
ax_scatter.xaxis.set_major_locator(plt.MultipleLocator(40))
ax_scatter.set_xlim([-5, 405])

# Set labels for scatter plot
ax_scatter.set_xlabel('Compound Number', fontsize=16)
ax_scatter.set_ylabel('Plant Sample', fontsize=16)
# ax_scatter.set_title('Scatter Plot of Compound Concentrations in Plants', fontsize=16)
ax_scatter.tick_params(axis='both', which='major', labelsize=14)
ax_scatter.xaxis.set_major_locator(plt.MultipleLocator(40))
ax_scatter.set_xlim([-5, 405])

ax_scatter.grid(visible=True, which='major', axis='x', linestyle='--', color='grey', alpha=0.6)
ax_scatter.grid(visible=True, which='major', axis='y', linestyle='--', color='grey', alpha=0.6)

ax_scatter.tick_params(axis='both', which='major', labelsize=14)

# 3. Right bar plot visualizing compound distribution across plant samples
# Process data for the bar plot
filtered_df = df[plant_samples + ['identification']]
results = {'Sample': [], 'Total Compounds': [], 'Possibly New Compounds': []}

for sample in plant_samples:
    total_compounds = (filtered_df[sample] != 0).sum()
    possibly_new_compounds = ((filtered_df[sample] != 0) & (filtered_df['identification'] == 'possibly undescribed')).sum()
    results['Sample'].append(sample)
    results['Total Compounds'].append(total_compounds)
    results['Possibly New Compounds'].append(possibly_new_compounds)

results_df = pd.DataFrame(results)

# Sort results to match plant sample order
results_df['Sample'] = pd.Categorical(results_df['Sample'], categories=plant_samples, ordered=True)
results_df = results_df.sort_values(by='Sample')

# Draw bar plot
ax_bar = fig.add_subplot(gs[1, 1])
ax_bar.barh(results_df['Sample'], results_df['Total Compounds'], color='lightgrey', edgecolor='white') #, label='Total Compounds')
ax_bar.barh(results_df['Sample'], results_df['Possibly New Compounds'], color='#e05759', edgecolor='white') #, label='Possibly New Compounds')

# Annotate each bar with compound counts
for index, (total, new) in enumerate(zip(results_df['Total Compounds'], results_df['Possibly New Compounds'])):
    ax_bar.text(total + 1, index, str(total), va='center', ha='left', color='black')
    ax_bar.text(new + 1, index, str(new), va='center', ha='left', color='black')

ax_bar.set_ylim(-padding, plant_samples_count)

ax_bar.set_yticks(range(plant_samples_count))
ax_bar.set_yticklabels(plant_samples, fontsize=12)

# Set labels for bar plot
ax_bar.set_xlabel('Counts', fontsize=16)
# ax_bar.set_title('Total and Possibly New Compounds in Each Plant', fontsize=16)
ax_bar.set_yticks(range(len(plant_samples)))
ax_bar.set_yticklabels(plant_samples)
ax_bar.tick_params(axis='both', which='major', labelsize=14)
ax_bar.legend(loc='lower right')
ax_bar.set_xlim([0, 130])

plt.tight_layout()

output_image_path = 'C:/Users/zhang/Desktop/Figure6adef-1.png'
plt.savefig(output_image_path, dpi=600, bbox_inches='tight')
print(f"图片已保存至: {output_image_path}")

plt.show()
