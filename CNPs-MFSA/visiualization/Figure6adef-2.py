import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FuncFormatter

input_csv_path = 'C:/Users/zhang/Desktop/Thyme56_Annotation Result.csv'
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

fig = plt.figure(figsize=(16, 16))
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 4], width_ratios=[4, 1], hspace=0.05, wspace=0.1)

# 1. Top scatter plot
ax_top = fig.add_subplot(gs[0, 0])

df['Total Peak Area'] = df[plant_samples].sum(axis=1)

colors = ['#e05759' if 'possibly undescribed' in id_ else 'grey' for id_ in df['identification']]
bars = ax_top.bar(df['Compound Number'], df['Total Peak Area'], color=colors)

ax_top.set_ylabel('Total Peak Area', fontsize=16)
ax_top.set_xticks(range(0, int(df['Compound Number'].max()) + 50, 50))
ax_top.tick_params(axis='both', which='major', labelsize=14)
ax_top.xaxis.set_major_locator(plt.MultipleLocator(40))
ax_top.set_xlim([-5, 405])

ax_top.set_ylim([0, 8.8e9])

ax_top.yaxis.set_major_locator(MultipleLocator(1e9))
ax_top.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x/1e9:.0f}E9'))

for i, bar in enumerate(bars):
    if df['Total Peak Area'].iloc[i] > 1e9:
        ax_top.text(
            df['Compound Number'].iloc[i],
            df['Total Peak Area'].iloc[i] + 1e8,
            str(df['Compound Number'].iloc[i]),
            fontsize=14, color='black', ha='center', va='bottom'
        )

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
ax_scatter.tick_params(axis='both', which='major', labelsize=14)
ax_scatter.xaxis.set_major_locator(plt.MultipleLocator(40))
ax_scatter.set_xlim([-5, 405])

ax_scatter.grid(visible=True, which='major', axis='x', linestyle='--', color='grey', alpha=0.6)
ax_scatter.grid(visible=True, which='major', axis='y', linestyle='--', color='grey', alpha=0.6)

ax_scatter.tick_params(axis='both', which='major', labelsize=14)

# 3. Right bar plot
ax_bar = fig.add_subplot(gs[1, 1])

total_peak_area = df[plant_samples].sum(axis=0)

possibly_undescribed = df[df['identification'] == 'possibly undescribed'][plant_samples].sum(axis=0)

ax_bar.barh(plant_samples, total_peak_area, color='lightgrey')

ax_bar.barh(plant_samples, possibly_undescribed, color='#e05759')

ax_bar.set_xlabel('Total Peak Area', fontsize=16)
ax_bar.legend(loc='best')

ax_bar.set_yticks(range(len(plant_samples)))
ax_bar.set_yticklabels(plant_samples, fontsize=14)
ax_bar.set_ylim(-padding, plant_samples_count)

ax_bar.set_yticks(range(plant_samples_count))
ax_bar.set_yticklabels(plant_samples, fontsize=12)
ax_bar.set_xlim([0, 6e9])

ax_bar.xaxis.set_major_locator(MultipleLocator(1e9))
ax_bar.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x/1e9:.0f}E9'))

ax_bar.tick_params(axis='both', which='major', labelsize=14)
ax_bar.set_xlim(0, total_peak_area.max() * 1.1)

plt.tight_layout()

output_image_path = 'C:/Users/zhang/Desktop/Figure6adef-2.png'
plt.savefig(output_image_path, dpi=600, bbox_inches='tight')
plt.show()