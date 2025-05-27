import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter, defaultdict

# 设置样式
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 14

# 修改为你自己的路径
base_dir = r"C:\Users\zhang\Desktop\limonoid\limonoid-MS2-ref"
desktop = os.path.join(os.path.expanduser("~"), "Desktop")

type_map = {
    "TypeI": "type1",
    "TypeII": "type2",
    "TypeIII": "type3",
    # "TypeIV": "type4"
}

mz_bin_size = 0.1
intensity_data = []
presence_data = []
labels = []

# 数据处理
for folder_name, label in type_map.items():
    folder_path = os.path.join(base_dir, folder_name)
    files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]

    for file in files:
        df = pd.read_csv(os.path.join(folder_path, file))
        if 'm/z' not in df.columns or 'Relative Intensity' not in df.columns:
            continue
        df['Relative Intensity'] = df['Relative Intensity'] / df['Relative Intensity'].max()
        df['mz_bin'] = (df['m/z'] / mz_bin_size).round(4) * mz_bin_size
        intensity = df.groupby('mz_bin')['Relative Intensity'].max()
        present = df['mz_bin'].unique()

        intensity_data.append(intensity)
        presence_data.append(set(present))
        labels.append(label)

df_intensity = pd.DataFrame(intensity_data).fillna(0)
df_intensity['type'] = labels
mean_intensity = df_intensity.drop(columns='type').groupby(df_intensity['type']).mean()

# 热图 1 数据
top20_mz_intensity = mean_intensity.mean().sort_values(ascending=False).head(20).index
heatmap1_data = mean_intensity[top20_mz_intensity]

# 重复率计算
type_counts = pd.Series(labels).value_counts().to_dict()
mz_type_counter = defaultdict(Counter)
for present, t in zip(presence_data, labels):
    mz_type_counter[t].update(present)

presence_matrix = pd.DataFrame(mz_type_counter).fillna(0).T
presence_matrix = presence_matrix.div(pd.Series(type_counts), axis=0)

# 热图 2 数据
top20_mz_presence = presence_matrix.mean().sort_values(ascending=False).head(20).index
heatmap2_data = presence_matrix[top20_mz_presence]

# 格式化列名
def format_mz_labels(df):
    df.columns = [f"{float(mz):.4f}" for mz in df.columns]
    return df

heatmap1_data = format_mz_labels(heatmap1_data)
heatmap2_data = format_mz_labels(heatmap2_data)

# 绘图函数并保存
def plot_and_save_heatmap(data, title, filename, cmap="viridis"):
    plt.figure(figsize=(14, 6))
    sns.heatmap(data, annot=True, fmt=".2f", cmap=cmap, cbar_kws={'label': 'Intensity or Frequency'})
    plt.title(title, fontsize=14)
    plt.xlabel("m/z", fontsize=14)
    plt.ylabel("Compound Type", fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    output_path = os.path.join(desktop, filename)
    plt.savefig(output_path, dpi=600)
    plt.close()
    print(f"✅ Saved to {output_path}")

# 绘制并导出
plot_and_save_heatmap(heatmap1_data, "Top 20 Intense m/z Features by Type", "heatmap1_intensity.png")
plot_and_save_heatmap(heatmap2_data, "Top 20 Frequent m/z Features by Type", "heatmap2_presence.png")
