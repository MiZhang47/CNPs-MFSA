import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = 'Arial'

categories = ['SIRIUS', 'CSNPs-MFSA']
methods = ['Annotated', 'Unannotated', 'Correct (Top1)', 'Error (Top 1']
data = {
    'SIRIUS': [98.2, 1.8, 25.9, 74.1],
    'CSNPs-MFSA': [96.5, 3.5, 74.1, 25.9]
}

bar_width = 0.3
method_gap = 0.2
x = np.arange(len(categories)) * 1.0

fig, ax = plt.subplots(figsize=(4, 4))

for i, method in enumerate(methods):
    if i == 0:
        ax.bar(x - method_gap, [data[cat][i] for cat in categories], width=bar_width, label=method, color='#66c2a5')
    elif i == 1:
        ax.bar(x - method_gap, [data[cat][i] for cat in categories],
               bottom=[data[cat][0] for cat in categories], width=bar_width, label=method, color='#fc8d62')
    elif i == 2:
        ax.bar(x + method_gap, [data[cat][i] for cat in categories], width=bar_width, label=method, color='#8da0cb')
    elif i == 3:
        ax.bar(x + method_gap, [data[cat][i] for cat in categories],
               bottom=[data[cat][2] for cat in categories], width=bar_width, label=method, color='grey')

for i, cat in enumerate(categories):
    ax.text(x[i] - method_gap, data[cat][0] / 2, f"{data[cat][0]}%", ha='center', va='center', color='white', fontsize=10)
    ax.text(x[i] - method_gap, data[cat][0] + data[cat][1] / 2, f"{data[cat][1]}%", ha='center', va='center', color='black', fontsize=10)
    ax.text(x[i] + method_gap, data[cat][2] / 2, f"{data[cat][2]}%", ha='center', va='center', color='black', fontsize=10)
    ax.text(x[i] + method_gap, data[cat][2] + data[cat][3] / 2, f"{data[cat][3]}%", ha='center', va='center', color='white', fontsize=10)

ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=2, fontsize=10, frameon=False)

ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=12)
ax.set_ylabel('Percentage (%)', fontsize=12)
ax.set_ylim(0, 110)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()

output_path = 'C:/Users/zhang/Desktop/Figure4a_Top1 Accuracy.png'
plt.savefig(output_path, dpi=600)
plt.close()