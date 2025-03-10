import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = 'Arial'

x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]  # Top k
csnps_1 = [74, 76, 81, 83, 91, 93, 95, 95, 95, 97]
sirius_1 = [26, 34, 38, 40, 43, 43, 43, 43, 43, 47]

csnps_2 = [86, 86, 86, 86, 93, 96, 96, 96, 96, 100]
sirius_2 = [25, 43, 50, 50, 57, 57, 57, 57, 57, 61]

csnps_3 = [63, 67, 77, 80, 90, 93, 93, 93, 93, 97]
sirius_3 = [27, 27, 27, 30, 30, 30, 30, 30, 33, 33]

fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharey=True)

axs[0].step(x, csnps_1, label="CSNPs-MFSA", where='post', linestyle='-', color='red')
axs[0].step(x, sirius_1, label="SIRIUS", where='post', linestyle='-', color='blue')
axs[0].set_title('58 in-house daphnane', fontsize=18)
axs[0].set_xlabel('Top k', fontsize=16)
axs[0].set_ylabel('Accuracy (%)', fontsize=16)
axs[0].tick_params(axis='x', labelsize=14)
axs[0].tick_params(axis='y', labelsize=14)
axs[0].set_xticks(x)
axs[0].legend()
axs[0].grid(True, linestyle='--', alpha=0.7)

axs[1].step(x, csnps_2, label="CSNPs-MFSA", where='post', linestyle='-', color='red')
axs[1].step(x, sirius_2, label="SIRIUS", where='post', linestyle='-', color='blue')
axs[1].set_title('28 normal daphnane', fontsize=18)
axs[1].set_xlabel('Top k', fontsize=16)
axs[1].tick_params(axis='x', labelsize=14)
axs[1].tick_params(axis='y', labelsize=14)
axs[1].set_xticks(x)
axs[1].legend()
axs[1].grid(True, linestyle='--', alpha=0.7)

axs[2].step(x, csnps_3, label="CSNPs-MFSA", where='post', linestyle='-', color='red')
axs[2].step(x, sirius_3, label="SIRIUS", where='post', linestyle='-', color='blue')
axs[2].set_title('30 macrocyclic daphnane', fontsize=18)
axs[2].set_xlabel('Top k', fontsize=16)
axs[2].tick_params(axis='x', labelsize=14)
axs[2].tick_params(axis='y', labelsize=14)
axs[2].set_xticks(x)
axs[2].legend()
axs[2].grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()

output_path = 'C:/Users/zhang/Desktop/Figure4b-TopK accuracy.png'
plt.savefig(output_path, dpi=600)
plt.close()