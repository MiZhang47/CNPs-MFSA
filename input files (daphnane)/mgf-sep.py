import os
import re
import pandas as pd
from pyteomics import mgf

# Paths to your input files
input_path_csv = os.path.join('C:/Users/zhang/Desktop/Thyme56_centroid_rerun_0507/',
                              'Thyme56-centroid_quant-re-mf-ad.csv')
input_path_mgf = os.path.join('C:/Users/zhang/Desktop/Thyme56_centroid_rerun_0507/',
                              'Thyme56-centroid-re.mgf')

# Read the CSV into a DataFrame
df = pd.read_csv(input_path_csv)

# Split the DataFrame by 'type' column
df_D  = df[df['type'] == 'D']
df_MD = df[df['type'] == 'MD']

# Output CSV filenames
output_csv_D  = os.path.join(os.path.dirname(input_path_csv), 'Thyme56-centroid-D_quant-re-mf-ad.csv')
output_csv_MD = os.path.join(os.path.dirname(input_path_csv), 'Thyme56-centroid-MD_quant-re-mf-ad.csv')

# Write out the split CSV files
df_D.to_csv(output_csv_D, index=False)
df_MD.to_csv(output_csv_MD, index=False)

# Collect the row IDs for each group (assuming 'row ID' column holds the matching IDs)
ids_D  = set(df_D['row ID'].astype(str))
ids_MD = set(df_MD['row ID'].astype(str))

# Prepare output MGF filenames
output_mgf_D  = os.path.join(os.path.dirname(input_path_mgf), 'Thyme56-centroid-D-re.mgf')
output_mgf_MD = os.path.join(os.path.dirname(input_path_mgf), 'Thyme56-centroid-MD-re.mgf')

# Open output files for writing
with open(output_mgf_D,  'w') as outD, \
     open(output_mgf_MD, 'w') as outMD:

    # Iterate through each spectrum in the input MGF
    for spectrum in mgf.read(input_path_mgf):
        title = spectrum['params'].get('title', '')
        # Extract the numeric ID from the TITLE string
        match = re.search(r'(\d+)', title)
        if not match:
            continue
        spec_id = match.group(1)

        # Write spectrum to the corresponding MGF
        if spec_id in ids_D:
            mgf.write([spectrum], outD)
        elif spec_id in ids_MD:
            mgf.write([spectrum], outMD)
        # else: ignore spectra whose IDs aren't in the CSV

print("Splitting complete:")
print(f"  CSV -> {output_csv_D}, {output_csv_MD}")
print(f"  MGF -> {output_mgf_D}, {output_mgf_MD}")