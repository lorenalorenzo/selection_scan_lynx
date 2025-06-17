import pandas as pd
import glob
import os

# Pattern to find the files
file_pattern = os.path.join("files", "*_repetitive_genomic_regions_annotated")

# Get the list of matching files
file_list = glob.glob(file_pattern)

print("Found files:", file_list)

if not file_list:
    print("⚠️ No matching files found. Check the directory and file names.")
else:
    # Read using space as separator
    df_all = pd.concat(
        (pd.read_csv(f, delim_whitespace=True) for f in file_list),
        ignore_index=True
    )

    # Export to Excel
    df_all.to_excel("merged_repetitive_genomic_regions.xlsx", index=False)
    print("✅ Merged and saved to Excel.")