# Load packages
# Python 3.12.10
import os
import pandas as pd  # Version 2.2.3
import numpy as np  # Version 2.2.5

# Set working directory
os.chdir(<USER/target_folder>)

# Load separate sheets with nuclei counts for each slide into dataframes
nuclei_counts_f1S1 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f1S1")
nuclei_counts_f1S2 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f1S2")
nuclei_counts_f1S3 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f1S3")
nuclei_counts_f1S4 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f1S4")
nuclei_counts_f1S5 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f1S5")
nuclei_counts_f1S6 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f1S6")
nuclei_counts_f1S7 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f1S7")
nuclei_counts_f1S8 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f1S8")
nuclei_counts_f1S9 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f1S9")
nuclei_counts_f1S10 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f1S10")
nuclei_counts_f1S11 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f1S11")
nuclei_counts_f1S12 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f1S12")
nuclei_counts_f2S1 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f2S1")
nuclei_counts_f2S2 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f2S2")
nuclei_counts_f2S3 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f2S3")
nuclei_counts_f2S4 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f2S4")
nuclei_counts_f2S5 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f2S5")
nuclei_counts_f2S6 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f2S6")
nuclei_counts_f2S7 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f2S7")
nuclei_counts_f2S8 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f2S8")
nuclei_counts_f2S9 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f2S9")
nuclei_counts_f2S10 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f2S10")
nuclei_counts_f2S11 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f2S11")
nuclei_counts_f2S12 = pd.read_excel("Nuclei_counts_all_slides.xlsx", sheet_name = "f2S12")

# Create dictionary with all dataframes
all_slides = { 
    "f1S1": nuclei_counts_f1S1, 
    "f1S2": nuclei_counts_f1S2, 
    "f1S3":nuclei_counts_f1S3, 
    "f1S4":nuclei_counts_f1S4, 
    "f1S5":nuclei_counts_f1S5, 
    "f1S6":nuclei_counts_f1S6, 
    "f1S7":nuclei_counts_f1S7, 
    "f1S8":nuclei_counts_f1S8, 
    "f1S9": nuclei_counts_f1S9, 
    "f1S10":nuclei_counts_f1S10, 
    "f1S11":nuclei_counts_f1S11, 
    "f1S12":nuclei_counts_f1S12, 
    "f2S1": nuclei_counts_f2S1, 
    "f2S2": nuclei_counts_f2S2, 
    "f2S3": nuclei_counts_f2S3, 
    "f2S4": nuclei_counts_f2S4, 
    "f2S5":nuclei_counts_f2S5, 
    "f2S6":nuclei_counts_f2S6,
    "f2S7":nuclei_counts_f2S7, 
    "f2S8":nuclei_counts_f2S8, 
    "f2S9":nuclei_counts_f2S9, 
    "f2S10":nuclei_counts_f2S10, 
    "f2S11":nuclei_counts_f2S11, 
    "f2S12":nuclei_counts_f2S12
}

# Add column with Slide and Square number for each scan
for sheet_name, df in all_slides.items():
    df["Square"] = df["Image_ID"].str[:9]
    all_slides[sheet_name] = df

# Generate the median, the mean and the standard deviation of all scans grouped by each region/square
stat_df = pd.DataFrame()

all_stats = {}

for sheet, df1 in all_slides.items():

    median = df1.groupby("Square")["nuclei_count"].median().reset_index(name="Median")
    mean   = df1.groupby("Square")["nuclei_count"].mean().reset_index(name="Mean")
    std    = df1.groupby("Square")["nuclei_count"].std().reset_index(name="Std")
    
    stat_df = median.merge(mean, on="Square").merge(std, on="Square")
    
    all_stats[sheet] = stat_df

# Combine output with original sheet name/imaging round and Slide number
all_stats_df = pd.concat(
    [df.assign(Slide=sheet) for sheet, df in all_stats.items()],
    ignore_index=True
)

# Make Slide the first column
cols = ["Slide"] + [c for c in all_stats_df.columns if c != "Slide"]
all_stats_df = all_stats_df[cols]

# Save output as Excel sheet
output_dir = <USER/target_folder>
output_file = "Nuclei_stats.xlsx"

output_path = os.path.join(output_dir, output_file)
all_stats_df.to_excel(output_path, index = False)