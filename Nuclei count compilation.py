# Load packages
# Python 3.12.10
import os
import pandas as pd  # Version 2.2.3


# Set working directory
os.chdir(<USER/Target_Folder>)

# Load nuclei measurements of all slides as dataframes
nuclei_f1S1 = pd.read_csv("First round/Slide 1/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f1S2 = pd.read_csv("First round/Slide 2/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f1S3 = pd.read_csv("First round/Slide 3/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f1S4 = pd.read_csv("First round/Slide 4/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f1S5 = pd.read_csv("First round/Slide 5/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f1S6 = pd.read_csv("First round/Slide 6/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f1S7 = pd.read_csv("First round/Slide 7/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f1S8 = pd.read_csv("First round/Slide 8/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f1S9 = pd.read_csv("First round/Slide 9/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f1S10 = pd.read_csv("First round/Slide 10/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f1S11 = pd.read_csv("First round/Slide 11/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f1S12 = pd.read_csv("First round/Slide 12/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f2S1 = pd.read_csv("Second round/Slide 1/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f2S2 = pd.read_csv("Second round/Slide 2/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f2S3 = pd.read_csv("Second round/Slide 3/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f2S4 = pd.read_csv("Second round/Slide 4/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f2S5 = pd.read_csv("Second round/Slide 5/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f2S6 = pd.read_csv("Second round/Slide 6/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f2S7 = pd.read_csv("Second round/Slide 7/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f2S8 = pd.read_csv("Second round/Slide 8/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f2S9 = pd.read_csv("Second round/Slide 9/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f2S10 = pd.read_csv("Second round/Slide 10/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f2S11 = pd.read_csv("Second round/Slide 11/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")
nuclei_f2S12 = pd.read_csv("Second round/Slide 12/Tif images/nuclei_segmentation/Total_nuclei_measurements.csv")

# Add all dataframes into a dictinionary
all_slides = { 
    "f1S1": nuclei_f1S1, 
    "f1S2": nuclei_f1S2, 
    "f1S3":nuclei_f1S3, 
    "f1S4":nuclei_f1S4, 
    "f1S5":nuclei_f1S5, 
    "f1S6":nuclei_f1S6, 
    "f1S7":nuclei_f1S7, 
    "f1S8":nuclei_f1S8, 
    "f1S9": nuclei_f1S9, 
    "f1S10":nuclei_f1S10, 
    "f1S11":nuclei_f1S11, 
    "f1S12":nuclei_f1S12, 
    "f2S1": nuclei_f2S1, 
    "f2S2": nuclei_f2S2, 
    "f2S3": nuclei_f2S3, 
    "f2S4": nuclei_f2S4, 
    "f2S5":nuclei_f2S5, 
    "f2S6":nuclei_f2S6,
    "f2S7":nuclei_f2S7, 
    "f2S8":nuclei_f2S8, 
    "f2S9":nuclei_f2S9, 
    "f2S10":nuclei_f2S10, 
    "f2S11":nuclei_f2S11, 
    "f2S12":nuclei_f2S12
}

# Generate empty list to store the results from for loop
all_slides_counts = []

# Run through all unique Image_IDs and according nuclei_counts for each region/square on each slide and combine them in one table
for name, df in all_slides.items():
    unique_combinations = df[["Image_ID", "nuclei_count"]].drop_duplicates()
    unique_combinations.insert(0, "source", name)
    all_slides_counts.append(unique_combinations)

nuclei_all_slides = pd.concat(all_slides_counts, ignore_index = True)    

print(nuclei_all_slides)

# Save the output as an Excel sheet
output_dir = <USER/Target_folder>
output_file = "Nuclei_counts_all_slides.xlsx"

output_path = os.path.join(output_dir, output_file)
nuclei_all_slides.to_excel(output_path, index = False)

file_path = <USER/target_folder>

# Create separate sheets containing the nuclei counts for all regions/squares of one slide
with pd.ExcelWriter(file_path, engine="openpyxl", mode="a") as writer:
    for source in nuclei_all_slides.iloc[:, 0].unique():
        subset = nuclei_all_slides[nuclei_all_slides.iloc[:, 0] == source]
        sheet_name = str(source)
        subset.to_excel(writer, sheet_name=sheet_name, index=False)


