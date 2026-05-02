# Load packages
# Python 3.12.10
import os
import glob
import cv2  # Version 4.11.0
import numpy as np  # Version 2.2.5
import pandas as pd  # Version 2.2.3
from skimage import measure  # Version 0.25.2
from skimage.measure import label  # Version 0.25.2

# Convert pixel size to microns according to metadata of image
pixels_to_um = 0.3442348

# Find all TIF files in working directory 
imaging_files = glob.glob("*.tif")

# Store all DataFrames and merge in the end for nuclei countings
all_dfs = []

# Batch Processing
for file_name in imaging_files:
    print(f"Processing: {file_name}")

    # Load imaging data and extract blue channel for nuclei
    img = cv2.imread(file_name)
    nuclei=img[:,:,0]

    # Set threshold to differentiate particles from background
    ret1, thresh = cv2.threshold(nuclei, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU) 

    # Set kernel and iterations to filter particles of certain size and intensity and get rid of debris particles
    kernel = np.ones((3,3), np.uint8)
    opening = cv2.morphologyEx(thresh, cv2.MORPH_OPEN, kernel, iterations = 2)

    # Define sure background of picture by dilating nuclei size 
    sure_bg = cv2.dilate(opening, kernel, iterations = 10)

    # Create distance map of center points of  nuclei
    dist_transform = cv2.distanceTransform(opening, cv2.DIST_L2, 5)

    # Define sure foreground (nuclei) by setting cutoff distance of maximum "density" of nuclei
    print(dist_transform.max())
    ret2, sure_fg = cv2.threshold(dist_transform, 0.25*dist_transform.max(), 255,0)

    sure_fg = np.uint8(sure_fg)
    unknown = cv2.subtract(sure_bg,sure_fg)

    # Adding markers and setting apart background 
    ret3, markers = cv2.connectedComponents(sure_fg)
    markers = markers + 10
    markers[unknown == 255] = 0

    # Create watershed mask
    markers = cv2.watershed(img, markers)
    img[markers == -1] = [0, 255, 255]
    
    # Relabel markers, so that counts start from 1
    markers[markers == -1] = 0
    markers = label(markers)

    # Measure properties of markers and nuclei in table (label, area, diameter, mean intensity, solidity, orientation and perimeter)
    props = measure.regionprops_table(markers, nuclei, 
                                  properties = ['label', 'area', 'equivalent_diameter', 'mean_intensity', 
                                              'solidity', 'orientation', 'perimeter', 'centroid'])
    df = pd.DataFrame(props)

    # Add nucleus area, diameter and x and y coordinates of each nucleus in microns to property dataframe
    df['area_sq_microns'] = df['area']* (pixels_to_um**2)
    df['equivalent_diameter_microns'] = df['equivalent_diameter']* (pixels_to_um**2)
    df["centroid_y_um"] = df["centroid-0"] * pixels_to_um
    df["centroid_x_um"] = df["centroid-1"] * pixels_to_um

    df['Image_ID'] = file_name
    df['nuclei_count'] = len(df)

    all_dfs.append(df)

    # Generate labeled nuclei outlines for each image
    props1 = measure.regionprops(markers)
    num_outline_img = img.copy()  

    for region in props1:
        y, x = region.centroid
        label_text = str(region.label)
        cv2.putText(num_outline_img, label_text, (int(x), int(y)),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.4, (255, 255, 255), 1, cv2.LINE_AA)

    # Save labeled pictures into subfolder
    output_dir = os.path.join(os.getcwd(), "nuclei_segmentation")
    os.makedirs(output_dir, exist_ok = True)
    
    base_name = os.path.splitext(file_name)[0]
    output_image_path = os.path.join(output_dir, f"{base_name}_labeled.tif")
    cv2.imwrite(output_image_path, num_outline_img)
    
    print(f" -> Found {len(df)} nuclei in {base_name}")

#Combine list of dataframes into one overview-dataframe and save as .csv and .xlsx
final_df = pd.concat(all_dfs, ignore_index = True)

csv_path = os.path.join(output_dir, "Total_nuclei_measurements.csv")
xlsx_path = os.path.join(output_dir, "Total_nuclei_measurements.xlsx")
final_df.to_csv(csv_path, index = False)
final_df.to_excel(xlsx_path, index = False)

# Summary table
summary = final_df.groupby("Image_ID")["nuclei_count"].first().reset_index()
print("\nSummary of nuclei counts:")
print(summary)
