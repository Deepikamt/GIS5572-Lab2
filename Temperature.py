import os
os.environ["CRYPTOGRAPHY_OPENSSL_NO_LEGACY"] = "1"

#importing modules
import arcgis
import psycopg2
from psycopg2 import sql
import json
import requests
import os
import warnings
from IPython.display import Image, display
import random
import zipfile
import io
from datetime import date


file_path = os.path.dirname(arcpy.mp.ArcGISProject('CURRENT').filePath)

# Change the current working directory to the extracted directory path
os.chdir(file_path)

# Set the workspace environment to the extracted directory path
arcpy.env.workspace = file_path

# Define a spatial reference with ID 26915 (likely UTM zone 15N)
spatial_ref = arcpy.SpatialReference(26915)

# Establish variables for project and map
project = arcpy.mp.ArcGISProject("CURRENT")
m = project.listMaps("Map")[0]

import arcpy
import os

# Define the SDE feature class path
sde_feature_class = "C:\\Users\\Deepika\\OneDrive\\Documents\\ArcGIS\\Projects\\ArcGIS-II-Lab2\\PostgreSQL-35-gis5572(postgres)(1).sde\\gis5572.postgres.weather_sde"

# Define the local geodatabase path
local_gdb = "C:\\Users\\Deepika\\OneDrive\\Documents\\ArcGIS\\Projects\\ArcGIS-II-Lab2\\ArcGIS-II-Lab2.gdb"

# Define the local feature class path
local_feature_class = os.path.join(local_gdb, "mn_clean_weather")

# Process: Copy Features
arcpy.CopyFeatures_management(sde_feature_class, local_feature_class)


import arcpy
import os
import pandas as pd
import random
from arcgis.features import GeoAccessor

def create_sample(input_feature_class, x_field, y_field, z_field, percent, gdb_path, out_feature_class):
    try:
        print('Creating feature class "' + out_feature_class + '"...')
        
        # Check if input feature class exists
        if not arcpy.Exists(input_feature_class):
            raise RuntimeError("Input feature class does not exist or is not accessible.")

        all_rows = []
        with arcpy.da.SearchCursor(in_table=input_feature_class, field_names=[z_field, x_field, y_field]) as cursor:
            for row in cursor:
                Z = row[0]
                X = row[1]
                Y = row[2]
                all_rows.append([Z, X, Y])

        total_rows = len(all_rows)
        sample_num = int(total_rows * percent)
        removed_num = total_rows - sample_num

        random.shuffle(all_rows)

        Z_list = []
        X_list = []
        Y_list = []

        for row in all_rows:
            Z_list.append(row[0])
            X_list.append(row[1])
            Y_list.append(row[2])

        df_list = [Z_list, X_list, Y_list]
        for l in df_list:
            del l[-removed_num:]

        rand_dict = {
            'Z': Z_list,
            'X': X_list,
            'Y': Y_list
        }

        random_sample_df = pd.DataFrame(rand_dict)

        # Convert table to sedf
        sedf = GeoAccessor.from_xy(df=random_sample_df, x_column="X", y_column="Y")

        # Convert sedf to feature class
        sedf.spatial.to_featureclass(location=os.path.join(gdb_path, out_feature_class))

        print('Done')
    except Exception as e:
        print("An error occurred:", str(e))

# Define the SDE feature class path
sde_feature_class = "C:\\Users\\Deepika\\OneDrive\\Documents\\ArcGIS\\Projects\\ArcGIS-II-Lab2\\PostgreSQL-35-gis5572(postgres)(1).sde\\gis5572.postgres.weather_sde"

# Define the local geodatabase path
local_gdb = "C:\\Users\\Deepika\\OneDrive\\Documents\\ArcGIS\\Projects\\ArcGIS-II-Lab2\\ArcGIS-II-Lab2.gdb"

# Define the local feature class path
local_feature_class = os.path.join(local_gdb, "mn_clean_weather")

# Process: Copy Features
arcpy.CopyFeatures_management(sde_feature_class, local_feature_class)

# Create a feature class with 20% of random points sampled
create_sample(
    input_feature_class=local_feature_class,
    x_field='x',
    y_field='y',
    z_field='max_tmpf',
    percent=0.2,
    gdb_path=local_gdb,
    out_feature_class='Random_Sample_20_temp',
)


# Define the well-known ID (WKID) of the spatial reference (e.g., EPSG code)
spatial_ref_wkid = 4326  

import arcpy
import os
arcpy.ga.ExploratoryInterpolation(
    in_features="Random_Sample_20_temp",
    value_field="z",
    out_cv_table=r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\ArcGIS-II-Lab2.gdb\ExploratoryInterpolation2",
    out_geostat_layer=None,
    interp_methods="ORDINARY_KRIGING;UNIVERSAL_KRIGING;EBK",
    comparison_method="SINGLE",
    criterion="ACCURACY",
    criteria_hierarchy="ACCURACY PERCENT #",
    weighted_criteria="ACCURACY 1",
    exclusion_criteria=None
)


with arcpy.EnvManager(scratchWorkspace=r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\ArcGIS-II-Lab2.gdb"):
    out_raster = arcpy.sa.Idw(
        in_point_features="gis5572.postgres.weather_sde",
        z_field="min_tmpf",
        cell_size=0.0218500823959998,
        power=2,
        search_radius="VARIABLE 12",
        in_barrier_polyline_features=None
    )
    out_raster.save(r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\ArcGIS-II-Lab2.gdb\Idw_weather_2")

# KRINGING - UNIVERSAL
with arcpy.EnvManager(mask="Clipped_dem", scratchWorkspace=r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\ArcGIS-II-Lab2.gdb"):
    out_surface_raster = arcpy.sa.Kriging(
        in_point_features="Random_Sample_20_temp",
        z_field="z",
        kriging_model="LinearDrift 0.021850 # # #",
        cell_size=0.0218500824,
        search_radius="VARIABLE 12",
        out_variance_prediction_raster=None
    )
    out_surface_raster.save(r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\ArcGIS-II-Lab2.gdb\Kriging_Universal")

# Empirical Bayesin Kriging

with arcpy.EnvManager(mask="Clipped_dem"):
    arcpy.ga.EmpiricalBayesianKriging(
        in_features="Random_Sample_20_temp",
        z_field="z",
        out_ga_layer=None,
        out_raster=r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\ArcGIS-II-Lab2.gdb\krig_temp",
        cell_size=3000,
        transformation_type="NONE",
        max_local_points=100,
        overlap_factor=1,
        number_semivariograms=100,
        search_neighborhood="NBRTYPE=StandardCircular RADIUS=2.32298810990939 ANGLE=0 NBR_MAX=15 NBR_MIN=10 SECTOR_TYPE=ONE_SECTOR",
        output_type="PREDICTION",
        quantile_value=0.5,
        threshold_type="EXCEED",
        probability_threshold=None,
        semivariogram_model_type="POWER"
    )

# Step 1: Downsample Raster
def downsample_raster(input_raster, output_raster, cell_size):
    arcpy.Resample_management(input_raster, output_raster, cell_size)

# Step 2: Convert Raster to Points
def raster_to_points(input_raster, output_points):
    arcpy.RasterToPoint_conversion(input_raster, output_points, "VALUE")

# Step 3: Upload Points to SDE
def upload_points_to_sde(input_points, output_sde_connection, output_sde_feature_class):
    arcpy.FeatureClassToFeatureClass_conversion(input_points, output_sde_connection, output_sde_feature_class)

# Paths and parameters for each raster
raster_paths = [
    r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\ArcGIS-II-Lab2.gdb\Kriging_Rand1",
    r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\ArcGIS-II-Lab2.gdb\Kriging_Universal",
    r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\ArcGIS-II-Lab2.gdb\krig_temp"
]

output_downsampled_rasters = [
    r'C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\downsampled_EBK.tif',
    r'C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\downsampled_KrigOrd.tif',
    r'C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\downsampled_KrigUni.tif'
]

output_points = [
    r'C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\EBK_Temp_points.shp',
    r'C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\KrigOrd_Temp_points.shp',
    r'C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\KrigUni_Temp_points.shp'
]

output_sde_connection = r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\PostgreSQL-35-gis5572(postgres)(1).sde"

output_sde_feature_classes = [
    'EBK_Temp_points',
    'KrigOrd_Temp_points',
    'KrigUni_Temp_points'
]


import arcpy

# Function to downsample raster
def downsample_raster(input_raster, output_raster, cell_size):
    try:
        arcpy.Resample_management(input_raster, output_raster, cell_size)
        return True
    except arcpy.ExecuteError as e:
        print(f"Error occurred while downsampling raster {input_raster}: {str(e)}")
        return False

# Loop through each raster
for i in range(len(raster_paths)):
    try:
        # Step 1: Downsample Raster
        success = downsample_raster(raster_paths[i], output_downsampled_rasters[i], "15000")
        if not success:
            continue  # Skip to the next raster if downsampling fails
        
        # Step 2: Convert Raster to Points
        arcpy.RasterToPoint_conversion(output_downsampled_rasters[i], output_points[i])
        
        # Step 3: Upload Points to SDE
        upload_points_to_sde(output_points[i], output_sde_connection, output_sde_feature_classes[i])
    except arcpy.ExecuteError as e:
        print(f"Error occurred during processing of raster {raster_paths[i]}: {str(e)}")
    except Exception as e:
        print(f"Unexpected error occurred during processing of raster {raster_paths[i]}: {str(e)}")

print("The entire 3 step process is completed.")


# Loop through each raster
for i in range(len(raster_paths)):
    # Step 1: Downsample Raster
    downsample_raster(raster_paths[i], output_downsampled_rasters[i], "15000")

    # Step 2: Convert Raster to Points
    raster_to_points(output_downsampled_rasters[i], output_points[i])

    # Step 3: Upload Points to SDE
    upload_points_to_sde(output_points[i], output_sde_connection, output_sde_feature_classes[i])

print("The entire 3 step process is completed.")
