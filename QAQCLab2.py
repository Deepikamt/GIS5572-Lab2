import os
os.environ["CRYPTOGRAPHY_OPENSSL_NO_LEGACY"] = "1"

#Libraries required

import pandas as pd
import arcpy
import json
import requests
import os
import arcgis
import warnings

# to suppress warnings during the execution
warnings.filterwarnings("ignore")

#path to local database
local_gdb = r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\ArcGIS-II-Lab2.gdb"


# Establish SDE Connection via PGAdmin & Catalog Pane in ArcGIS Pro
db = r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\PostgreSQL-35-gis5572(postgres)(1).sde"

db = {
    "server": "35.224.213.125",
    "database": "gis5572",
    "user": "postgres",
    "password": "********",
    "port": "5432"
}

# defining a function to check the data quality
def check_raster(file_path, categorical=True, expected_cell_size=None, expected_srid=None, minnesota_boundary=None):
    """
    A function to check the quality of a raster dataset.
    """
    # Check for Null Values
    null_values = arcpy.management.GetRasterProperties(file_path, "ANYNODATA").getOutput(0)

    if null_values == "1":
        print("Null values exist.")
    else:
        print("Null values do not exist.")

    # Check if Cell Size is Correct
    x_size = float(arcpy.management.GetRasterProperties(file_path, "CELLSIZEX").getOutput(0))
    y_size = float(arcpy.management.GetRasterProperties(file_path, "CELLSIZEY").getOutput(0))

    if x_size == expected_cell_size and y_size == expected_cell_size:
        print("Actual spatial resolution matches expected spatial resolution.")
    else:
        print("Actual spatial resolution does not match expected spatial resolution.")

    # If Dataset is not Categorical, Check if there are Outliers
    if not categorical:
        mean_val = float(arcpy.management.GetRasterProperties(file_path, "MEAN").getOutput(0))
        std_val = float(arcpy.management.GetRasterProperties(file_path, "STD").getOutput(0))

        max_val = float(arcpy.management.GetRasterProperties(file_path, "MAXIMUM").getOutput(0))
        min_val = float(arcpy.management.GetRasterProperties(file_path, "MINIMUM").getOutput(0))

        # Check if Min < Mean - 3 Std Devs or if Max > Mean + 3 Std Devs
        if min_val < (mean_val - (3 * std_val)) or max_val > (mean_val + (3 * std_val)):
            print("Outliers exist within the dataset. Values exist outside of +- 3 standard deviations of the mean.")
        else:
            print("Outliers do not exist within the dataset. No values +- 3 standard deviations of the mean.")
    else:
        print("Raster is categorical. Not checking for outliers.")

    # Check CRS of Raster
    sr = arcpy.Describe(file_path).spatialReference

    if expected_srid is None:
        print(f"Coordinate system of the raster is: {sr}")
    else:
        arcpy_expected_sr = arcpy.SpatialReference(expected_srid)

        if arcpy_expected_sr.factoryCode == sr.factoryCode:
            print("Actual coordinate system matches expected coordinate system.")
        else:
            print("Actual coordinate system does not match expected coordinate system.")
            print(f"Coordinate system of the raster is: {sr.factoryCode}")

    # Check if Raster is within Minnesota Boundary
    if minnesota_boundary:
        left = float(arcpy.management.GetRasterProperties(file_path, "LEFT").getOutput(0))
        bottom = float(arcpy.management.GetRasterProperties(file_path, "BOTTOM").getOutput(0))
        right = float(arcpy.management.GetRasterProperties(file_path, "RIGHT").getOutput(0))
        top = float(arcpy.management.GetRasterProperties(file_path, "TOP").getOutput(0))
        
        # Check if raster extent overlaps with Minnesota boundary
        if right > minnesota_boundary[0] and left < minnesota_boundary[2] and top > minnesota_boundary[1] and bottom < minnesota_boundary[3]:
            print("Raster overlaps with Minnesota boundary.")
        else:
            print("Raster does not overlap with Minnesota boundary.")
    else:
        print("Minnesota boundary not provided. Skipping boundary check.")


# NCLD Land Cover (2019), Minnesota
landcover_path = r"D:\4th sem_minnesota\ArcGIS_II\Lab 2\NLCD_2019_Land_Cover.tif\NLCD_2019_Land_Cover.tif"
mn_boundary_gdb = r"D:\4th sem_minnesota\ArcGIS_II\Lab 2\Projected\MinnesotaBoundry.shp"


# Output name for the clipped raster within the geodatabase
output_name = "NLCD_2019_Land_Cover_clipped"

# Output path for the clipped raster within the geodatabase
output_raster = arcpy.os.path.join(r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\ArcGIS-II-Lab2.gdb", output_name)

try:
    # Perform the clip operation
    arcpy.management.Clip(landcover_path, mn_boundary_gdb, output_raster)
    print("Clipping completed successfully.")
except arcpy.ExecuteError:
    print(arcpy.GetMessages(2))

# Set the paths to the input raster and the geodatabase
mnlandcover_path = r"D:\4th sem_minnesota\ArcGIS_II\Lab 2\NLCD_2019_Land_Cover.tif\NLCD_2019_Land_Cover.tif"
output_gdb = r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\ArcGIS-II-Lab2.gdb"   
    
# Extract the file name from the input raster path
raster_name = arcpy.Describe(mnlandcover_path).baseName

# Sanitize the raster name to remove invalid characters
raster_name = arcpy.ValidateTableName(raster_name, output_gdb)

# Build the full output path within the geodatabase
output_raster = arcpy.os.path.join(output_gdb, raster_name)

# Use the Copy Raster tool to copy the raster to the geodatabase
arcpy.CopyRaster_management(dem_path, output_raster)

print(f"Raster '{raster_name}' copied to geodatabase.")

def check_raster(file_path, categorical=True, expected_cell_size=None, expected_srid=None, xmin=None, ymin=None, xmax=None, ymax=None):
    """
    Check the quality of a raster dataset.
    """
    # Null Values Check
    null_values = arcpy.management.GetRasterProperties(file_path, "ANYNODATA").getOutput(0)
    null_msg = "Detected: Null values present. Initiating data quality check." if null_values == "1" else "No null values found."
    print(null_msg)

    # Spatial Resolution Check
    x_size = float(arcpy.management.GetRasterProperties(file_path, "CELLSIZEX").getOutput(0))
    y_size = float(arcpy.management.GetRasterProperties(file_path, "CELLSIZEY").getOutput(0))
    resolution_msg = "Spatial resolution matches expected specs." if x_size == expected_cell_size and y_size == expected_cell_size else "Spatial resolution does not meet expected specifications."
    print(resolution_msg)

    # Outliers Check
    if not categorical:
        mean_val = float(arcpy.management.GetRasterProperties(file_path, "MEAN").getOutput(0))
        std_val = float(arcpy.management.GetRasterProperties(file_path, "STD").getOutput(0))
        max_val = float(arcpy.management.GetRasterProperties(file_path, "MAXIMUM").getOutput(0))
        min_val = float(arcpy.management.GetRasterProperties(file_path, "MINIMUM").getOutput(0))
        outliers_msg = "Warning: Outliers detected. Verify data integrity." if min_val < (mean_val - (3 * std_val)) or max_val > (mean_val + (3 * std_val)) else "No outliers found. Data integrity confirmed."
        print(outliers_msg)
    else:
        print("Categorical raster. Outliers not applicable.")

    # CRS Check
    sr = arcpy.Describe(file_path).spatialReference
    crs_msg = f"Coordinate system matches the expected specs: {sr}" if expected_srid is None else "Coordinate system matches expected specification." if arcpy.SpatialReference(expected_srid).factoryCode == sr.factoryCode else "Coordinate system does not match expected specification."
    print(crs_msg)

    # Bounding Box Check
    if None not in [xmin, ymin, xmax, ymax]:
        left = float(arcpy.management.GetRasterProperties(file_path, "LEFT").getOutput(0))
        bottom = float(arcpy.management.GetRasterProperties(file_path, "BOTTOM").getOutput(0))
        right = float(arcpy.management.GetRasterProperties(file_path, "RIGHT").getOutput(0))
        top = float(arcpy.management.GetRasterProperties(file_path, "TOP").getOutput(0))
        bbox_msg = "Caution: Raster extends slightly beyond specified boundary. May require closer inspection." if left < xmin or bottom < ymin or right > xmax or top > ymax else "Raster fits entirely within specified boundary." if left >= xmin and bottom >= ymin and right <= xmax and top <= ymax else "Raster partially outside specified boundary."
        print(bbox_msg)
    else:
        print("Bounding box check skipped.")

# Path to Land Cover Dataset
geodatabase_path = r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\ArcGIS-II-Lab2.gdb"
NLCD_raster = "NLCD_2019_Land_Cover_clipped"
NLCD_raster_path = arcpy.os.path.join(geodatabase_path, NLCD_raster)

# Define Minnesota boundary extent
minnesota_boundary = [-97.5, 43.0, -89.00, 49.5]  # xmin, ymin, xmax, ymax


#To check the Land Cover Raster
check_raster(NLCD_raster_path, categorical=True, expected_cell_size=30, expected_srid=26915, xmin=-97.5, ymin=43.0, xmax=-89.00, ymax=49.5)

# Step 1: Downsample Raster
def downsample_raster(input_raster, output_raster, cell_size):
    arcpy.Resample_management(input_raster, output_raster, cell_size)

# Step 2: Convert Raster to Points
def raster_to_points(input_raster, output_points):
    arcpy.RasterToPoint_conversion(input_raster, output_points, "VALUE")

# Step 3: Upload Points to SDE
def upload_points_to_sde(input_points, output_sde_connection, output_sde_feature_class):
    arcpy.FeatureClassToFeatureClass_conversion(input_points, output_sde_connection, output_sde_feature_class)

# Paths and parameters
input_raster = r'C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\ArcGIS-II-Lab2.gdb\NLCD_2019_Land_Cover_clipped'
output_downsampled_raster = r"D:\4th sem_minnesota\ArcGIS_II\Lab 2\downsampled_nlcdraster.tif"
output_points = r'D:\4th sem_minnesota\ArcGIS_II\Lab 2\nlcd_points.shp'
output_sde_connection = r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\PostgreSQL-35-gis5572(postgres)(1).sde"

output_sde_feature_class = 'nlcd_points_in_sde'

# Step 1: Downsample Raster
downsample_raster(input_raster, output_downsampled_raster, "15000")

# Step 2: Convert Raster to Points
raster_to_points(output_downsampled_raster, output_points)

# Step 3: Upload Points to SDE
upload_points_to_sde(output_points, output_sde_connection, output_sde_feature_class)


# Set the paths to the input raster and the geodatabase
dem_path = r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\mn_dem\elev_30m_digital_elevation_model.gdb\digital_elevation_model_30m"
output_gdb = r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\ArcGIS-II-Lab2.gdb"

# Extract the file name from the input raster path
raster_name = arcpy.Describe(dem_path).baseName

# Sanitize the raster name to remove invalid characters
raster_name = arcpy.ValidateTableName(raster_name, output_gdb)

# Build the full output path within the geodatabase
output_raster = arcpy.os.path.join(output_gdb, raster_name)

# Use the Copy Raster tool to copy the raster to the geodatabase
arcpy.CopyRaster_management(dem_path, output_raster)

print(f"Raster '{raster_name}' copied to geodatabase.")

# Path to DEM Dataset
geodatabase_path = r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\ArcGIS-II-Lab2.gdb"
DEM_raster = "digital_elevation_model_30m"
DEM_raster_path = arcpy.os.path.join(geodatabase_path, DEM_raster)

#To check the DEM Raster
check_raster(DEM_raster_path, False, 30, 26915, -97.5, 43.0, -89.00, 49.5)

# Step 1: Downsample Raster
def downsample_raster(input_raster, output_raster, cell_size):
    arcpy.Resample_management(input_raster, output_raster, cell_size)

# Step 2: Convert Raster to Points
def raster_to_points(input_raster, output_points):
    arcpy.RasterToPoint_conversion(input_raster, output_points, "VALUE")

# Step 3: Upload Points to SDE
def upload_points_to_sde(input_points, output_sde_connection, output_sde_feature_class):
    arcpy.FeatureClassToFeatureClass_conversion(input_points, output_sde_connection, output_sde_feature_class)

# Paths and parameters
input_raster = r'C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\ArcGIS-II-Lab2.gdb\digital_elevation_model_30m'
output_downsampled_raster = r"D:\4th sem_minnesota\ArcGIS_II\Lab 2\downsampled_demraster.tif"
output_points = r'D:\4th sem_minnesota\ArcGIS_II\Lab 2\dem_points.shp'
output_sde_connection = r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\PostgreSQL-35-gis5572(postgres)(1).sde"

output_sde_feature_class = 'dem_points_in_sde'

# Step 1: Downsample Raster
downsample_raster(input_raster, output_downsampled_raster, "15000")

# Step 2: Convert Raster to Points
raster_to_points(output_downsampled_raster, output_points)

# Step 3: Upload Points to SDE
upload_points_to_sde(output_points, output_sde_connection, output_sde_feature_class)


weather_url = r"https://mesonet.agron.iastate.edu/api/1/daily.geojson?network=MN_RWIS&month=2&year=2024"

# Send a GET request to the URL
response = requests.get(weather_url)

# Specify a different file path
output_file_path = r'D:\4th sem_minnesota\ArcGIS_II\Lab 2\weather_data_feb24.geojson'

# Save the GeoJSON data to the new file path
with open(output_file_path, 'w') as file:
    json.dump(weather_data, file)

print(f"GeoJSON data has been successfully saved as '{output_file_path}'")


# Define the URL for retrieving GeoJSON weather data
weather_url = r"https://mesonet.agron.iastate.edu/api/1/daily.geojson?network=MN_RWIS&month=2&year=2024"

# Fetch weather data and create a DataFrame
weather_response = requests.get(weather_url)
weather_json = weather_response.json()["features"]
weather_df_raw = pd.DataFrame.from_records(weather_json)

# Function to extract properties from dictionaries as DataFrame columns
def extractToCol(field):
    weather_df_raw[field] = weather_df_raw["properties"].apply(lambda x: dict(x).get(field))

# Define the list of properties to extract
weather_props = ["date", "station", "name", "min_tmpf", "max_tmpf", "precip"]

# Extract properties from dictionaries and create DataFrame columns
for i in weather_props:
    extractToCol(i)


# Extract geometry information (longitude and latitude) and create new columns in the DataFrame.
weather_df_raw["x"] = weather_df_raw["geometry"].apply(lambda x: dict(x)["coordinates"][0])
weather_df_raw["y"] = weather_df_raw["geometry"].apply(lambda x: dict(x)["coordinates"][1])


# Copying Relevant Columns to a New DataFrame:
weather_df = weather_df_raw[["date", "station", "name", "min_tmpf", "max_tmpf", "precip", "x", "y"]].copy()

# Display summary statistics of the original DataFrame
print("Original DataFrame Summary:")
print(weather_df.describe())

# Fill NaN values in 'precip' with 0
weather_df["precip"].fillna(0, inplace=True)

# Display summary statistics after filling null values
print("\nDataFrame After Filling Nulls:")
print(weather_df.describe())

# Drop rows with null latitude or longitude values
weather_df = weather_df.dropna(subset=["x", "y"])

# Display summary statistics after dropping null values
print("\nDataFrame After Dropping Nulls:")
print(weather_df.describe())

# Convert data types
# "station" to string, "name" to string, and "date" to datetime.
weather_df["station"] = weather_df["station"].astype(str)
weather_df["name"] = weather_df["name"].astype(str)
weather_df["date"] = weather_df["date"].astype('datetime64[ns]')


# Display summary statistics after data type conversion
print("\nDataFrame After Data Type Conversion:")
print(weather_df.describe())

# Drop rows where 'precip' is less than 0
weather_df = weather_df.loc[weather_df["precip"] >= 0]

# Display summary statistics after dropping rows where 'precip' is less than 0
print("\nDataFrame After Dropping Rows with 'precip' < 0:")
print(weather_df.describe())

# Drop outliers for 'max_tmpf'
# Assuming temp. values follow a perfect normal distribution
# approx. 99.7% of the data falls within 3 SDs of the mean
mxtmp_mn = weather_df["max_tmpf"].mean()
mxtmp_std = weather_df["max_tmpf"].std()
weather_df = weather_df.loc[weather_df["max_tmpf"].between(mxtmp_mn - 3 * mxtmp_std, mxtmp_mn + 3 * mxtmp_std)]
# basically removes outliers in the "max_tmpf" column using a threshold of 3 SDs from the mean


# Display summary statistics after dropping max_temp outliers
print("\nDataFrame After Dropping Max Temp Outliers:")
print(weather_df.describe())

# Display summary statistics after commenting out min_temp and precip outlier removal
print("\nDataFrame After Commenting Out Min Temp and Precip Outliers:")
print(weather_df.describe())

# Comment out bounding box filtering
# weather_df = weather_df.loc[(weather_df["x"] > -97.5) & (weather_df["x"] < -89.0) & (weather_df["y"] > 43.0) & (weather_df["y"] < 49.5)]

# Display summary statistics after commenting out bounding box filtering
print("\nDataFrame After Commenting Out Bounding Box Filtering:")
print(weather_df.describe())

# Display the final DataFrame
print("\nFinal DataFrame:")
print(weather_df)


# Specify the complete file path
save_path = r'D:\4th sem_minnesota\ArcGIS_II\Lab 2\mn_temperature\weather_data_feb24.csv'

# Save the final DataFrame to a CSV file at the specified location
weather_df.to_csv(save_path, index=False)

print(f"Final DataFrame has been successfully saved as '{save_path}'")


import arcpy

# Function to upload CSV to SDE as feature class
def upload_csv_to_sde(csv_path, sde_connection, output_feature_class):
    try:
        # Convert CSV to table
        arcpy.TableToTable_conversion(csv_path, arcpy.env.workspace, "temp_table")
        
        # Create XY event layer from table
        arcpy.MakeXYEventLayer_management("temp_table", "x", "y", "temp_event_layer")
        
        # Copy features to output feature class
        arcpy.FeatureClassToFeatureClass_conversion("temp_event_layer", sde_connection, output_feature_class)
        
        print(f"Feature class '{output_feature_class}' uploaded successfully to SDE.")
    except Exception as e:
        print(f"Error uploading feature class to SDE: {str(e)}")

# Paths and parameters
csv_path = r'D:\4th sem_minnesota\ArcGIS_II\Lab 2\mn_temperature\weather_data_feb24.csv'
sde_connection = r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\PostgreSQL-35-gis5572(postgres)(1).sde"
output_feature_class = 'weather_sde'

# Upload CSV to SDE as feature class
upload_csv_to_sde(csv_path, sde_connection, output_feature_class)


import os
import requests

def data_download(data_url, target_folder_path, file_name):
    # Create the target folder if it does not exist
    if not os.path.exists(target_folder_path):
        os.makedirs(target_folder_path)
    
    # Download the data
    response = requests.get(data_url)
    
    # Construct the full path for the target file
    target_file_path = os.path.join(target_folder_path, file_name)

    # Save the data to file
    with open(target_file_path, 'wb') as file:
        file.write(response.content)

# Define the URL for the corn data
corn_data_url = r"https://ndawn.ndsu.nodak.edu/table.csv?ttype=cogdd&station=78&station=111&station=98&station=162&station=174&station=142&station=164&station=138&station=161&station=9&station=160&station=159&station=10&station=118&station=56&station=165&station=11&station=12&station=58&station=13&station=84&station=55&station=179&station=7&station=186&station=87&station=14&station=15&station=96&station=191&station=16&station=201&station=137&station=124&station=143&station=17&station=85&station=140&station=134&station=18&station=136&station=65&station=104&station=99&station=192&station=19&station=129&station=20&station=101&station=166&station=178&station=81&station=21&station=97&station=22&station=75&station=184&station=2&station=172&station=139&station=158&station=23&station=157&station=62&station=86&station=24&station=89&station=126&station=167&station=93&station=183&station=90&station=25&station=83&station=107&station=156&station=77&station=26&station=155&station=70&station=127&station=144&station=27&station=173&station=132&station=28&station=195&station=185&station=29&station=30&station=154&station=31&station=187&station=102&station=32&station=119&station=4&station=80&station=33&station=59&station=153&station=105&station=82&station=34&station=198&station=72&station=135&station=35&station=76&station=120&station=141&station=109&station=36&station=79&station=193&station=71&station=37&station=38&station=189&station=39&station=130&station=73&station=188&station=40&station=41&station=54&station=69&station=194&station=145&station=113&station=128&station=42&station=43&station=103&station=171&station=116&station=196&station=88&station=114&station=3&station=163&station=200&station=64&station=115&station=168&station=67&station=175&station=146&station=170&station=197&station=44&station=133&station=106&station=100&station=121&station=45&station=46&station=61&station=66&station=181&station=74&station=60&station=199&station=125&station=176&station=177&station=8&station=180&station=204&station=47&station=122&station=108&station=5&station=152&station=48&station=151&station=147&station=68&station=169&station=49&station=50&station=91&station=182&station=117&station=63&station=150&station=51&station=6&station=52&station=92&station=112&station=131&station=123&station=95&station=53&station=203&station=190&station=57&station=149&station=148&station=202&station=110&year=2024&begin_date=2023-05-01&end_date=2023-08-31"

# Specify the target folder path
target_folder_path = r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2"

# Download the corn data CSV file
data_download(corn_data_url, target_folder_path, "corn_data.csv")



import os
import pandas as pd
import arcpy
from arcgis import GeoAccessor

# Step 1: Download Corn Data
ndawn_url = r"https://ndawn.ndsu.nodak.edu/table.csv?ttype=cogdd&station=78&station=111&station=98&station=162&station=174&station=142&station=164&station=138&station=161&station=9&station=160&station=159&station=10&station=118&station=56&station=165&station=11&station=12&station=58&station=13&station=84&station=55&station=179&station=7&station=186&station=87&station=14&station=15&station=96&station=191&station=16&station=201&station=137&station=124&station=143&station=17&station=85&station=140&station=134&station=18&station=136&station=65&station=104&station=99&station=192&station=19&station=129&station=20&station=101&station=166&station=178&station=81&station=21&station=97&station=22&station=75&station=184&station=2&station=172&station=139&station=158&station=23&station=157&station=62&station=86&station=24&station=89&station=126&station=167&station=93&station=183&station=90&station=25&station=83&station=107&station=156&station=77&station=26&station=155&station=70&station=127&station=144&station=27&station=173&station=132&station=28&station=195&station=185&station=29&station=30&station=154&station=31&station=187&station=102&station=32&station=119&station=4&station=80&station=33&station=59&station=153&station=105&station=82&station=34&station=198&station=72&station=135&station=35&station=76&station=120&station=141&station=109&station=36&station=79&station=193&station=71&station=37&station=38&station=189&station=39&station=130&station=73&station=188&station=40&station=41&station=54&station=69&station=194&station=145&station=113&station=128&station=42&station=43&station=103&station=171&station=116&station=196&station=88&station=114&station=3&station=163&station=200&station=64&station=115&station=168&station=67&station=175&station=146&station=170&station=197&station=44&station=133&station=106&station=100&station=121&station=45&station=46&station=61&station=66&station=181&station=74&station=60&station=199&station=125&station=176&station=177&station=8&station=180&station=204&station=47&station=122&station=108&station=5&station=152&station=48&station=151&station=147&station=68&station=169&station=49&station=50&station=91&station=182&station=117&station=63&station=150&station=51&station=6&station=52&station=92&station=112&station=131&station=123&station=95&station=53&station=203&station=190&station=57&station=149&station=148&station=202&station=110&year=2024&begin_date=2023-05-01&end_date=2023-08-31"
data_download(ndawn_url, target_folder_path=r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\ArcGIS-II-Lab2.gdb", file_name="corn_data.csv")

# Step 2: Read and Filter the Data
corn_data_path = "corn_data.csv"
corn_df = pd.read_csv(corn_data_path, skiprows=[0, 1, 2, 4])
columns_to_print = ['Station Name', 'Latitude', 'Longitude', 'Elevation', 'Year', 'Month', 'Day', 'Max Temp', 'Min Temp', 'Corn Accumulated Growing Degree Days']

# Define geographical boundaries
north = 45.0
south = 44.0
east = -92.0
west = -94.0

minnesota_corn_df = corn_df[(corn_df['Latitude'] <= north) & (corn_df['Latitude'] >= south) & (corn_df['Longitude'] <= east) & (corn_df['Longitude'] >= west)]
minnesota_corn_df = minnesota_corn_df[columns_to_print]

# Step 3: Convert Data to Spatial DataFrame
corn_spatial_df = GeoAccessor.from_xy(minnesota_corn_df, "Longitude", "Latitude")

# Step 4: Upload Data to Local Geodatabase
local_gdb = r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\ArcGIS-II-Lab2.gdb"
corn_feature_class_name = "Corn_AGDD"
corn_spatial_df.spatial.to_featureclass(location=os.path.join(local_gdb, corn_feature_class_name))

print(f"Feature class '{corn_feature_class_name}' uploaded successfully to local geodatabase.")

# Step 5: Upload Data to PostgreSQL
output_sde_connection = r"C:\Users\Deepika\OneDrive\Documents\ArcGIS\Projects\ArcGIS-II-Lab2\PostgreSQL-35-gis5572(postgres)(1).sde"
output_sde_feature_class = 'corn_agdd_in_sde'
arcpy.FeatureClassToFeatureClass_conversion(corn_feature_class_name, output_sde_connection, output_sde_feature_class)
print(f"Feature class '{output_sde_feature_class}' uploaded successfully to PostgreSQL.")



