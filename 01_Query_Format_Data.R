# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# Query data for ACCS Nelchina project
# Author: Timm Nawrocki, Amanda Droghini, Alaska Center for Conservation Science
# Last Updated: 2025-03-28
# Usage: Script should be executed in R 4.4.3+.
# Description: "Query data for ACCS Nelchina project" pulls data for the ACCS Nelchina project from the AKVEG Database for all available tables.
# ---------------------------------------------------------------------------

# Import required libraries
library(dplyr)
library(fs)
library(janitor)
library(lubridate)
library(readr)
library(readxl)
library(RPostgres)
library(sf)
library(stringr)
library(terra)
library(tibble)
library(tidyr)

#### SET UP DIRECTORIES AND FILES
####------------------------------

# Set root directory (modify to your folder structure)
root_folder = 'ACCS_Work'

# Define input folders
database_repository = path('C:', root_folder, 'Repositories/akveg-database')
credentials_folder = path('C:', root_folder, 'Credentials/akveg_private_read')
project_folder = path('C:', root_folder, 'Projects/VegetationEcology/USNVC_SouthernBoreal_Alpine/Data')
input_folder = path(project_folder, 'Data_Input/region_data')
output_folder = path(project_folder, 'Data_Input/plot_data')

# Define raster files
raster_root = 'D:/ACCS_Work/Data'
elevation_input = path(raster_root, 'topography/Alaska_Composite_DTM_10m/integer/Elevation_10m_3338.tif')
aspect_input = path(raster_root, 'topography/Alaska_Composite_DTM_10m/integer/RadiationAspect_10m_3338.tif')
slope_input = path(raster_root, 'topography/Alaska_Composite_DTM_10m/integer/Slope_10m_3338.tif')
heatload_input = path(raster_root, 'topography/Alaska_Composite_DTM_10m/integer/HeatLoad_10m_3338.tif')
position_input = path(raster_root, 'topography/Alaska_Composite_DTM_10m/integer/Position_10m_3338.tif')
relief_input = path(raster_root, 'topography/Alaska_Composite_DTM_10m/integer/Relief_10m_3338.tif')
rough_input = path(raster_root, 'topography/Alaska_Composite_DTM_10m/integer/Roughness_10m_3338.tif')
wetness_input = path(raster_root, 'hydrography/processed/Wetness_10m_3338.tif')

# Define input files
zone_input = path(input_folder, 'AlaskaYukon_VegetationZones_v1.1_3338.shp') # Use EPSG:3338 for Alaska

# Define output files
taxa_output = path(output_folder, '00_taxonomy.csv')
site_visit_output = path(output_folder, '03_site_visit.csv')
site_points_output = path(output_folder, '03_site_visit_3338.shp') # Spatial output in EPSG:3338 for Alaska
vegetation_output = path(output_folder, '05_vegetation.csv')
abiotic_output = path(output_folder, '06_abiotic_top_cover.csv')
tussock_output = path(output_folder, '07_whole_tussock_cover.csv')
ground_output = path(output_folder, '08_ground_cover.csv')
structural_output = path(output_folder, '09_structural_group_cover.csv')
shrub_output = path(output_folder, '11_shrub_structure.csv')
environment_output = path(output_folder, '12_environment.csv')
soilmetrics_output = path(output_folder, '13_soil_metrics.csv')

# Define queries
taxa_file = path(database_repository, '05_queries/analysis/00_taxonomy.sql')
site_visit_file = path(database_repository, '05_queries/analysis/03_site_visit.sql')
vegetation_file = path(database_repository, '05_queries/analysis/05_vegetation.sql')
abiotic_file = path(database_repository, '05_queries/analysis/06_abiotic_top_cover.sql')
tussock_file = path(database_repository, '05_queries/analysis/07_whole_tussock_cover.sql')
ground_file = path(database_repository, '05_queries/analysis/08_ground_cover.sql')
structural_file = path(database_repository, '05_queries/analysis/09_structural_group_cover.sql')
shrub_file = path(database_repository, '05_queries/analysis/11_shrub_structure.sql')
environment_file = path(database_repository, '05_queries/analysis/12_environment.sql')
soilmetrics_file = path(database_repository, '05_queries/analysis/13_soil_metrics.sql')

# Read local data
zone_shape = st_read(zone_input)
elevation_raster = rast(elevation_input)
aspect_raster = rast(aspect_input)
slope_raster = rast(slope_input)
heatload_raster = rast(heatload_input)
position_raster = rast(position_input)
relief_raster = rast(relief_input)
rough_raster = rast(rough_input)
wetness_raster = rast(wetness_input)

#### QUERY AKVEG DATABASE
####------------------------------

# Import database connection function
connection_script = path(database_repository, 'package_DataProcessing', 'connect_database_postgresql.R')
source(connection_script)

# Create a connection to the AKVEG PostgreSQL database
authentication = path(credentials_folder, 'authentication_akveg_private.csv')
database_connection = connect_database_postgresql(authentication)

# Read taxonomy standard from AKVEG Database
taxa_query = read_file(taxa_file)
taxa_data = as_tibble(dbGetQuery(database_connection, taxa_query))

# Get geometry for intersection
intersect_geometry = st_geometry(zone_shape[zone_shape$zone == 'Boreal Southern' | zone_shape$zone == 'Boreal Central',])

# Read site visit data from AKVEG Database
site_visit_query = read_file(site_visit_file)
site_visit_data = as_tibble(dbGetQuery(database_connection, site_visit_query)) %>%
  # Convert geometries to points with EPSG:4269
  st_as_sf(x = ., coords = c('long_dd', 'lat_dd'), crs = 4269, remove = FALSE) %>%
  # Reproject coordinates to EPSG 3338
  st_transform(crs = st_crs(3338)) %>%
  # Add EPSG:3338 centroid coordinates
  mutate(cent_x = st_coordinates(.$geometry)[,1],
         cent_y = st_coordinates(.$geometry)[,2]) %>%
  # Subset points to those within the target zone
  st_intersection(intersect_geometry) %>% # Intersect with named feature
  # Extract raster data to points
  mutate(elevation_m = terra::extract(elevation_raster, ., raw=TRUE)[,2]) %>%
  mutate(aspect = terra::extract(aspect_raster, ., raw=TRUE)[,2]) %>%
  mutate(slope_deg = terra::extract(slope_raster, ., raw=TRUE)[,2]) %>%
  mutate(heat_load = terra::extract(heatload_raster, ., raw=TRUE)[,2]) %>%
  mutate(position = terra::extract(position_raster, ., raw=TRUE)[,2]) %>%
  mutate(relief = terra::extract(relief_raster, ., raw=TRUE)[,2]) %>%
  mutate(roughness = terra::extract(rough_raster, ., raw=TRUE)[,2]) %>%
  mutate(wetness = terra::extract(wetness_raster, ., raw=TRUE)[,2]) %>%
  # Drop geometry
  st_zm(drop = TRUE, what = "ZM") %>%
  # Filter to Nelchina project
  filter(prjct_cd == 'accs_nelchina_2023') %>%
  # Filter year to include 2022 & 2023
  filter(year(obs_date) >= 2022) %>%
  # Select columns
  dplyr::select(st_vst, prjct_cd, st_code, data_tier, obs_date, scp_vasc, scp_bryo, scp_lich, perspect,
                cvr_mthd, strc_class, elevation_m, aspect, slope_deg, heat_load, position, relief, roughness,
                wetness, homogeneous, plt_dim_m, lat_dd, long_dd, cent_x, cent_y, geometry)

# Export site visit data to shapefile
st_write(site_visit_data, site_points_output, append = FALSE) # Optional to check point selection in a GIS

# Write where statement for site visits
input_sql = site_visit_data %>%
  # Drop geometry
  st_drop_geometry() %>%
  select(st_vst) %>%
  # Format site visit codes
  mutate(st_vst = paste('\'', st_vst, '\'', sep = '')) %>%
  # Collapse rows
  summarize(st_vst=paste(st_vst,collapse=", ")) %>%
  # Pull result out of dataframe
  pull(st_vst)
where_statement = paste('\r\nWHERE site_visit.site_visit_code IN (',
                        input_sql,
                        ');',
                        sep = '')

# Read vegetation cover data from AKVEG Database for selected site visits
vegetation_query = read_file(vegetation_file) %>%
  # Modify query with where statement
  str_replace(., ';', where_statement)
vegetation_data = as_tibble(dbGetQuery(database_connection, vegetation_query))

# Read abiotic top cover data from AKVEG Database for selected site visits
abiotic_query = read_file(abiotic_file) %>%
  # Modify query with where statement
  str_replace(., ';', where_statement)
abiotic_data = as_tibble(dbGetQuery(database_connection, abiotic_query))

# Read whole tussock cover data from AKVEG Database for selected site visits
tussock_query = read_file(tussock_file) %>%
  # Modify query with where statement
  str_replace(., ';', where_statement)
tussock_data = as_tibble(dbGetQuery(database_connection, tussock_query))

# Read ground cover data from AKVEG Database for selected site visits
ground_query = read_file(ground_file) %>%
  # Modify query with where statement
  str_replace(., ';', where_statement)
ground_data = as_tibble(dbGetQuery(database_connection, ground_query))

# Read structural group cover data from AKVEG Database for selected site visits
structural_query = read_file(structural_file) %>%
  # Modify query with where statement
  str_replace(., ';', where_statement)
structural_data = as_tibble(dbGetQuery(database_connection, structural_query))

# Read shrub structure data from AKVEG Database for selected site visits
shrub_query = read_file(shrub_file) %>%
  # Modify query with where statement
  str_replace(., ';', where_statement)
shrub_data = as_tibble(dbGetQuery(database_connection, shrub_query))

# Read environment data from AKVEG Database for selected site visits
environment_query = read_file(environment_file) %>%
  # Modify query with where statement
  str_replace(., ';', where_statement)
environment_data = as_tibble(dbGetQuery(database_connection, environment_query))

# Read soil metrics data from AKVEG Database for selected site visits
soilmetrics_query = read_file(soilmetrics_file) %>%
  str_replace(., ';', where_statement)
soilmetrics_data = as_tibble(dbGetQuery(database_connection, soilmetrics_query))

# Check number of cover observations per project
project_data = vegetation_data %>%
  left_join(site_visit_data, join_by('st_vst')) %>%
  group_by(prjct_cd) %>%
  summarize(obs_n = n())

# Export data to csv files
taxa_data %>%
  write.csv(., file = taxa_output, fileEncoding = 'UTF-8', row.names = FALSE)
site_visit_data %>%
  st_drop_geometry() %>%
  write.csv(., file = site_visit_output, fileEncoding = 'UTF-8', row.names = FALSE)
vegetation_data %>%
  write.csv(., file = vegetation_output, fileEncoding = 'UTF-8', row.names = FALSE)
abiotic_data %>%
  write.csv(., file = abiotic_output, fileEncoding = 'UTF-8', row.names = FALSE)
tussock_data %>%
  write.csv(., file = tussock_output, fileEncoding = 'UTF-8', row.names = FALSE)
ground_data %>%
  write.csv(., file = ground_output, fileEncoding = 'UTF-8', row.names = FALSE)
structural_data %>%
  write.csv(., file = structural_output, fileEncoding = 'UTF-8', row.names = FALSE)
shrub_data %>%
  write.csv(., file = shrub_output, fileEncoding = 'UTF-8', row.names = FALSE)
environment_data %>%
  write.csv(., file = environment_output, fileEncoding = 'UTF-8', row.names = FALSE)
soilmetrics_data %>%
  write.csv(., file = soilmetrics_output, fileEncoding = 'UTF-8', row.names = FALSE)
