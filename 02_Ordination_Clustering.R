# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# Clustering for ACCS Nelchina project
# Author: Timm Nawrocki, Alaska Center for Conservation Science
# Last Updated: 2025-03-28
# Usage: Script should be executed in R 4.4.3+.
# Description: "Clustering for ACCS Nelchina project" assigns clusters to data collected in the eastern Talkeetna Mountains and Nelchina Basin in 2022-2023.
# ---------------------------------------------------------------------------

# Import required libraries
library(dplyr)
library(fs)
library(ggplot2)
library(metR)
library(janitor)
library(lubridate)
library(readr)
library(readxl)
library(writexl)
library(stringr)
library(tibble)
library(tidyr)
library(cluster)
library(vegan)
library(vegan3d)
library(vegclust)
library(rgl)
library(indicspecies)
library(viridis)

#### SET UP DIRECTORIES AND FILES
####------------------------------

# Set random seed
set.seed(314)

# Set root directory (modify to your folder structure)
root_folder = 'ACCS_Work'

# Define input folders (modify to your folder structure)
database_repository = path('C:', root_folder, 'Repositories/akveg-database')
credentials_folder = path('C:', root_folder, 'Credentials/akveg_private_read')
project_folder = path('C:', root_folder, 'Projects/VegetationEcology/USNVC_SouthernBoreal_Alpine/Data')
input_folder = path(project_folder, 'Data_Input/plot_data')
output_folder = path(project_folder, 'Data_Output')
plot_folder = path(output_folder, 'plots')

# Define input files
taxa_input = path(input_folder, '00_taxonomy.csv')
site_visit_input = path(input_folder, '03_site_visit.csv')
vegetation_input = path(input_folder, '05_vegetation.csv')
abiotic_input = path(input_folder, '06_abiotic_top_cover.csv')
tussock_input = path(input_folder, '07_whole_tussock_cover.csv')
ground_input = path(input_folder, '08_ground_cover.csv')
structural_input = path(input_folder, '09_structural_group_cover.csv')
shrub_input = path(input_folder, '11_shrub_structure.csv')
environment_input = path(input_folder, '12_environment.csv')
soilmetrics_input = path(input_folder, '13_soil_metrics.csv')
reassignment_input = path(output_folder, 'manual_reassignments.xlsx')

# Define output files
comparison_output = path(output_folder, 'cluster_comparison.xlsx')
correlation_output = path(output_folder, 'correlation_comparison.xlsx')
constancy_output = path(output_folder, 'constancy_cover.xlsx')
cluster_output = path(output_folder, 'cluster_assignments.xlsx')

#### DEFINE FUNCTIONS
####------------------------------

# Define a function to compare variance and silhouette widths across multiple cluster numbers for fuzzy clustering with noise
fuzzy_nc_compare = function(vegetation_matrix, initial_n, final_n) {
  # Calculate bray dissimilarity
  bray_distance = vegdist(vegetation_matrix, method="bray")
  
  # Create empty data frame to store results
  cluster_results = tibble(
    cluster_n = numeric(),
    cluster = character(),
    variance = double(),
    avg_sil = double()
  )
  
  # Loop through cluster values to test cluster variance
  while (initial_n <= final_n) {
    print(paste('Conducting clustering with', {initial_n}, 'clusters...', sep = ' '))
    
    # Conduct clustering with n clusters
    nc_results = vegclust(x = vegetation_matrix,
                          mobileCenters = initial_n, 
                          method = 'NC',
                          m = 1.2,
                          dnoise = 0.8,
                          nstart = 50)
    
    # Assign hard clusters
    nc_clusters = tibble(defuzzify(nc_results)$cluster) %>%
      rename(cluster = `defuzzify(nc_results)$cluster`) %>%
      mutate(cluster_int = case_when(cluster == 'N' ~ -1.0,
                                     cluster != 'N' ~ as.double(str_replace(cluster, 'M', '')),
                                     TRUE ~ -999))
    
    # Calculate silhouette widths
    sil_widths = silhouette(nc_clusters$cluster_int, bray_distance)
    sil_width = tibble(sil_widths[, "sil_width"]) %>%
      rename(sil_width = `sil_widths[, "sil_width"]`)
    sil_data = cbind(nc_clusters, sil_width) %>%
      group_by(cluster) %>%
      summarize(avg_sil = mean(sil_width))
    
    # Calculate cluster variance
    clust_var = enframe(clustvar(bray_distance, nc_clusters$cluster)) %>%
      rename(cluster = name, variance = value) %>%
      mutate(cluster_n = initial_n) %>%
      left_join(sil_data, by = 'cluster')
    
    # Bind rows to cluster results
    cluster_results = rbind(cluster_results, clust_var)
    
    # Increase counter
    initial_n = initial_n + 1
  }
  
  # Return cluster results
  return(cluster_results)
  
}

# Define a function to calculate correlations of environmental variables with NMDS axes
calculate_nmds_correlations = function(nmds_results, site_data, variable_list) {
  # Create empty data frame to store results
  correlation_results = tibble(
    nmds_axis = numeric(),
    variable = character(),
    correlation = double(),
    p_value = double()
  )
  
  # Loop through all axis & variable combinations
  i = 1
  # Loop through all NMDS axes
  while (i <= nmds_results[['ndim']]) {
    # Extract NMDS axis scores
    nmds_axis = scores(nmds_results, display = "sites")[, i]
    # Loop through all variables
    for (variable_name in variable_list) {
      # Perform a correlation test between the axis and the variable
      correlation_test = cor.test(nmds_axis, site_data[[variable_name]], method = "pearson")
      # Store results in a data frame
      correlation_row = data.frame(nmds_axis=c(i),
                                   variable=c(variable_name),
                                   correlation = c(correlation_test[['estimate']][['cor']]),
                                   p_value = c(correlation_test[['p.value']]))
      # Append row to data frame of results
      correlation_results = rbind(correlation_results, correlation_row)
    }
    
    # Increase counter
    i = i + 1
    
  }
  
  # Adjust p values for multiple comparisons
  p_value_adjusted = p.adjust(as.vector(correlation_results$p_value),
                              method = 'bonferroni',
                              n = length(correlation_results$p_value))
  
  # Add adjusted p values to data frame
  correlation_results = correlation_results %>%
    mutate(p_value_adj = p_value_adjusted) %>%
    arrange(nmds_axis, p_value_adj)
  
  # Return correlation results
  return(correlation_results)
  
}

# Define a function to calculate a constancy and cover table
calculate_constancy = function(vegetation_data, site_data, cluster_variable,
                               cover_threshold=3, constancy_threshold=60) {
  # Format cluster data
  cluster_data = site_data %>%
    # Select columns
    select(st_vst, !!as.name(cluster_variable))
  
  # Count number of sites per cluster
  cluster_count = cluster_data %>%
    group_by(!!as.name(cluster_variable)) %>%
    summarize(cluster_count = n()) %>%
    ungroup()
  
  # Join vegetation data
  constancy_data = cluster_data %>%
    # Join vegetation data
    left_join(vegetation_data, by = 'st_vst') %>%
    # Group by site
    group_by(!!as.name(cluster_variable), taxon_code) %>%
    # Summarize constancy and cover
    summarize(cvr_mean = mean(cvr_pct),
              cvr_std = sd(cvr_pct),
              cvr_min = min(cvr_pct),
              cvr_max = max(cvr_pct),
              cvr_median=quantile(cvr_pct,probs=0.5),
              taxon_count = n()) %>%
    # Calculate constancy
    left_join(cluster_count, by = cluster_variable) %>%
    mutate(constancy = round((taxon_count/cluster_count)*100)) %>%
    # Filter to constancy and cover thresholds
    filter(cvr_mean >= cover_threshold) %>%
    filter(constancy >= constancy_threshold) %>%
    # Arrange rows
    arrange(!!as.name(cluster_variable), desc(constancy), desc(cvr_mean)) %>%
    ungroup()
    
    return(constancy_data)
    
}

# Define a function to calculate and plot an ordination surface
calculate_ordination_surface = function(nmds_results, environment_data, plot_data, cluster_variable,
                                        x_axis = 1, y_axis = 2, file_output, plot_title, variable_label,
                                        cb_palette, shape_scale, width = 6, height = 8) {
  # Calculate ordination surface
  surface_result = ordisurf(nmds_normalized, environment_data, choices = c(x_axis, y_axis), plot=FALSE)
  
  # Store summary
  surface_summary = summary(surface_result)
  
  # Create surface data
  grid_result = surface_result$grid
  surface_data = expand.grid(x = grid_result$x, y = grid_result$y)
  surface_data$z = as.vector(grid_result$z)
  surface_data = data.frame(na.omit(surface_data))
  
  # Store results in vector
  result_vector = c(surface_data, surface_summary)
  
  # Create output plot
  plot_output = ggplot() +
    metR::geom_contour2(data=surface_data, aes(x=x, y=y, z=z, label=after_stat(level))) +
    geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
    geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
    geom_point(data=plot_data, 
               aes(x=!!as.name(paste('nmds_axis_', x_axis, sep = '')),
                   y=!!as.name(paste('nmds_axis_', y_axis, sep='')),
                   color=factor(!!as.name(cluster_variable)),
                   shape=factor(!!as.name(cluster_variable))),
               size=2) +
    xlab(paste('NMDS Axis', toString(x_axis), sep = ' ')) +
    ylab(paste('NMDS Axis', toString(y_axis), sep = ' ')) +
    labs(title=plot_title,
         color=variable_label,
         shape=variable_label) +
    theme_light() +
    theme(legend.position='bottom', legend.direction='vertical') +
    scale_color_manual(values = cb_palette) +
    scale_shape_manual(values = shape_scale)
  ggsave(file_output,
         plot = plot_output,
         device = 'jpeg',
         path = NULL,
         scale = 1,
         width = width,
         height = height,
         units = 'in',
         dpi = 600,
         limitsize = TRUE)
  
  return(result_vector)
}

#### FORMAT DATA
####------------------------------

# Omit canyon plot (430)
omit_plots = c('NLC_430_20230723')

# Read data into dataframes
taxa_data = read_csv(taxa_input)
site_visit_data = read_csv(site_visit_input) %>%
  filter(!st_vst %in% omit_plots)
vegetation_data = read_csv(vegetation_input)%>%
  # Standardize taxa
  mutate(name_accepted = case_when(code_accepted == 'betnan' ~ 'Betula nana ssp. exilis',
                                   code_accepted == 'calthpal' ~ 'Caltha palustris ssp. radicans',
                                   code_accepted == 'camlas' ~ 'Campanula lasiocarpa ssp. latisepala',
                                   code_accepted == 'dryaja' ~ 'Dryas ajanensis ssp. beringensis',
                                   code_accepted == 'rhotom' ~ 'Rhododendron tomentosum ssp. decumbens',
                                   TRUE ~ name_accepted)) %>%
  mutate(code_accepted = case_when(code_accepted == 'betnan' ~ 'betnansexi',
                                   code_accepted == 'calthpal' ~ 'calpalsrad',
                                   code_accepted == 'camlas' ~ 'camlasslat',
                                   code_accepted == 'dryaja' ~ 'dryajasber',
                                   code_accepted == 'rhotom' ~ 'rhotomsdec',
                                   TRUE ~ code_accepted)) %>%
  # Convert dead vegetation to unique names
  mutate(taxon_code = case_when(dead_status == TRUE ~ paste(code_accepted, 'dead', sep='#'),
                                TRUE ~ code_accepted)) %>%
  # Convert trace values to 0.1%
  mutate(cvr_pct = case_when(cvr_pct == 0 ~ 0.1,
                             TRUE ~ cvr_pct)) %>%
  # Ensure no duplicate records from name changes
  group_by(st_vst, code_accepted, name_accepted, dead_status, cvr_type, taxon_code) %>%
  summarize(cvr_pct = sum(cvr_pct)) %>%
  # Round cover values to nearest 0.1%
  mutate(cvr_pct = round(cvr_pct, 1)) %>%
  # Omit plots
  filter(!st_vst %in% omit_plots) %>%
  # Order by site visit
  arrange(st_vst) %>%
  # Ungroup data
  ungroup()
abiotic_data = read_csv(abiotic_input) %>%
  filter(!st_vst %in% omit_plots)
tussock_data = read_csv(tussock_input) %>%
  filter(!st_vst %in% omit_plots)
ground_data = read_csv(ground_input) %>%
  filter(!st_vst %in% omit_plots)
structural_data = read_csv(structural_input) %>%
  filter(!st_vst %in% omit_plots)
shrub_data = read_csv(shrub_input)
environment_data = read_csv(environment_input) %>%
  filter(!st_vst %in% omit_plots)
soilmetrics_data = read_csv(soilmetrics_input) %>%
  filter(!st_vst %in% omit_plots)

# Summarize vegetation names and codes
vegetation_name_code = vegetation_data %>%
  distinct(code_accepted, name_accepted) %>%
  arrange(name_accepted)

# Add vegetation summaries
betshr_data = vegetation_data %>%
  filter(dead_status == FALSE) %>%
  filter(code_accepted %in% c('betnansexi', 'betgla')) %>%
  group_by(st_vst, dead_status, cvr_type) %>%
  summarize(cvr_pct = sum(cvr_pct)) %>%
  mutate(code_accepted = 'betshr',
         name_accepted = 'Betula shrubs',
         taxon_code = 'betshr') %>%
  ungroup()
dryas_data = vegetation_data %>%
  filter(dead_status == FALSE) %>%
  filter(code_accepted %in% c('dryala', 'dryajasber', 'dryhoo', 'dryint')) %>%
  group_by(st_vst, dead_status, cvr_type) %>%
  summarize(cvr_pct = sum(cvr_pct)) %>%
  mutate(code_accepted = 'dryas',
         name_accepted = 'Dryas dwarf shrubs',
         taxon_code = 'dryas') %>%
  ungroup()
dsalix_data = vegetation_data %>%
  filter(dead_status == FALSE) %>%
  filter(code_accepted %in% c('salarc', 'salpol', 'salret')) %>%
  group_by(st_vst, dead_status, cvr_type) %>%
  summarize(cvr_pct = sum(cvr_pct)) %>%
  mutate(code_accepted = 'dsalix',
         name_accepted = 'Salix dwarf shrubs',
         taxon_code = 'dsalix') %>%
  ungroup()
ndsalix_data = vegetation_data %>%
  filter(dead_status == FALSE) %>%
  filter(code_accepted %in% c('salala', 'salarb', 'salbarc', 'salbart', 'salfus', 'salgla',
                              'salhas', 'salpul', 'salric')) %>%
  group_by(st_vst, dead_status, cvr_type) %>%
  summarize(cvr_pct = sum(cvr_pct)) %>%
  mutate(code_accepted = 'ndsalix',
         name_accepted = 'Salix non-dwarf shrubs',
         taxon_code = 'ndsalix') %>%
  ungroup()
nerishr_data = vegetation_data %>%
  filter(dead_status == FALSE) %>%
  filter(code_accepted %in% c('castet', 'empnig')) %>%
  group_by(st_vst, dead_status, cvr_type) %>%
  summarize(cvr_pct = sum(cvr_pct)) %>%
  mutate(code_accepted = 'nerishr',
         name_accepted = 'Needleleaf ericaceous dwarf shrubs',
         taxon_code = 'nerishr') %>%
  ungroup()
erishr_data = vegetation_data %>%
  filter(dead_status == FALSE) %>%
  left_join(taxa_data, by = join_by('code_accepted' == 'code_akveg')) %>%
  filter(taxon_family == 'Ericaceae') %>%
  group_by(st_vst, dead_status, cvr_type) %>%
  summarize(cvr_pct = sum(cvr_pct)) %>%
  mutate(code_accepted = 'erishr',
         name_accepted = 'Ericaceous shrubs',
         taxon_code = 'erishr') %>%
  ungroup()
wetsed_data = vegetation_data %>%
  filter(dead_status == FALSE) %>%
  filter(code_accepted %in% c('caraqu', 'carlac', 'carmem', 'carrar', 'carrot', 'carsax', 'carutr', 'carvag')) %>%
  group_by(st_vst, dead_status, cvr_type) %>%
  summarize(cvr_pct = sum(cvr_pct)) %>%
  mutate(code_accepted = 'wetsed',
         name_accepted = 'Wetland sedges',
         taxon_code = 'wetsed') %>%
  ungroup()
sphagn_data = vegetation_data %>%
  filter(dead_status == FALSE) %>%
  left_join(taxa_data, by = join_by('code_accepted' == 'code_akveg')) %>%
  filter(taxon_genus == 'Sphagnum') %>%
  group_by(st_vst, dead_status, cvr_type) %>%
  summarize(cvr_pct = sum(cvr_pct)) %>%
  mutate(code_accepted = 'sphagn',
         name_accepted = 'Sphagnum',
         taxon_code = 'sphagn') %>%
  ungroup()
moss_data = vegetation_data %>%
  filter(dead_status == FALSE) %>%
  left_join(taxa_data, by = join_by('code_accepted' == 'code_akveg')) %>%
  filter(taxon_habit == 'moss') %>%
  group_by(st_vst, dead_status, cvr_type) %>%
  summarize(cvr_pct = sum(cvr_pct)) %>%
  mutate(code_accepted = 'moss',
         name_accepted = 'Moss',
         taxon_code = 'moss') %>%
  ungroup()
lichen_data = vegetation_data %>%
  filter(dead_status == FALSE) %>%
  left_join(taxa_data, by = join_by('code_accepted' == 'code_akveg')) %>%
  filter(taxon_habit == 'lichen') %>%
  group_by(st_vst, dead_status, cvr_type) %>%
  summarize(cvr_pct = sum(cvr_pct)) %>%
  mutate(code_accepted = 'lichen',
         name_accepted = 'Lichen',
         taxon_code = 'lichen') %>%
  ungroup()
nonvascular_data = vegetation_data %>%
  filter(dead_status == FALSE) %>%
  left_join(taxa_data, by = join_by('code_accepted' == 'code_akveg')) %>%
  filter(taxon_habit %in% c('lichen', 'moss', 'liverwort')) %>%
  group_by(st_vst, dead_status, cvr_type, taxon_genus) %>%
  summarize(cvr_pct = sum(cvr_pct)) %>%
  rename(name_accepted = taxon_genus) %>%
  left_join(taxa_data, by = join_by('name_accepted' == 'taxon_name')) %>%
  filter(!taxon_category %in% c('unknown', 'functional group')) %>%
  mutate(code_accepted = code_akveg,
         taxon_code = code_akveg) %>%
  select(st_vst, code_accepted, name_accepted, dead_status, cvr_type, taxon_code, cvr_pct) %>%
  ungroup()
shrub_data = vegetation_data %>%
  filter(dead_status == FALSE) %>%
  left_join(taxa_data, by = join_by('code_accepted' == 'code_akveg')) %>%
  filter(taxon_habit %in% c('dwarf shrub', 'dwarf shrub, shrub', 'dwarf shrub, shrub, tree',
                            'shrub', 'shrub, deciduous tree')) %>%
  group_by(st_vst, dead_status, cvr_type) %>%
  summarize(cvr_pct = sum(cvr_pct)) %>%
  mutate(code_accepted = 'shrub',
         name_accepted = 'Shrub',
         taxon_code = 'shrub') %>%
  ungroup()
graminoid_data = vegetation_data %>%
  filter(dead_status == FALSE) %>%
  left_join(taxa_data, by = join_by('code_accepted' == 'code_akveg')) %>%
  filter(taxon_habit == 'graminoid') %>%
  group_by(st_vst, dead_status, cvr_type) %>%
  summarize(cvr_pct = sum(cvr_pct)) %>%
  mutate(code_accepted = 'graminoid',
         name_accepted = 'Graminoid',
         taxon_code = 'graminoid') %>%
  ungroup()
forb_data = vegetation_data %>%
  filter(dead_status == FALSE) %>%
  left_join(taxa_data, by = join_by('code_accepted' == 'code_akveg')) %>%
  filter(taxon_habit == 'forb') %>%
  group_by(st_vst, dead_status, cvr_type) %>%
  summarize(cvr_pct = sum(cvr_pct)) %>%
  mutate(code_accepted = 'forb',
         name_accepted = 'Forb',
         taxon_code = 'forb') %>%
  ungroup()
bare_data = abiotic_data %>%
  filter(abiotic_element %in% c('bedrock (exposed)', 'rock fragments', 'soil')) %>%
  group_by(st_vst) %>%
  summarize(cvr_pct = sum(cvr_pct)) %>%
  mutate(code_accepted = 'bare',
         name_accepted = 'Bare ground',
         taxon_code = 'bare',
         dead_status = FALSE,
         cvr_type = 'top abiotic cover') %>%
  ungroup()

#### ADD BARE GROUND
# Combine vegetation data
vegetation_combined = vegetation_data %>%
  left_join(taxa_data, by = join_by('code_accepted' == 'code_akveg')) %>%
  filter(taxon_genus != taxon_name) %>%
  filter(!taxon_category %in% c('unknown', 'functional group')) %>%
  filter(!taxon_habit %in% c('lichen', 'moss', 'liverwort')) %>%
  select(st_vst, code_accepted, name_accepted, dead_status, cvr_type, taxon_code, cvr_pct) %>%
  rbind(betshr_data,
        dryas_data,
        dsalix_data,
        ndsalix_data,
        nerishr_data,
        erishr_data,
        wetsed_data,
        sphagn_data,
        moss_data,
        lichen_data,
        nonvascular_data,
        shrub_data,
        graminoid_data,
        forb_data,
        bare_data) %>%
  group_by(st_vst, taxon_code) %>%
  summarize(cvr_pct = sum(cvr_pct)) %>%
  ungroup()

# Convert vegetation data to matrix
vegetation_matrix = vegetation_combined %>%
  # Convert to wide format
  select(st_vst, taxon_code, cvr_pct) %>%
  pivot_wider(names_from = taxon_code, values_from = cvr_pct) %>%
  # Convert NA values to zero
  replace(is.na(.), 0) %>%
  # Convert st_vst column to row names
  column_to_rownames(var='st_vst')

# Collapse soil metric data to single measurement closest to soil surface
soil_minimum = soilmetrics_data %>%
  group_by(st_vst) %>%
  slice(which.min(measure_depth_cm))

# Format site data
site_data = site_visit_data %>%
  # Join environment data
  left_join(environment_data, by = 'st_vst') %>%
  # Join soils data from minimum measurement depth
  left_join(soil_minimum, by = 'st_vst') %>%
  # Convert moisture regime to a numeric relative scale
  mutate(moisture_quant = case_when(moisture_regime == 'xeric' ~ 0,
                                    moisture_regime == 'mesic-xeric heterogenous' ~ 1,
                                    moisture_regime == 'mesic' ~ 2,
                                    moisture_regime == 'hygric-mesic heterogenous' ~ 3,
                                    moisture_regime == 'hygric' ~ 4,
                                    moisture_regime == 'hydric-hygric heterogenous' ~ 5,
                                    moisture_regime == 'hydric' ~ 6,
                                    moisture_regime == 'aquatic-hydric heterogenous' ~ 7,
                                    TRUE ~ NA)) %>%
  # Update no data in quantitative fields
  mutate(depth_15_percent_coarse_fragments_cm = case_when(depth_15_percent_coarse_fragments_cm == -999 ~ 80,
                                                       TRUE ~ depth_15_percent_coarse_fragments_cm)) %>%
  mutate(temperature_deg_c = case_when(temperature_deg_c == -999 ~ NA,
                                       TRUE ~ temperature_deg_c)) %>%
  mutate(depth_moss_duff_cm = case_when(depth_moss_duff_cm == -999 ~ NA,
                                        TRUE ~ depth_moss_duff_cm)) %>%
  mutate(depth_restrictive_layer_cm = case_when(depth_restrictive_layer_cm == -999 ~ 100,
                                                TRUE ~ depth_restrictive_layer_cm)) %>%
  mutate(ph = case_when(ph == -999 ~ NA,
                        TRUE ~ ph)) %>%
  mutate(conductivity_mus = case_when(conductivity_mus == -999 ~ NA,
                                      TRUE ~ conductivity_mus)) %>%
  # Order by site visit
  arrange(st_vst)

#### CONDUCT EXPLORATORY ORDINATION
####------------------------------

# Normalize vegetation matrix
vegetation_normalized = decostand(vegetation_matrix, method='normalize')

# Calculate bray dissimilarity on normalized vegetation matrix
bray_raw = vegdist(vegetation_matrix, method="bray")
bray_normalized = vegdist(vegetation_normalized, method="bray")

# Calculate NMDS ordination on dissimilarity matrix
nmds_raw = metaMDS(bray_raw,
                   distance = "bray",
                   k = 3,
                   maxit = 50, 
                   trymax = 100,
                   wascores = TRUE)
nmds_normalized = metaMDS(bray_normalized,
                          distance = "bray",
                          k = 3,
                          maxit = 50, 
                          trymax = 100,
                          wascores = TRUE)


# Compare stress plots
stressplot(nmds_raw, main="Stress Plot for Raw Cover/Bray-Curtis")
stressplot(nmds_normalized, main="Stress Plot for Normalized Cover/Bray-Curtis")
#### MANUALLY EVALUATE STRESS PLOTS TO SELECT RAW COVER OR ONE OF TWO NORMALIZATIONS

# Assign results
nmds_results = nmds_normalized
vegetation_select = vegetation_normalized

# Test correlations amond NMDS axes and variables
variable_list = c('elevation_m', 'aspect', 'slope_deg', 'heat_load', 'position', 'relief', 'roughness',
                  'wetness', 'moisture_quant', 'depth_moss_duff_cm', 'depth_restrictive_layer_cm',
                  'depth_15_percent_coarse_fragments_cm', 'ph', 'conductivity_mus', 'temperature_deg_c')
correlation_results = calculate_nmds_correlations(nmds_results, site_data, variable_list)
#### ADD SNOW-OFF DATE OR NUMBER OF ANNUAL SNOW-COVER DAYS IN FUTURE

# Export correlation results
write_xlsx(correlation_results,
           path = correlation_output,
           col_names=TRUE,
           format_headers = FALSE)

# Plot an interactive 3D ordination
nmds_sites = site_data %>%
  column_to_rownames(var='st_vst')
ordination_rgl = ordirgl(nmds_results, radius = .01, color="black")

# Prepare species scores
nmds_results.spp_scrs <- 
  sppscores(nmds_results) <- vegetation_select
species_count = colSums(vegetation_select > 0)
species_cover = colSums(vegetation_select)
species_adjusted = species_cover/species_count

# Plot NMDS axes with species scores
x11()
plot(nmds_results, type = "n", choices = c(1, 2))
points(nmds_results, display = "sites", cex = 0.7, pch=21, col="red", bg="yellow")
orditorp(nmds_results, "sp", priority = species_adjusted, pch="+", pcol="blue")

#### CONDUCT FUZZY CLUSTERING WITH NOISE
####------------------------------

# Compare noise clustering with different cluster numbers
cluster_results = fuzzy_nc_compare(vegetation_select, 2, 20)

# Format cluster results
cluster_variance = cluster_results %>%
  select(cluster, variance, cluster_n) %>%
  pivot_wider(names_from = cluster, values_from = variance)
cluster_comparison = cluster_results %>%
  filter(cluster != 'N') %>%
  group_by(cluster_n) %>%
  summarize(mean_variance = mean(variance),
            mean_sil = mean(avg_sil)) %>%
  ungroup() %>%
  left_join(cluster_variance, by = 'cluster_n') %>%
  arrange(desc(mean_sil))

# Export cluster comparison
write_xlsx(cluster_comparison,
           path = comparison_output,
           col_names=TRUE,
           format_headers = FALSE)
#### MANUALLY EVALUATE CLUSTER COMPARISON AND SELECT CLUSTER_N

# Conduct clustering with n clusters
nc_results = vegclust(x = vegetation_normalized,
                      mobileCenters = 19, 
                      method = 'NC',
                      m = 1.2,
                      dnoise = 0.8,
                      nstart = 50)

# Extract noise cluster membership
nc_membership = tibble(nc_results$memb) %>%
  mutate_if(is.numeric, round, 2)

# Assign hard clusters to sites
nc_clusters = tibble(defuzzify(nc_results)$cluster) %>%
  rename(cluster = `defuzzify(nc_results)$cluster`)

# Bind clusters and membership to site data
analysis_data = site_data %>%
  mutate(nmds_axis_1 = scores(nmds_normalized, display = "sites")[, 1]) %>%
  mutate(nmds_axis_2 = scores(nmds_normalized, display = "sites")[, 2]) %>%
  mutate(nmds_axis_3 = scores(nmds_normalized, display = "sites")[, 3]) %>%
  cbind(nc_clusters) %>%
  cbind(nc_membership) %>%
  mutate(cluster_int = case_when(cluster == 'N' ~ -1.0,
                                 cluster != 'N' ~ as.double(str_replace(cluster, 'M', '')),
                                 TRUE ~ -999))

# Export table with reassignments
reassignment_data = read_xlsx(reassignment_input, sheet='reassignment')
vegetation_join = vegetation_matrix %>%
  rownames_to_column('st_vst')
cluster_assignments = analysis_data %>%
  select(st_vst, cluster_int) %>%
  left_join(reassignment_data, by = 'st_vst') %>%
  left_join(vegetation_join, by = 'st_vst') %>%
  select(st_vst, cluster_int, alliance_v1, alliance_n_v1, comment,
         betshr, nerishr, erishr, dryas, dsalix, ndsalix, graminoid, forb) %>%
  arrange(cluster_int, alliance_v1, st_vst)
write_xlsx(cluster_assignments,
           path = cluster_output,
           col_names=TRUE,
           format_headers = FALSE)

# Prepare data for final analyses and plotting
final_data = cluster_assignments %>%
  select(-cluster_int) %>%
  left_join(analysis_data, by = 'st_vst')

# Indicator species analysis
preliminary_indicators = multipatt(vegetation_matrix,
                                   final_data$alliance_n_v1,
                                   control = how(nperm=100)) 
summary(preliminary_indicators, indvalcomp=TRUE)

# Calculate constancy and cover
constancy_data = calculate_constancy(vegetation_combined, final_data, 'alliance_v1')

# Export constancy and cover table
write_xlsx(constancy_data,
           path = constancy_output,
           col_names=TRUE,
           format_headers = FALSE)

# Set color-blind friendly palette
shape_scale = c(15,15,15,15,15,17,17,17,17,17,19,19,19,19,19)
cb_palette_15 <- c(
  "#E69F00", # Orange
  "#56B4E9", # Sky Blue
  "#009E73", # Bluish Green
  "#F0E442", # Yellow
  "#0072B2", # Dark Blue
  "#D55E00", # Vermillion (Red-Orange)
  "#CC79A7", # Reddish Purple
  "#999999", # Gray
  "#882255", # Dark Red
  "#44AA99", # Teal
  "#117733", # Dark Green
  "#DDCC77", # Light Yellow-Brown
  "#AA4499", # Pinkish Purple
  "#661100", # Dark Brown
  "#88CCEE"  # Light Blue
)

# Export plot for moisture (axes 1 & 2)
surface_moisture_1_2 = calculate_ordination_surface(nmds_results,
                                                    site_data$moisture_quant,
                                                    final_data,
                                                    'alliance_v1',
                                                    1, 2,
                                                    path(plot_folder, 'moisture_1_2.jpg'),
                                                    'Moisture Gradient Relative to NMDS Axes 1 & 2',
                                                    'Clusters',
                                                    cb_palette_15,
                                                    shape_scale,
                                                    6.5, 9)
surface_moisture_1_3 = calculate_ordination_surface(nmds_results,
                                                    site_data$moisture_quant,
                                                    final_data,
                                                    'alliance_v1',
                                                    1, 3,
                                                    path(plot_folder, 'moisture_1_3.jpg'),
                                                    'Moisture Gradient Relative to NMDS Axes 1 & 3',
                                                    'Clusters',
                                                    cb_palette_15,
                                                    shape_scale,
                                                    6.5, 9)
surface_wetness_1_2 = calculate_ordination_surface(nmds_results,
                                                   site_data$wetness,
                                                   final_data,
                                                   'alliance_v1',
                                                   1, 2,
                                                   path(plot_folder, 'wetness_1_2.jpg'),
                                                   'Topographic Wetness Relative to NMDS Axes 1 & 2',
                                                   'Clusters',
                                                   cb_palette_15,
                                                   shape_scale,
                                                   6.5, 9)
surface_wetness_1_3 = calculate_ordination_surface(nmds_results,
                                                   site_data$wetness,
                                                   final_data,
                                                   'alliance_v1',
                                                   1, 3,
                                                   path(plot_folder, 'wetness_1_3.jpg'),
                                                   'Topographic Wetness Relative to NMDS Axes 1 & 3',
                                                   'Clusters',
                                                   cb_palette_15,
                                                   shape_scale,
                                                   6.5, 9)
surface_coarsefrag_1_2 = calculate_ordination_surface(nmds_results,
                                                      site_data$depth_15_percent_coarse_fragments_cm,
                                                      final_data,
                                                      'alliance_v1',
                                                      1, 2,
                                                      path(plot_folder, 'coarsefrag_1_2.jpg'),
                                                      'Depth to 15% Coarse Fragments Relative to NMDS Axes 1 & 2',
                                                      'Clusters',
                                                      cb_palette_15,
                                                      shape_scale,
                                                      6.5, 9)
surface_coarsefrag_1_3 = calculate_ordination_surface(nmds_results,
                                                   site_data$depth_15_percent_coarse_fragments_cm,
                                                   final_data,
                                                   'alliance_v1',
                                                   1, 3,
                                                   path(plot_folder, 'coarsefrag_1_3.jpg'),
                                                   'Depth to 15% Coarse Fragments Relative to NMDS Axes 1 & 3',
                                                   'Clusters',
                                                   cb_palette_15,
                                                   shape_scale,
                                                   6.5, 9)
surface_restrict_1_2 = calculate_ordination_surface(nmds_results,
                                                    site_data$depth_restrictive_layer_cm,
                                                    final_data,
                                                    'alliance_v1',
                                                    1, 2,
                                                    path(plot_folder, 'restrictive_1_2.jpg'),
                                                    'Depth Restrictive Layer Relative to NMDS Axes 1 & 2',
                                                    'Clusters',
                                                    cb_palette_15,
                                                    shape_scale,
                                                    6.5, 9)
surface_restrict_1_3 = calculate_ordination_surface(nmds_results,
                                                    site_data$depth_restrictive_layer_cm,
                                                    final_data,
                                                    'alliance_v1',
                                                    1, 3,
                                                    path(plot_folder, 'restrictive_1_3.jpg'),
                                                    'Depth Restrictive Layer Relative to NMDS Axes 1 & 3',
                                                    'Clusters',
                                                    cb_palette_15,
                                                    shape_scale,
                                                    6.5, 9)
surface_elevation_2_1 = calculate_ordination_surface(nmds_results,
                                              site_data$elevation_m,
                                              final_data,
                                              'alliance_v1',
                                              2, 1,
                                              path(plot_folder, 'elevation_2_1.jpg'),
                                              'Elevation Relative to NMDS Axes 2 & 1',
                                              'Clusters',
                                              cb_palette_15,
                                              shape_scale,
                                              6.5, 9)
surface_elevation_2_3 = calculate_ordination_surface(nmds_results,
                                                     site_data$elevation_m,
                                                     final_data,
                                                     'alliance_v1',
                                                     2, 3,
                                                     path(plot_folder, 'elevation_2_3.jpg'),
                                                     'Elevation Relative to NMDS Axes 2 & 3',
                                                     'Clusters',
                                                     cb_palette_15,
                                                     shape_scale,
                                                     6.5, 9)
surface_slope_3_1 = calculate_ordination_surface(nmds_results,
                                                 site_data$slope_deg,
                                                 final_data,
                                                 'alliance_v1',
                                                 3, 1,
                                                 path(plot_folder, 'slope_3_1.jpg'),
                                                 'Slope (deg) Relative to NMDS Axes 3 & 1',
                                                 'Clusters',
                                                 cb_palette_15,
                                                 shape_scale,
                                                 6.5, 9)
surface_slope_3_2 = calculate_ordination_surface(nmds_results,
                                                 site_data$slope_deg,
                                                 final_data,
                                                 'alliance_v1',
                                                 3, 2,
                                                 path(plot_folder, 'slope_3_2.jpg'),
                                                 'Slope (deg) Relative to NMDS Axes 3 & 2',
                                                 'Clusters',
                                                 cb_palette_15,
                                                 shape_scale,
                                                 6.5, 9)