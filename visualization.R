library(sf)
library(data.table)
library(ggplot2)
library(dplyr)

# uncorrected standard deviation
usd <- function(data) {
  n <- length(data)
  return (sd(data) * sqrt((n-1)/n))
}

# removes data points whose value is greater or less than mean +- standard deviation 
clean_sd <- function(data, var, sd_no){
  mean_ <- mean(var, na.rm = TRUE)
  sd_ <- sd(var, na.rm = TRUE)
  min_ <- (mean_ - sd_no*sd_)
  max_ <- (mean_ + sd_no*sd_)
  data$clean_col <- ifelse(data$VRYIELDVOL < min_, "low", 
                           ifelse(data$VRYIELDVOL > max_, "high",
                                  "right"))
  clean_data <- subset(data, clean_col == "right")
  outliers <- subset(data, clean_col != "right")
  return(list(clean_data, outliers))
}

# SETS THE PATH TO THE SHAPE FILE ####
setwd("Documents/Agriculture Project/Analysis")
yield_2017_path <- "../04 Kolby/01 Farm Data/yieldHistory 2017-2020/2017/Jennings Far_Bailey Farm_42 ac_Harvest_2017-08-25_00.shp"
yield_2018_path <- "../04 Kolby/01 Farm Data/yieldHistory 2017-2020/2018/Jennings Far_Bailey Farm_42 ac_Harvest_2018-06-06_00.shp"
yield_2019_path <- "../04 Kolby/01 Farm Data/yieldHistory 2017-2020/2019/Jennings Far_Bailey Farm_42 ac_Harvest_2019-09-04_00.shp"
yield_2020_path <- "../04 Kolby/01 Farm Data/yieldHistory 2017-2020/2020/Jennings Far_Bailey Farm_42 ac_Harvest_2020-10-16_00.shp"

boundary_path <- "../04 Kolby/01 Farm Data/boundaries/Jennings Farms_Kolby Chrz_Bailey Farm_42 ac.shp"

# reading the shape files and transforming the unit of zone to m
boundary <- st_read(boundary_path)
boundary <- st_transform(boundary, "+proj=utm +zone=11 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# all_yield_paths_data
all_yield_paths <- list(
  list(yield_2017_path, 2017),
  list(yield_2018_path, 2018),
  list(yield_2019_path, 2019),
  list(yield_2020_path, 2020)
  )

# transform and clean the data
data <- list()
for (i in seq(from = 1, to = length(all_yield_paths), by = 1)) {
  read_data <- st_read(all_yield_paths[[i]][[1]])
  transformed_data <- st_transform(read_data, st_crs(boundary))
  cleandata_outliers <- clean_sd(transformed_data, transformed_data$VRYIELDVOL, sd_no = 3)
  data[[i]] <- list(transformed_data, cleandata_outliers[[1]], all_yield_paths[[i]][[2]], cleandata_outliers[[2]])
}

# making 10m by 10m grids
grid <- st_make_grid(boundary, cellsize = c(10,10))
boundary_with_cells <- st_intersection(boundary$geometry, grid)

# for each yield point, cell_id column contains which cell does it belong to. cell_id value corresponds to boundary_with_cells index
get_cell_map = function (yield_sf, cell_boundary_sf) {
  intersected_grids <- st_within(st_as_sf(yield_sf), cell_boundary_sf)
  yield_with_cells <- cbind(yield_sf, cell_id = 0)
  cells_yield <- vector("numeric", length(cell_boundary_sf))
  cells_yield_count <- vector("numeric", length(cell_boundary_sf))
  for (point in 1:length(intersected_grids)) {
    cell_ind = as.integer(intersected_grids[point])
    cells_yield[cell_ind] <- cells_yield[cell_ind] + yield_sf$VRYIELDVOL[point]
    cells_yield_count[cell_ind] <- cells_yield_count[cell_ind] + 1
    setDT(yield_with_cells)[point, cell_id:=cell_ind]
  }
  
  # map that stores the cell map for each cell
  cell_map <- data.frame(
    geometry = cell_boundary_sf,
    yield_sum = cells_yield,
    yield_count = cells_yield_count
  )
  return (st_as_sf(cell_map))
}

# create a cell map for clean and unclean data
cell_maps <- list()
unclean_index = 1
clean_index = 2
year_index = 3
for (i in seq(from = 1, to = length(data), by = 1)) {
  unclean_map <- get_cell_map(data[[i]][[unclean_index]], boundary_with_cells)
  clean_map <- get_cell_map(data[[i]][[clean_index]], boundary_with_cells)
  year <- data[[i]][[year_index]]
  cell_maps[[i]] <- list(unclean_map, clean_map, year)
}

# cell_map_yearly_list
yield_sum_map <- data.frame(
  geometry = cell_maps[[1]][[1]]$geometry
)
for (i in seq(from = 1, to = length(cell_maps), by = 1)) {
  unclean_col_name <- paste('unclean_yield_sum_', cell_maps[[i]][[year_index]], sep = "")
  clean_col_name <- paste('clean_yield_sum_', cell_maps[[i]][[year_index]], sep = "")
  yield_sum_map[, unclean_col_name] <- cell_maps[[i]][[unclean_index]]$yield_sum
  yield_sum_map[, clean_col_name] <- cell_maps[[i]][[unclean_index]]$yield_sum
}

# spatial mean and variability
get_spatial_mean <- function(var) {
  return (mean(var))
}

get_spatial_variability <- function(var) {
  return (sd(var))
}

spatial_mean_variability_list <- list()
clean_columns <- names(yield_sum_map)[grep("^clean", names(yield_sum_map))]
unclean_columns <- names(yield_sum_map)[grep("^unclean", names(yield_sum_map))]

for (i in 1:length(data)) {
  mean_ <- get_spatial_mean(yield_sum_map[, unclean_columns[[i]]])
  sd_ <- get_spatial_variability(yield_sum_map[, unclean_columns[[i]]])
  spatial_mean_variability_list[[i]] <- list(mean_, sd_)
}

# temporal mean and variability
get_temporal_mean_and_variability <- function(yearly_yields) {
  temporal_mean <- vector('numeric', length = length(yearly_yields[[1]]))
  temporal_sd <- vector('numeric', length = length(yearly_yields[[1]]))
  for (i in 1:length(yearly_yields[[1]])) {
    yearly_data_row <- vector('numeric', length = length(yearly_yields))
    for (j in 1:length(yearly_yields)) {
      yearly_data_row[j] <- yearly_yields[[j]][[i]]
    }
    temporal_mean[i] <- mean(yearly_data_row)
    temporal_sd[i] <- sd(yearly_data_row)
  }
  return(data.frame(
    temporal_mean = temporal_mean, 
    temporal_sd = temporal_sd
    ))
}

temp <- yield_sum_map[, unclean_columns]
temporal_yield_list <- list()
for (i in 1:length(unclean_columns)) {
  col <- unclean_columns[[i]]
  temp_yield <- yield_sum_map[, col]
  temp_yield <- temp_yield - spatial_mean_variability_list[[i]][[1]]
  temporal_yield_list[[i]] <- temp_yield
}

temporal_analysis <- get_temporal_mean_and_variability(temporal_yield_list)
temporal_analysis <- cbind(temporal_analysis, geometry = yield_sum_map$geometry)
temporal_analysis <- st_as_sf(temporal_analysis)

plot(st_as_sf(temporal_analysis[, c('temporal_mean', 'geometry')]), key.pos = 1, key.length = 1)

############### outliers ############
# gets the outliers 
get_outliers <- function(unclean_data, clean_data) {
  outliers <- anti_join(as.data.frame(unclean_data), as.data.frame(clean_data), by = "geometry")
  return (st_as_sf(outliers))
}

# plots the outliers
plot_outliers <- function(outliers, year) {
  return (ggplot() +
    ggtitle(paste('Outliers ', as.character(year)), ) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_sf(data = st_as_sf(boundary_with_cells), fill = "lightblue") +
    geom_sf(data = st_as_sf(outliers), aes(color = outliers[["clean_col"]]), size = 2) +
    scale_color_manual(values = c("red", "blue"))
      )
}

get_directory_path <- function() {
  day <- format(Sys.time(), "%Y_%m_%d")
  path <- paste(day, "/", sep = "")
  if(!file.exists(path)) {
    dir.create(path)
  }
  return (path)
}

get_filename_with_timestamp <- function(filename, device) {
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  filename <- paste(c(filename, "_created_at_", timestamp, ".", device), collapse = "")
  return (filename)
}

# save the plots
for (i in 1:length(data)) {
  path <- get_directory_path()
  dev <- "png"
  filename <- get_filename_with_timestamp(paste("outliers_", data[[i]][[3]], sep = ""), dev)
  plt <- plot_outliers(data[[i]][[4]], data[[i]][[3]])
  ggsave(
    filename = filename,
    plot = plt,
    device = dev,
    dpi = 300,
    path = path
  )
}


# compressed value z-statistic
yield_z_stat_val <- data.frame(
  geometry = yield_sum_map[,'geometry']
)
for (i in 1:length(data)) {
  col_name_ <- unclean_columns[[i]]
  yield_ <- yield_sum_map[, col_name_]
  mean_ <- get_spatial_mean(yield_)
  sd_ <- get_spatial_variability(yield_)
  z_stat_ <- (yield_ - mean_) / sd_
  yield_z_stat_val[col_name_] <- z_stat_
}
yield_z_stat_val_sf <- st_as_sf(yield_z_stat_val)


# plotting z_stat for yield_sum
breaks <- c(-Inf, -2, 0, 2, (2+1e-10), Inf)
labels <- c('Q4', 'Q3', 'Q1', 'Q1', 'Q2')
plots <- list()
for (i in 1:length(data)) {
  col_name_ <- unclean_columns[[i]]
  bin_col_name_ <- paste('bin_', col_name_, sep = '')
  col_ <- yield_z_stat_val[, bin_col_name_]
  yield_z_stat_val[bin_col_name_] <- cut(yield_z_stat_val[,col_name_], breaks = breaks, labels = labels, right = FALSE, include.lowest = TRUE)
}
get_temporal_quartiles = function(row) {
  bin_quartiles <- list()
  q5 <- FALSE
  for (i in 1:length(data)) {
    col_name_ <- unclean_columns[[i]]
    bin_col_name_ <- paste('bin_', col_name_, sep = '')
    bin_quartiles[[i]] <- (row[bin_col_name_])
  }
  bin_qs <- unlist(bin_quartiles)[1]
  for (i in 1:length(bin_quartiles)) {
    bin <- unlist(bin_quartiles)[i]
    if (bin != bin_qs) {
      q5 <- TRUE
    }
  }
  
  if (q5) {
    return ("Q5")
  }
  
  return (bin_qs)
}


yield_z_stat_val <- as.data.frame(yield_z_stat_val)
yield_z_stat_val['temporal_quartiles'] <- apply(yield_z_stat_val, 1, get_temporal_quartiles)

col_ <- yield_z_stat_val$temporal_quartiles
print(ggplot() +
    geom_sf(data = yield_z_stat_val$geometry, aes(fill = col_)) + 
    scale_fill_manual(values = c('#ff0000ff', '#0000ffff', '#00ffffff', 'black')) +
    labs(title = bin_col_name_, fill="Yield cut off"))
