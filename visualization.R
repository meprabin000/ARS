library(sf)
library(data.table)

# uncorrected standard deviation
usd <- function(data) {
  n <- length(data)
  return (sd(data) * sqrt((n-1)/n))
}

# removes data points whose value is greater or less than mean +- standard deviation 
clean_sd <- function(data, var, sd_no){
  data <- subset(data, var >= mean(var, na.rm = TRUE) - sd_no*sd(var, na.rm = TRUE) &  var <= mean(var, na.rm = TRUE) + sd_no*sd(var, na.rm = TRUE))
  return(data)
}

# SETS THE PATH TO THE SHAPE FILE ####
setwd("Documents/Agriculture Project/Analysis")
yield_2017_path <- "../04 Kolby/01 Farm Data/yieldHistory 2017-2020/2017/Jennings Far_Bailey Farm_42 ac_Harvest_2017-08-25_00.shp"
yield_2018_path <- "../04 Kolby/01 Farm Data/yieldHistory 2017-2020/2018/Jennings Far_Bailey Farm_42 ac_Harvest_2018-06-06_00.shp"
yield_2019_path <- "../04 Kolby/01 Farm Data/yieldHistory 2017-2020/2019/Jennings Far_Bailey Farm_42 ac_Harvest_2019-09-04_00.shp"
yield_2020_path <- "../04 Kolby/01 Farm Data/yieldHistory 2017-2020/2020/Jennings Far_Bailey Farm_42 ac_Harvest_2020-10-16_00.shp"

boundary_path <- "../04 Kolby/01 Farm Data/boundaries/Jennings Farms_Kolby Chrz_Bailey Farm_42 ac.shp"

# reading the shape files
boundary <- st_read(boundary_path)

# reading the yield shape files
yield_2017_sf <- st_read(yield_2017_path)
yield_2018_sf <- st_read(yield_2018_path)
yield_2019_sf <- st_read(yield_2019_path)
yield_2020_sf <- st_read(yield_2020_path)

# transforming the shape file
boundary <- st_transform(boundary, "+proj=utm +zone=11 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
yield_2017_sf <- st_transform(yield_2017_sf, st_crs(boundary))
yield_2018_sf <- st_transform(yield_2018_sf, st_crs(boundary))
yield_2019_sf <- st_transform(yield_2019_sf, st_crs(boundary))
yield_2020_sf <- st_transform(yield_2020_sf, st_crs(boundary))

# cleaned data using standard deviation greater than 3
yield_2017_sf_cl <- clean_sd(yield_2017_sf, yield_2017_sf$VRYIELDVOL, sd_no = 3)
yield_2018_sf_cl <- clean_sd(yield_2018_sf, yield_2017_sf$VRYIELDVOL, sd_no = 3)
yield_2019_sf_cl <- clean_sd(yield_2019_sf, yield_2017_sf$VRYIELDVOL, sd_no = 3)
yield_2020_sf_cl <- clean_sd(yield_2020_sf, yield_2017_sf$VRYIELDVOL, sd_no = 3)

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

# uncleaned cell_map
cell_map_2017 = get_cell_map(yield_2017_sf, boundary_with_cells)
cell_map_2018 = get_cell_map(yield_2018_sf, boundary_with_cells)
cell_map_2019 = get_cell_map(yield_2019_sf, boundary_with_cells)
cell_map_2020 = get_cell_map(yield_2020_sf, boundary_with_cells)

# cleaned cell_map
cell_map_2017_cl = get_cell_map(yield_2017_sf_cl, boundary_with_cells)
cell_map_2018_cl = get_cell_map(yield_2018_sf_cl, boundary_with_cells)
cell_map_2019_cl = get_cell_map(yield_2019_sf_cl, boundary_with_cells)
cell_map_2020_cl = get_cell_map(yield_2020_sf_cl, boundary_with_cells)

cell_map_yearly_list = data.frame(
  geometry = cell_map_2017_cl$geometry,
  yield_2017 = cell_map_2017_cl$yield_sum,
  yield_2018 = cell_map_2018_cl$yield_sum,
  yield_2019 = cell_map_2019_cl$yield_sum,
  yield_2020 = cell_map_2020_cl$yield_sum,
  yield_count_2017 = cell_map_2017_cl$yield_count,
  yield_count_2018 = cell_map_2018_cl$yield_count,
  yield_count_2019 = cell_map_2019_cl$yield_count,
  yield_count_2020 = cell_map_2020_cl$yield_count
)

# spatial mean and variability
get_spatial_mean <- function(var) {
  return (mean(var))
}

get_spatial_variability <- function(var) {
  return (usd(var))
}

spatial_mean_2017_cl <- get_spatial_mean(cell_map_yearly_list$yield_2017)
spatial_mean_2018_cl <- get_spatial_mean(cell_map_yearly_list$yield_2018)
spatial_mean_2019_cl <- get_spatial_mean(cell_map_yearly_list$yield_2019)
spatial_mean_2020_cl <- get_spatial_mean(cell_map_yearly_list$yield_2020)

spatial_variability_2017_cl <- get_spatial_variability(cell_map_yearly_list$yield_2017)
spatial_variability_2018_cl <- get_spatial_variability(cell_map_yearly_list$yield_2018)
spatial_variability_2019_cl <- get_spatial_variability(cell_map_yearly_list$yield_2019)
spatial_variability_2020_cl <- get_spatial_variability(cell_map_yearly_list$yield_2020)

# temporal mean and variability
get_temporal_mean_and_variability <- function(yearly_yields) {
  temporal_mean <- vector('numeric', length = length(yearly_yields[[1]]))
  temporal_sd <- vector('numeric', length = length(yearly_yields[[1]]))
  for (i in 1:length(yearly_yields[[1]])) {
    yearly_data_row <- vector('numeric', length = length(yearly_yields))
    for (j in 1:length(yearly_yields)) {
      yearly_data_row[j] <- yearly_yields[[j]][[i]]
    }
    print(yearly_data_row)
    temporal_mean[i] <- mean(yearly_data_row)
    temporal_sd[i] <- usd(yearly_data_row)
  }
  return(data.frame(
    temporal_mean = temporal_mean, 
    temporal_sd = temporal_sd
    ))
}

yearly_yield_list <- list(
  as.list(cell_map_yearly_list$yield_2017),
  as.list(cell_map_yearly_list$yield_2018),
  as.list(cell_map_yearly_list$yield_2019),
  as.list(cell_map_yearly_list$yield_2020)
)
temporal_analysis <- get_temporal_mean_and_variability(yearly_yield_list)
temporal_analysis <- cbind(temporal_analysis, geometry = cell_map_yearly_list$geometry)
temporal_analysis <- st_as_sf(temporal_analysis)

plotBreaks <- seq(min(0), max(3000), by = 50) 
plot(st_as_sf(cell_map_yearly_list)['yield_2020'], key.pos = 1, key.length = 1)

