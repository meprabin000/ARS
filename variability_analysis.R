library(sf)
library(data.table)
library(ggplot2)
library(dplyr)
library(writexl)
library(tidyverse)
library(gstat)

#------------------------------------------------------------------------------------------------

# uncorrected standard deviation
usd <- function(data) {
  n <- length(data)
  return (sd(data) * sqrt((n-1)/n))
}

#------------------------------------------------------------------------------------------------

# removes data points whose value is greater or less than mean +- standard deviation 
label_outliers <- function(data, var, sd_no){
  mean_ <- mean(var, na.rm = TRUE)
  sd_ <- sd(var, na.rm = TRUE)
  min_ <- (mean_ - sd_no*sd_)
  max_ <- (mean_ + sd_no*sd_)
  data$clean_col <- ifelse(data$VRYIELDVOL < min_, "low", 
                           ifelse(data$VRYIELDVOL > max_, "high",
                                  "clean"))
  return(data)
}

get_outliers <- function(data) {
  return (subset(data, clean_col == "unclean"))
}
#------------------------------------------------------------------------------------------------

###############KOLBY FIELD######################
# SETS THE PATH TO THE SHAPE FILE ####
setwd("/Users/prabin/Documents/Agriculture Project/Analysis")
kolby_yield_2017_path <- "../04 Kolby/01 Farm Data/yieldHistory 2017-2020/2017/Jennings Far_Bailey Farm_42 ac_Harvest_2017-08-25_00.shp"
kolby_yield_2018_path <- "../04 Kolby/01 Farm Data/yieldHistory 2017-2020/2018/Jennings Far_Bailey Farm_42 ac_Harvest_2018-06-06_00.shp"
kolby_yield_2019_path <- "../04 Kolby/01 Farm Data/yieldHistory 2017-2020/2019/Jennings Far_Bailey Farm_42 ac_Harvest_2019-09-04_00.shp"
kolby_yield_2020_path <- "../04 Kolby/01 Farm Data/yieldHistory 2017-2020/2020/Jennings Far_Bailey Farm_42 ac_Harvest_2020-10-16_00.shp"

kolby_boundary_path <- "../04 Kolby/01 Farm Data/boundaries/Jennings Farms_Kolby Chrz_Bailey Farm_42 ac.shp"

# all_yield_paths_data
kolby_yield_paths <- list(
  list(kolby_yield_2017_path, 2017),
  list(kolby_yield_2017_path, 2018),
  list(kolby_yield_2019_path, 2019),
  list(kolby_yield_2020_path, 2020)
)
#------------------------------------------------------------------------------------------------

boundary_path = kolby_boundary_path
all_yield_paths = kolby_yield_paths
field_name = "Kolby"
#------------------------------------------------------------------------------------------------

# reading the shape files and transforming the unit of zone to m
boundary <- st_read(boundary_path)
boundary <- st_transform(boundary, "+proj=utm +zone=11 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
#------------------------------------------------------------------------------------------------

# transform and clean the data
data <- list()
for (i in seq(from = 1, to = length(all_yield_paths), by = 1)) {
  read_data <- st_read(all_yield_paths[[i]][[1]])
  transformed_data <- st_transform(read_data, st_crs(boundary))
  outlier_labeled_data <- label_outliers(transformed_data, transformed_data$VRYIELDVOL, sd_no = 3)
  data[[i]] <- list(outlier_labeled_data, all_yield_paths[[i]][[2]])
}
#------------------------------------------------------------------------------------------------

# making 10m by 10m grids
grid <- st_make_grid(boundary, cellsize = c(20,20))
boundary_with_cells <- st_intersection(boundary$geometry, grid)

#------------------------------------------------------------------------------------------------
# for each yield point, identify which cell it belongs to and label in column cell_id
map_yield_points_to_cell = function(yield_sf, cell_boundary_sf) {
  intersected_grids <- st_within(st_as_sf(yield_sf), cell_boundary_sf)
  yield_with_cells <- cbind(yield_sf, cell_id = 0)
  for (point in 1:length(intersected_grids)) {
    cell_ind = as.integer(intersected_grids[point])
    setDT(yield_with_cells)[point, cell_id:=cell_ind]
  }
  return (yield_with_cells)
}

# mapping data to yield points
data_index <- 1
for (i in 1:length(data)) {
  d <- map_yield_points_to_cell(data[[i]][[data_index]], boundary_with_cells)
  data[[i]][[data_index]] <- d
}
#------------------------------------------------------------------------------------------------
# impute the outliers using idw method, 
# Takes a long time to finish run: Could try paralleling the process
year_index <- 2
for (i in 1:length(data)) {
  wd <- data[[i]][[data_index]] # get the raw yield
  clean_data <- subset(wd, clean_col == "clean") # only use the clean data to impute the outliers
  print(paste(c('Imputing data of year ', data[[i]][[year_index]])))
  for( point in 1:nrow(wd)) {
    if (wd$clean_col[point] != "clean") {
      clean_points <- subset(clean_data, cell_id == wd$cell_id[point]) # get clean neighboring yield points
      if (nrow(clean_points) == 0) { # if no clean point is found for that specific cell, mark it as 'unclean'
        wd$clean_col[point] = 'unclean'
        wd$ImputedYield[point] <- wd$VRYIELDVOL[point]
        next # skip to next iteration
      }
      wd$clean_col[point] = 'imputed' # distinguish imputed data from original clean data
      gs <- gstat(formula = VRYIELDVOL~1, locations = st_as_sf(clean_points)) # idw model
      p <- predict(gs, st_as_sf(wd[point,]), debug.level = 0) # get the value for a given point
      wd$ImputedYield[point] = p[[1]]
    } else {
      wd$ImputedYield[point] <- wd$VRYIELDVOL[point]
    }
  }
  data[[i]][[data_index]] <- wd
}

#------------------------------------------------------------------------------------------------
# for each yield point, cell_id column contains which cell does it belong to. cell_id value corresponds to boundary_with_cells index
get_cell_map = function (yield_sf, cell_boundary_sf) {
  cell_ids <- sort(unique(yield_sf$cell_id))
  cells_yield <- vector("numeric", length(cell_boundary_sf))
  cells_yield_count <- vector("numeric", length(cell_boundary_sf))
  for (cell_ind in cell_ids) {
    cell_points <- subset(yield_sf, cell_id == cell_ind)
    cells_yield[cell_ind] <- sum(cell_points$ImputedYield)
    cells_yield_count[cell_ind] <- nrow(cell_points)
  }

  # map that stores the cell map for each cell
  cell_map <- data.frame(
    cell_id = (1:length(cell_boundary_sf)),
    geometry = cell_boundary_sf,
    yield_sum = cells_yield,
    yield_count = cells_yield_count
  )
  return (st_as_sf(cell_map))
}

# create a cell map for clean data
cell_maps <- list()
for (i in 1:length(data)) {
  unclean_map <- get_cell_map(data[[i]][[data_index]], boundary_with_cells)
  year <- data[[i]][[year_index]]
  cell_maps[[i]] <- list(unclean_map, year)
}

#------------------------------------------------------------------------------------------------
# remove cells that contain na values and outlier data
# yield_count_thresold <- 1
# for (i in 1:length(data)) {
#   filtered_data <- data[[i]][[1]][!is.na(data[[i]][[1]]$cell_id)]
#   cell_ids_with_outliers <- subset(filtered_data, clean_col == "unclean")$cell_id %>% unique()
#   cell_maps[[i]][[1]] <- cell_maps[[i]][[1]][!cell_maps[[i]][[1]]$cell_id %in% cell_ids_with_outliers, ]
#   cell_maps[[i]][[1]] <- subset(cell_maps[[i]][[1]], yield_count > yield_count_thresold)
# }
# plot(cell_maps[[4]][[1]])
#------------------------------------------------------------------------------------------------

# cell_map_yearly_list
yield_sum_map <- data.frame(
  geometry = cell_maps[[1]][[1]]$geometry
)
for (i in seq(from = 1, to = length(cell_maps), by = 1)) {
  clean_col_name <- paste('clean_yield_sum_', cell_maps[[i]][[year_index]], sep = "")
  yield_sum_map[, clean_col_name] <- cell_maps[[i]][[data_index]]$yield_sum
}

#------------------------------------------------------------------------------------------------



# spatial mean and variability
get_spatial_mean <- function(var) {
  return (mean(var))
}

get_spatial_variability <- function(var) {
  return (sd(var))
}

spatial_mean_variability_list <- list()
clean_columns <- names(yield_sum_map)[grep("^clean", names(yield_sum_map))]

for (i in 1:length(data)) {
  mean_ <- get_spatial_mean(yield_sum_map[, clean_columns[[i]]])
  sd_ <- get_spatial_variability(yield_sum_map[, clean_columns[[i]]])
  spatial_mean_variability_list[[i]] <- list(mean_, sd_)
}

# # temporal mean and variability
# get_temporal_mean_and_variability <- function(yearly_yields) {
#   temporal_mean <- vector('numeric', length = length(yearly_yields[[1]]))
#   temporal_sd <- vector('numeric', length = length(yearly_yields[[1]]))
#   for (i in 1:length(yearly_yields[[1]])) {
#     yearly_data_row <- vector('numeric', length = length(yearly_yields))
#     for (j in 1:length(yearly_yields)) {
#       yearly_data_row[j] <- yearly_yields[[j]][[i]]
#     }
#     temporal_mean[i] <- mean(yearly_data_row)
#     temporal_sd[i] <- sd(yearly_data_row)
#   }
#   return(data.frame(
#     temporal_mean = temporal_mean,
#     temporal_sd = temporal_sd
#   ))
# }
# 
# temp <- yield_sum_map[, unclean_columns]
# temporal_yield_list <- list()
# for (i in 1:length(unclean_columns)) {
#   col <- unclean_columns[[i]]
#   temp_yield <- yield_sum_map[, col]
#   temp_yield <- temp_yield - spatial_mean_variability_list[[i]][[1]]
#   temporal_yield_list[[i]] <- temp_yield
# }
# 
# temporal_analysis <- get_temporal_mean_and_variability(temporal_yield_list)
# temporal_analysis <- cbind(temporal_analysis, geometry = yield_sum_map$geometry)

############### outliers ############

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

#------------------------------------------------------------------------------------------------

# saving a file
get_directory_path <- function(field_name = "") {
  day <- format(Sys.time(), "%Y_%m_%d")
  path <- paste(day, "/", sep = "")
  if (field_name != "") {
    path <- paste(c(path, field_name, "/"), collapse = "")
  }
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

# # save the plots
# for (i in 1:length(data)) {
#   path <- get_directory_path(field_name)
#   dev <- "png"
#   filename <- get_filename_with_timestamp(paste("outliers_", data[[i]][[3]], sep = ""), dev)
#   plt <- plot_outliers(data[[i]][[5]], data[[i]][[3]])
#   ggsave(
#     filename = filename,
#     plot = plt,
#     device = dev,
#     dpi = 300,
#     path = path
#   )
# }


# compressed value z-statistic
yield_z_stat_val <- data.frame(
  geometry = yield_sum_map[,'geometry']
)
for (i in 1:length(data)) {
  col_name_ <- clean_columns[[i]]
  yield_ <- yield_sum_map[, col_name_]
  mean_ <- get_spatial_mean(yield_)
  sd_ <- get_spatial_variability(yield_)
  z_stat_ <- (yield_ - mean_) / sd_
  yield_z_stat_val[col_name_] <- z_stat_
}
yield_z_stat_val_sf <- st_as_sf(yield_z_stat_val)


# plotting z_stat for yield_sum
breaks <- c(-Inf, 0, Inf)
labels <- c('lessThan0', 'greaterEqualThan0')
plots <- list()
for (i in 1:length(data)) {
  col_name_ <- clean_columns[[i]]
  bin_col_name_ <- paste('bin_', col_name_, sep = '')
  yield_z_stat_val[bin_col_name_] <- cut(yield_z_stat_val[,col_name_], breaks = breaks, labels = labels, right = FALSE, include.lowest = TRUE)
}

# compute z's
compute_zs <- function(row) {
  unique_ <- unique(row)
  if (length(unique_) == 1) {
    result <- unique_
  } else {
    result <- 'Mixed'
  }
  return (result)
}
yield_z_stat_val['Zs'] <- apply(yield_z_stat_val[,bin_columns], 1, compute_zs)


# identifies the cell quartiles
get_del_val = function(zjs) {
  del_val <- abs(max(zjs) - min(zjs))
  return (del_val)
}

yield_z_stat_val <- as.data.frame(yield_z_stat_val)
yield_z_stat_val['del_val'] <- apply(yield_z_stat_val[, clean_columns], 1, get_del_val)

delta <- 1
plot_hs <- function(delta) {
  breaks = c(-Inf, delta, Inf)
  labels = c('lessEqualThanDelta', 'greaterThanDelta')
  bin_columns <- unlist(lapply(clean_columns, function (col) { paste('bin_', col, sep = '') }))
  yield_z_stat_val['Hs'] <- cut(yield_z_stat_val[,'del_val'], breaks = breaks, labels = labels, right = TRUE, include.lowest = TRUE)
  
  final_Hs <- function(row) {
    if (row['Hs'] == "lessEqualThanDelta") {
      if (row['Zs'] == "lessThan0") {
        return ('H3')
      } else if(row['Zs'] == "greaterEqualThan0") {
        return ('H1')
      } else {
        return ('H5')
      }
    } else {
      if (row['Zs'] == "lessThan0") {
        return ('H4')
      } else if(row['Zs'] == "greaterEqualThan0") {
        return ('H2')
      } else {
        return ('H6')
      }
    }
  }
  yield_z_stat_val['Final_Hs'] <- apply(yield_z_stat_val[, c('Hs', 'Zs')], 1, final_Hs)
  
  return_plot <- ggplot(st_as_sf(yield_z_stat_val[, c('Final_Hs', 'geometry')])) %>%
    
  
  return (return_plot)
}
plot_hs(2)


# # visualize the slider
# library(shiny)
# 
# ui <- fluidPage(
#   sliderInput("slider", "delta", min = min(yield_z_stat_val$del_val), max = max(yield_z_stat_val$del_val), value = min(yield_z_stat_val$del_val)),
#   plotOutput('plot')
# )
# 
# server <- function(input, output) {
#   observeEvent(input$slider, {
#     output$plot <- renderPlot({
#       plot_hs(input$slider)
#     })
#   })
# }
# 
# shinyApp(ui, server)
