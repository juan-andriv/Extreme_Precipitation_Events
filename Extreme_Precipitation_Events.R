### CHARACTERIZATION OF EXTREME PRECIPITATION EVENTS ###
### Juan Andrade Rivera
### linkedin.com/in/juan-andriv
### 16.03.2024

### input Data
### Daily precipitation data (e.g. MSWEP) in this example named "Cropped_daily_P"
### Shapefile for geographic extent (check http://tapiquen-sig.jimdo.com) in this example Europe was used

## Libraries
library(SpatIndex) # To calculate indices, download: https://github.com/obaezvil/SpatIndex
library(terra)
library(raster)
library(sf)
library(sp)
library(remotes)
library(rgdal)
library(ggplot2)

# set WD, create an "Output" folder inside to save all products
setwd("C:/Users/Juan/Documents/0_RStudio/Extreme_Precipitation_Events")


# Calculate desired index. SpatIndex include the following

# Abbreviation   |   Full Name
# -------------------------------------------
# Rx1day         |   Maximum 1-day precipitation amount
# Rx5day         |   Maximum 5-day precipitation amount
# SDII           |   Simple Daily Intensity Index
# R10mm          |   Number of days with precipitation >= 10 mm
# R20mm          |   Number of days with precipitation >= 20 mm
# Rnnmm          |   Number of days with precipitation >= nn mm (custom value)
# CDD            |   Consecutive dry days
# CWD            |   Consecutive wet days
# R95p           |   95th percentile of daily precipitation amounts
# R99p           |   99th percentile of daily precipitation amounts
# R95pTOT        |   95th percentile of total precipitation amounts
# R99pTOT        |   99th percentile of total precipitation amounts
# PRCPTOT        |   Total precipitation
# -------------------------------------------

# In this example the "Rx1day" index was used. 

Rx1day <- Rx1day("Cropped_daily_P", vct = NULL, start_date = 11, end_date = 20, init_month = NULL)
writeRaster(Rx1day, "Output/Rx1day.tif", overwrite = TRUE)
plot(Rx1day) # just to check

## Identify trends with the Non-parametric Mann-Kendall trend analysis
Rx1day <- "Output/Rx1day.tif"
Rx1day <- rast(Rx1day)
resultRx1day <- mk_spatial(Rx1day, vct = NULL, conf_level = 0.95, p_thres = 0.05) # run MK test

## save MK outputs
terra::writeCDF(resultRx1day$gridded, "Output/Rx1day_MK.nc", overwrite = TRUE) # save MK results as CDF
Rx1dayras <- terra::writeRaster(resultRx1day$gridded, "Output/Rx1day_MK.tif", overwrite = TRUE) # save MK results as tif
Rx1dayvec <- terra::writeVector(resultRx1day$sig_points, "Output/Rx1day_MK.shp", overwrite = TRUE) # save significant ponits as shp

## plot to check
plot(resultRx1day$gridded) # full MK results
plot(resultRx1day$sig_points) # just the significant points

## extract necessary infrmation 
# p value for the probability of occurence
p.value.Rx1day <- raster::subset(Rx1dayras, 2)
# Sen's slope for the trend (positive, negative, none)
s.slope.Rx1day <- raster::subset(Rx1dayras, 6)

## Identify significant changes
sig_changes.Rx1day <- p.value.Rx1day # add p data to new object
# keep only p values below 0.05, and normalize all to 1
sig_changes.Rx1day[sig_changes.Rx1day > 0.05] <- NA # all p values above 0.05 are NA
sig_changes.Rx1day[sig_changes.Rx1day > 0] <- 1 # reamaining p values are all 1
# keep only Sen's slope positive or negative trend
unit_changes.Rx1day <- s.slope.Rx1day / abs(s.slope.Rx1day) # reduces data only to directionality (+1 or -1)
# integrate
sig_changes.Rx1day  <- unit_changes.Rx1day * sig_changes.Rx1day

## Plotting ##
## Creating a map theme
theme_map <- function(){
  theme_bw() +
    theme(panel.grid = element_blank(),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = "cm"),
          plot.title = element_text(size = 20, vjust = 1, hjust = 0.5),
          legend.position = "right")
}

# Set geographic shapefile (Europe used as example)
shapefile.europe <- rgdal::readOGR("Europe/Europe.shp") # read shp
data.frame(shapefile.europe) # convert to dataframe
europe.vec <- terra::vect(shapefile.europe) # convert to vector

# Mapping trend Significance (p value < 0.05 means its statistically significant)
Rx1dayshp <- rgdal::readOGR("Output/Rx1day_MK.shp") # read MK shp
cropRx1day <- intersect(Rx1dayshp, shapefile.europe) # crop to world
data.frame(cropRx1day)
plot(cropRx1day) # just to check

# Visualization specs
mapRx1day <- ggplot(data.frame(cropRx1day)) +
  geom_point(data = NULL, mapping = aes(x = x, y = y, color = Significan),
    size = 0.25) +
  geom_polygon(
    data = shapefile.europe, aes(x = long, y = lat, group = group),
    color = "black",fill = NA) +
  scale_color_gradient(low = "red", high = "blue")+
  ggtitle("Trend analysis") + theme_map()

# Final view
plot(mapRx1day)
ggsave("Maps/mapRx1day_sig.png", width = 10, height = 8, dpi = 300)

# Mapping Trend Direction
sig_changes.Rx1day <- raster::mask(sig_changes.Rx1day, europe.vec)
png("Maps/mapRx1day_dir.png", width = 10, height = 8, units = "cm", res = 300) # to save as png
plot(sig_changes.Rx1day, col= c('red', 'blue'), main= "Trend direction") # plotting
lines(europe.vec) # adding Europe
dev.off() # closing 
