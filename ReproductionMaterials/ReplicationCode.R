############################################
############################################
############################################
## Installation software
## Please follow the instructions on the GitHub site (https://github.com/soilsensingmonitoring/ml-hrsm_1.0) to install the "mlhrsm" package

install.packages(c('raster', 'rgee', 'sf', 'tidyverse', 'viridis', 'FedData', 'RColorBrewer', 'caret', 'chillR', 'leaflet', 'hydroGOF', 'quantregForest', 'randomForest', 'reshape2', 'rgdal', 'sp', 'lubridate', 'geojsonio', 'stars'))
install.packages(c("devtools","R.rsp"))

library(devtools)
install_version("rgdal", version = "1.6-7", repos = "http://cran.us.r-project.org")
library(rgdal)
options("rgdal_show_exportToProj4_warnings"="none")   ## Please run this line so that rgdal package can still be used without confusion with sf and terra packages

devtools::install_github("soilsensingmonitoring/ml-hrsm_1.0", build_vignettes=T)

## Set up a Google Earth Engine account (https://earthengine.google.com/signup/) 
## Install the gcloud CLI (available from https://dl.google.com/dl/cloudsdk/channels/rapid/GoogleCloudSDKInstaller.exe) 

library(mlhrsm)
ee_install()  ## Only needs to be run once

## Once the installatin is complete, run the following lines
library(mlhrsm)
ee_Initialize("the email address used to create the Google Earth Engine account", drive = T)




############################################
############################################
############################################
## Example code 
## the user can reproduce the results shown in the paper by running the following lines



## Split a large ROI to small subregions
split_region("Grant.shp", "Grant2021")   ## These files are available from https://github.com/soilsensingmonitoring/ml-hrsm_1.0/tree/main/data
sub_grant <- read_sf("Grant2021/sub_regions.shp")$geometry
plot(sub_grant)


## Case study 
## Downloading satellite data and apply the machine learning models to produce soil moisture maps
## VWC_map("WI.shp", "2020-06-15", "2020-09-15", 30, TRUE, "WI_region")
## Note: Downloading maps at the 30-m resolution could be time-consuming (~ 2.5 hours) 
## The user can run the example with 100-m resolution maps for illustration with similar results 
## This takes about 45 min mostly due to the processing of satellite data in Google Earth Engine Cloud 
VWC_map("WI.shp", "2020-06-15", "2020-09-15", 100, TRUE, "WI_region")

## Plot the soil moisture map at one depth on a specific date
plot_map("2020-08-01", 5, TRUE, project="WI_region")


## Spatial and temporal analysis examples
## Area-based/zonal functions
area_sum(5, project="WI_region")
head(read.csv("WI_region/VWC_5_ts_data.csv"))

## Display the zonal statistics maps
pixelwise_sum(100, project="WI_region")

## Re-calculate soil moisture from existing maps
aggregate_interval(5, "2020-06-15", "2020-08-15", frequency=7, project="WI_region")
head(read.csv("WI_region/VWC_aggregation/2020-06-15_2020-08-15/7days_5cm/7_day_mean.csv"))

## Plot the new zonal statistics map
plot_aggregated_VWC(5, "2020-06-15", "2020-08-15", date="2020-06-29", frequency=7, project="WI_region")



## Compare daily vs. weekly soil moisture maps
## define the dates to plot
dates <- c("2020-07-06", "2020-07-13", "2020-07-20", "2020-07-27")

## define the files (with specific VWC statistical characteristics) to plot
pattern <- paste(paste0("VWC_5_mean_", dates), sep="", collapse="|") 
maps <- stack(list.files("WI_region/VWC", pattern, full.names=T))
## set at same scale, daily maps
plot(maps, zlim=c(0.15, 0.35)) 

## Plot weekly maps
WD <- "WI_region/VWC_aggregation/2020-06-15_2020-08-15/7days_5cm"
pattern <- paste(paste0("mean_", dates), sep="", collapse="|")
maps <- stack(list.files(paste0(WD, "/aggregated_VWC"), pattern, full.names=T))
plot(maps, zlim=c(0.15, 0.3))
               

## Extract soil moisture data at individual sites
## Before running the following code, first copy and paste the following shapefiles ("WI_SM.shp", "WI_SM.dbf", "WI_SM.prj", "WI_SM.shx", "WI_SM.rda") into the Project's working directory, i.e., /WI_region/
point_extraction("WI_SM.shp", TRUE, project="WI_region")
head(read.csv("WI_region/VWC_point_data.csv"))

## Plot time series of soil moisture values (need to change the working directory to the Project's working directory)
setwd("WI_region/")
plots <- site_variation(5, TRUE, project="WI_point")
plots[[1]] # returns leaflet plot showing site locations
plots[[2]] # returns time series plot showing point variation over-time

## Plot times series of soil moisture with 90% CI
## Before running the following code, first change the working directory to the parent directory
setwd("../")
point_CI("Huges_VWC.csv", 5, project="WI_region")

## Plot the measured vs. predicted soil moisture at the selected sites
point_performance("Huges_VWC.csv", 5, TRUE, project="WI_region")
