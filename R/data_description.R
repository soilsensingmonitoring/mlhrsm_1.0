#' Existing VWC data from all sites
#'
#' A dataset containing the existing VWC data obtained from sites all around the US,
#' originally used for the machine learning model.
#'
#' @format A data frame with 172887 rows and 41 variables:
#' \describe{
#'   \item{Longitude}{x coordinate of the site location}
#'   \item{Latitude}{y coordinate of the site location}
#'   \item{landcover}{the type of landcover at the site}
#'   \item{elevation}{Elevation data downloaded from \href{https://developers.google.com/earth-engine/datasets/catalog/USGS_3DEP_10m}{USGS DEM}}
#'   \item{slope}{Slope data downloaded from USGS DEM}
#'   \item{aspect}{Aspect data downloaded from USGS DEM}
#'   \item{hillshade}{Hillshade data downloaded from USGS DEM}
#'   \item{clay_5}{Clay content data at 0-5 cm downloaded from Polaris (Chaney et al. (2019))}
#'   \item{sand_5}{Sand content data at 0-5 cm downloaded from Polaris}
#'   \item{bd_5}{Bulk density data at 0-5 cm downloaded from Polaris}
#'   \item{clay_100}{Clay content data at 0-1 m downloaded from Polaris}
#'   \item{sand_100}{Sand content data at 0-1 m downloaded from Polaris}
#'   \item{bd_100}{Bulk density data at 0-1 m downloaded from Polaris}
#'   \item{ssm}{Surface soil moisture data downloaded from \href{https://developers.google.com/earth-engine/datasets/catalog/NASA_USDA_HSL_SMAP10KM_soil_moisture}{SMAP}}
#'   \item{susm}{Subsurface soil moisture data downloaded from SMAP}
#'   \item{vv}{Backscatter at VV (vertical transmit/vertical receive) mode data downloaded from \href{https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S1_GRD}{Sentinel-1}}
#'   \item{vh}{Backscatter at VH (vertical transmit/horizontal receive) mode data downloaded from Sentinel-1}
#'   \item{angle}{Incidence angle data downloaded from Sentinel-1}
#'   \item{LST}{Land surface temperature data downloaded from \href{https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD11A1}{MODIS}}
#'   \item{LS_B5}{Surface reflectance of band 5 (near-infrared) data downloaded from \href{https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C02_T1_L2}{Landsat 8}}
#'   \item{LS_B6}{Surface reflectance of band 6 (shortwave infrared-1) data downloaded from Landsat 8}
#'   \item{LS_B7}{Surface reflectance of band 7 (shortwave infrared-2) data downloaded from Landsat 8}
#'   \item{LS_B10}{Surface reflectance of band 10 (thermal) data downloaded from Landsat 8}
#'   \item{LS_NDVI}{Surface reflectance of NDVI (Normalized Difference Vegetation Index) data calculated from band 4}
#'   \item{LS_NDWI}{Surface reflectance of NDWI (Normalized Difference Water Index) data calculated from band 4}
#'   \item{landcover}{Land cover type data downloaded from \href{https://developers.google.com/earth-engine/datasets/catalog/USGS_NLCD_RELEASES_2016_REL}{NCLD}}
#'   ...
#' }
"all"



#' ML model for surface VWC prediction
#'
#' A dataset containing the soil surface level ML model based on the existing
#' VWC data obtained from sites all around the US.
#'
#' @format A list with 23 objects
"model_surface"



#' ML model for rootzone VWC prediction based on existing data
#'
#' A dataset containing the soil root zone ML model based on the existing
#' VWC data obtained from sites all around the US.
#'
#' @format A list with 23 objects
"model_rz"



#' Polygon for testing package regional mapping
#'
#' A dataset of a cropland in Rock County, WI.
#' Note that the data can't be directly input into the functions. To use the data,
#' users should first save it into a shapefile for the mapping function to access.
#'
#' @usage data(WI)
#'
#' @format A polygon with area of 70-ha
"WI"



#' Points for testing package point mapping
#'
#' A dataset consisting of 12 points within the cropland in Rock County, WI.
#' Note that the data can't be directly input into the functions. To use the data,
#' users should first save it into a shapefile for the mapping function to access.
#'
#' @usage data(WI_SM)
#'
#' @format 12 points
#'
#' @references TEROS 12, Meter Group, Inc.
"WI_SM"
