data(all, envir=environment())
data(model_surface, envir=environment())
data(model_rz, envir=environment())

#' Regional mapping for small areas
#'
#' The major function that predicts the VWC across a region. It mostly works for small regions (preferably less than 1000 km2)
#' and/or low resolutions. If the downloaded maps are too large that GEE automatically splits them, the function will ask users
#' to consider splitting the region and mosaicking it afterwards. (direct to function: splitting())
#'
#' @param area character. Gives the name of the file that stores the study area (in shapefile format). It can be a file within the working directory or a detailed path specifying the location of the file.
#' @param start_date character. Defines the start date of the prediction period. The start date should be no earlier than 2016-01-01.
#' @param end_date character. Defines the end date of the prediction period.
#' @param resolution numeric. Defines the resolution of mapping region. The resolution should be smaller than or equal to 500m.
#' @param percentile logical. If TRUE, the upper and lower bounds of the 90\% Confidence Interval of the prediction will also be mapped.
#' @param project character. Defines the name of the folder to be created and save the maps. If not defined, the maps will be saved in the working directory.
#'
#' @return None
#'
#' @import FedData rgee randomForest quantregForest caret sp sf raster rgdal leaflet tidyverse reshape2
#'
#' @export
#'
#' @examples
#' \dontrun{ VWC_map("WI.shp", "2020-06-01", "2020-09-01", 30, TRUE, "WI_region")}
VWC_map <- function(area, start_date, end_date, resolution, percentile=FALSE,
                       project=NULL, large_area=FALSE){

  setwd(getwd())
  WD <- getwd()

  #####################################################################
  #####################################################################
  #####################################################################
  ##### Note: You can skip Part 1 by loading the pre-saved model into R and update the working directory to the current one
  # load("Existing_model.RData")
  # load("Existing_model_all_sites.rdata")
  setwd(getwd())
  WD <- getwd()
  if(!is.null(project)){
    if (file.exists(project)){
      return("Duplicated project folder name.")
    } else {
      dir.create(project)
      WD <- paste0(WD, "/", project)
    }
  }

  bunk <- function(i, message, pb, j){
    cat("0\14")
    writeLines(c(paste0("Tasks done: ", i, "/10"),message))
    setTxtProgressBar(pb, j)
  }
  #####################################################################
  #####################################################################

  ############################ Part 2  Mapping soil moisture at surface (0-5 cm) and rootzone (0-1 m)
  {

    roi_new <-  st_transform(st_read(area),crs=4326)

    print(paste0("site = ", area))

    boundary <- as(roi_new , Class = "Spatial")

    roi <- sf_as_ee(roi_new$geometry)

    area <- area(boundary)/1000000 # in km2

    if (!large_area) {
      if (area > 1000 & resolution <=50)
        return("ROI is too large. Please split the region of interest. Set large_area=TRUE to bypass the message.")
    }

    ## Define the spatial resolution (meters) of final soil moisture maps
    resolution <- resolution

    ## Define starting and end dates of soil moisture maps
    doy_frame <- c(as.Date(start_date), as.Date(end_date))


    ############################ Part 2  Download remote sensing data and land surface parameters from Google Earth Engine & Cloud Storage

    ## Download remote sensing data and land surface parameters from GEE and Cloud Storage
    {
      ## DEM data and terrain parameters
      dem <- ee$Image("USGS/SRTMGL1_003")
      terrain <- ee$Terrain$products(dem)$toFloat()$reproject(crs = "EPSG:4326", scale = resolution)$clip(roi)

      ## Land cover from USGS NLCD 2016
      LC = ee$ImageCollection("USGS/NLCD_RELEASES/2016_REL")$filterDate('2016-01-01', '2016-12-31')$select('landcover')$first()$rename("LC")$reproject(crs = "EPSG:4326", scale = resolution)$clip(roi)

      cat("\014")
      writeLines(c("Tasks done: 0/10","Downloading static covariates..."))
      pb <- txtProgressBar(min=0, max=4, style=3)
      ## Start downloading static covariates within the ROI
      terrain_raster <- ee_as_raster(image = terrain, region = roi, dsn = file.path(WD, "terrain"), scale = resolution, quiet =T, maxPixels = 1e+13)
      if (length(list.files(WD, "terrain"))>1) return("File size is too large. Please split the region of interest.")
      bunk(0,"Downloading static covariates...",pb,1)
      LC_raster <- ee_as_raster(image = LC, region = roi, dsn = file.path(WD, "LC"), scale = resolution, quiet =T, maxPixels = 1e+13)
      if (length(list.files(WD, "LC"))>1) return("File size is too large. Please split the region of interest.")
      bunk(0,"Downloading static covariates...",pb,2)

      ## Input and plot static covariates
      terrain_raster <- stack(paste0(WD,"/terrain.tif"))
      names_terrain <- names(terrain_raster)
      LC_raster <- raster(paste0(WD, "/LC.tif"))


      clay_5 <- ee$Image("users/xuehaiwuya8/ML-HRSM/polaris_clay_5")
      sand_5 <- ee$Image("users/xuehaiwuya8/ML-HRSM/polaris_sand_5")
      bd_5 <- ee$Image("users/xuehaiwuya8/ML-HRSM/polaris_bd_5")
      clay_5_15 <- ee$Image("users/xuehaiwuya8/ML-HRSM/polaris_clay_15")
      sand_5_15 <- ee$Image("users/xuehaiwuya8/ML-HRSM/polaris_sand_15")
      bd_5_15 <- ee$Image("users/xuehaiwuya8/ML-HRSM/polaris_bd_15")
      clay_15_30 <- ee$Image("users/xuehaiwuya8/ML-HRSM/polaris_clay_30")
      sand_15_30 <- ee$Image("users/xuehaiwuya8/ML-HRSM/polaris_sand_30")
      bd_15_30 <- ee$Image("users/xuehaiwuya8/ML-HRSM/polaris_bd_30")
      clay_30_60 <- ee$Image("users/xuehaiwuya8/ML-HRSM/polaris_clay_60")
      sand_30_60 <- ee$Image("users/xuehaiwuya8/ML-HRSM/polaris_sand_60")
      bd_30_60 <- ee$Image("users/xuehaiwuya8/ML-HRSM/polaris_bd_60")
      clay_60_100 <- ee$Image("users/xuehaiwuya8/ML-HRSM/polaris_clay_100")
      sand_60_100 <- ee$Image("users/xuehaiwuya8/ML-HRSM/polaris_sand_100")
      bd_60_100 <- ee$Image("users/xuehaiwuya8/ML-HRSM/polaris_bd_100")

      clay_100 <- (clay_5$multiply(5)$add(clay_5_15$multiply(10))$add(clay_15_30$multiply(15))$add(clay_30_60$multiply(30))$add(clay_60_100$multiply(40)))$multiply(0.01)
      sand_100 <- (sand_5$multiply(5)$add(sand_5_15$multiply(10))$add(sand_15_30$multiply(15))$add(sand_30_60$multiply(30))$add(sand_60_100$multiply(40)))$multiply(0.01)
      bd_100 <- (bd_5$multiply(5)$add(bd_5_15$multiply(10))$add(bd_15_30$multiply(15))$add(bd_30_60$multiply(30))$add(bd_60_100$multiply(40)))$multiply(0.01)

      soil_5 <- clay_5$rename('clay_5')$addBands(sand_5$rename('sand_5'))$addBands(bd_5$rename('bd_5'))$reproject(crs = "EPSG:4326", scale = resolution)$clip(roi)
      soil_100 <- clay_100$rename('clay_100')$addBands(sand_100$rename('sand_100'))$addBands(bd_100$rename('bd_100'))$reproject(crs = "EPSG:4326", scale = resolution)$clip(roi)

      soil_5_raster <- ee_as_raster(image = soil_5, region = roi, dsn = file.path(WD, "soil_5"), scale = resolution, quiet =T, maxPixels = 1e+13)
      if (length(list.files(WD, "soil_5"))>1) return("File size is too large. Please split the region of interest.")
      bunk(0,"Downloading static covariates...",pb,3)
      soil_100_raster <- ee_as_raster(image = soil_100, region = roi, dsn = file.path(WD, "soil_100"), scale = resolution, quiet =T, maxPixels = 1e+13)
      if (length(list.files(WD, "soil_100"))>1) return("File size is too large. Please split the region of interest.")
      bunk(0,"Downloading static covariates...",pb,4)
      close(pb)

      ## Input and plot static covariates
      soil_5_raster <- stack(paste0(WD, "/soil_5.tif"))
      soil_100_raster <- stack(paste0(WD, "/soil_100.tif"))

      names(soil_5_raster) <- c("clay_5", "sand_5", "bd_5")
      names(soil_100_raster) <- c("clay_100", "sand_100", "bd_100")


      constant_raster <- stack(terrain_raster, soil_5_raster, soil_100_raster)
    }


    ## Temporal gap-filling function from https://gist.github.com/johnbaums/10465462
    {
      interpolateTemporal <- function(s, xin, xout, outdir, prefix, progress,
                                      writechange, returnstack, ...) {

        if(missing(outdir)) stop('Please specify outdir')
        if(missing(prefix)) stop('Please specify prefix')
        if(nlayers(s) != length(xin)) stop('Length of xin must equal the number of layers in s.')
        if(nlayers(s) < 2) stop('stack s must have at least 2 layers.')
        if(!all(findInterval(xout, range(xin), rightmost.closed=TRUE) == 1)) {
          if(any(xout < min(xin))) {
            stop('This function does not extrapolate backwards (i.e. below the earliest element in xin). All elements of xout must be greater that min(xin).')
          } else {
            warning('Some values of xout require extrapolation beyond the range of xin.\nThe rate of change for extrapolation is assumed to be equal to that for the period between the last and second-to-last elements of xin (after temporal ordering).')
          }
        }
        outdir <- normalizePath(sub('/$|\\\\$', '', outdir), winslash='/',
                                mustWork=FALSE)
        if(!file.exists(outdir)) dir.create(outdir, recursive=TRUE)
        xout <- unique(xout)
        if(is.unsorted(xin)) {
          s <- s[[order(xin)]]
          xin <- sort(xin)
        }
        len <- diff(xin)
        base <- findInterval(xout, xin)
        lower <- unique(base[base < nlayers(s)])
        s.change <- stack(sapply(if(length(lower) > 0) lower else nlayers(s) - 1,
                                 function(x) {
                                   message(sprintf('Calculating change grid for %s to %s.', xin[x], xin[x+1]))
                                   overlay(s[[x]], s[[x+1]], fun=function(x1, x2) (x2-x1)/len[x],
                                           filename=ifelse(writechange,
                                                           file.path(outdir, sprintf('changegrid_%s_%s', xin[x], xin[x+1])),
                                                           ''), recycle=FALSE, format='GTiff', ...)
                                 }))

        multi <- xout - xin[base]
        chg.ind <- ifelse(base > nlayers(s.change), nlayers(s.change), base)
        message('Calculating grids for years specified in xout...')
        if(progress) pb <- txtProgressBar(0, length(xout), style=3)
        invisible(sapply(seq_along(xout), function(x) {
          out.rast <- if(xout[x] %in% xin) {
            s[[base[x]]]
          } else {
            overlay(s[[base[x]]], s.change[[chg.ind[x]]],
                    fun=function(x1, x2) x1 + (x2*multi[x]))
          }
          writeRaster(out.rast,
                      filename=file.path(outdir, sprintf('%s_%s', prefix, xout[x])),
                      format='GTiff', ...)
          if(progress) setTxtProgressBar(pb, x)
        }))
        if(isTRUE(returnstack)) {
          f <- file.path(outdir, paste0(prefix, '_', xout, '.tif'))
          return(stack(f[order(as.numeric(sub('.*_(\\d+)\\.tif$', '\\1', f)))]))
        }
      }
    }

    ## Masking functions
    {
      ## Function for clip image collections from GEE
      clip_fun = function(img){
        img <- img$reproject(crs = "EPSG:4326", scale = resolution)$clip(roi)
        return(img)
      }

      ## Function for pre-processing sentinel-1 data
      preprocess_vv = function(img){
        vv_masked <- img$mask(img$gt(-20)$And(img$lt(-5)))
        img <- vv_masked$convolve(ee$Kernel$gaussian(3))
        return(img)
      }

      preprocess_vh = function(img){
        vh_masked <- img$mask(img$gt(-30)$And(img$lt(-10)))
        img <- vh_masked$convolve(ee$Kernel$gaussian(3))
        return(img)
      }

      ## Cloud masking function for Landsat-8_L2
      maskL8sr <- function(image) {
        cloudShadowBitMask <- bitwShiftL(1, 4)
        cloudsBitMask <- bitwShiftL(1, 3)
        snowBitMask <- bitwShiftL(1, 5)
        qa <- image$select('QA_PIXEL')
        mask <- qa$bitwiseAnd(cloudShadowBitMask)$eq(0)$And(qa$bitwiseAnd(cloudsBitMask)$eq(0))$And(qa$bitwiseAnd(snowBitMask)$eq(0))

        image$updateMask(mask)$copyProperties(image, list("system:time_start"))
      }

      applyScaleFactors <- function(image) {
        opticalBands <- image$select('SR_B.')$multiply(0.0000275)$add(-0.2)
        thermalBands <- image$select('ST_B.*')$multiply(0.00341802)$add(149.0)

        image$addBands(opticalBands, NULL, TRUE)$addBands(thermalBands, NULL, TRUE)
      }

      ## Cloud masking function for MODIS
      maskMODIS <- function(image) {
        MandatoryBitMask <- bitwShiftL(1, 1)
        DataBitMask <- bitwShiftL(1, 3)
        qa <- image$select('QC_Day')
        mask <- qa$bitwiseAnd(MandatoryBitMask)$eq(0)$And(qa$bitwiseAnd(DataBitMask)$eq(0))

        image$updateMask(mask)$copyProperties(image, list("system:time_start"))
      }


    }

    ## Download dynamic covariates
    {


      dir.create(paste0(WD, "/covariates_temp"))

      ## NASA_USDA Enhanced SMAP data - image collection from GEE
      ssm = ee$ImageCollection("NASA_USDA/HSL/SMAP10KM_soil_moisture")$select('ssm')$filterBounds(roi)$filterDate(as.character(doy_frame[1] - 6), as.character(doy_frame[2] + 6))
      susm = ee$ImageCollection("NASA_USDA/HSL/SMAP10KM_soil_moisture")$select('susm')$filterBounds(roi)$filterDate(as.character(doy_frame[1] - 6), as.character(doy_frame[2] + 6))
      ssm = ssm$map(clip_fun)$toBands()
      susm = susm$map(clip_fun)$toBands()

      ## MODIS Land surface temperature data - image collection from GEE
      LST = ee$ImageCollection("MODIS/006/MOD11A1")$filterBounds(roi)$filterDate(as.character(doy_frame[1] - 24), as.character(doy_frame[2] + 24))$map(maskMODIS)$select('LST_Day_1km')
      LST = LST$map(clip_fun)$toBands()

      ## Mosaic S1 and LS-8 as they are not global datasets
      dates_dynamics <- seq(doy_frame[1], doy_frame[2], 1 )

      ## Sentinel-1 data - image collection from GEE
      dates_length_vv <- ceiling(length(dates_dynamics)/6) + 48/6
      dates_dynamics_select <- dates_dynamics[1] - 24

      vv = ee$ImageCollection('COPERNICUS/S1_GRD')$filter(ee$Filter$eq('instrumentMode', 'IW'))$filter(ee$Filter$listContains('transmitterReceiverPolarisation', 'VV'))$select('VV')$filterDate(as.character(dates_dynamics_select), as.character(dates_dynamics_select + 6))$sort('SLC_Processing_start')$filterBounds(roi)$map(clip_fun)$map(preprocess_vv)$mean()$rename(as.character(dates_dynamics_select))
      vh = ee$ImageCollection('COPERNICUS/S1_GRD')$filter(ee$Filter$eq('instrumentMode', 'IW'))$filter(ee$Filter$listContains('transmitterReceiverPolarisation', 'VH'))$select('VH')$filterDate(as.character(dates_dynamics_select), as.character(dates_dynamics_select + 6))$sort('SLC_Processing_start')$filterBounds(roi)$map(clip_fun)$map(preprocess_vh)$mean()$rename(as.character(dates_dynamics_select))
      angle = ee$ImageCollection('COPERNICUS/S1_GRD')$filter(ee$Filter$eq('instrumentMode', 'IW'))$filter(ee$Filter$listContains('transmitterReceiverPolarisation', 'VV'))$select('angle')$filterDate(as.character(dates_dynamics_select), as.character(dates_dynamics_select + 6))$sort('SLC_Processing_start')$filterBounds(roi)$map(clip_fun)$mean()$rename(as.character(dates_dynamics_select))

      for (l in 2:dates_length_vv)
      {

        dates_dynamics_select <- dates_dynamics[1] - 24 +(l-1)*6
        ## Sentinel-1 data - image collection from GEE
        vv_temp = ee$ImageCollection('COPERNICUS/S1_GRD')$filter(ee$Filter$eq('instrumentMode', 'IW'))$filter(ee$Filter$listContains('transmitterReceiverPolarisation', 'VV'))$select('VV')$filterDate(as.character(dates_dynamics_select), as.character(dates_dynamics_select + 6))$sort('SLC_Processing_start')$filterBounds(roi)$map(clip_fun)$map(preprocess_vv)$mean()
        vh_temp = ee$ImageCollection('COPERNICUS/S1_GRD')$filter(ee$Filter$eq('instrumentMode', 'IW'))$filter(ee$Filter$listContains('transmitterReceiverPolarisation', 'VH'))$select('VH')$filterDate(as.character(dates_dynamics_select), as.character(dates_dynamics_select + 6))$sort('SLC_Processing_start')$filterBounds(roi)$map(clip_fun)$map(preprocess_vh)$mean()
        angle_temp = ee$ImageCollection('COPERNICUS/S1_GRD')$filter(ee$Filter$eq('instrumentMode', 'IW'))$filter(ee$Filter$listContains('transmitterReceiverPolarisation', 'VV'))$select('angle')$filterDate(as.character(dates_dynamics_select), as.character(dates_dynamics_select + 6))$sort('SLC_Processing_start')$filterBounds(roi)$map(clip_fun)$mean()

        vv <- vv$addBands(vv_temp$rename(as.character(dates_dynamics_select)))
        vh <- vh$addBands(vh_temp$rename(as.character(dates_dynamics_select)))
        angle <- angle$addBands(angle_temp$rename(as.character(dates_dynamics_select)))

      }

      ## Landsat-8 bands - image collection from GEE

      dates_length_LS <- ceiling(length(dates_dynamics)/16) + 192*2/16
      dates_dynamics_select <- dates_dynamics[1] - 192

      LS_B4 = ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")$filterBounds(roi)$filterDate(as.character(dates_dynamics_select), as.character(dates_dynamics_select + 16))$map(maskL8sr)$map(applyScaleFactors)$select('SR_B4')$map(clip_fun)$mean()$rename(as.character(dates_dynamics_select))
      LS_B5 = ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")$filterBounds(roi)$filterDate(as.character(dates_dynamics_select), as.character(dates_dynamics_select + 16))$map(maskL8sr)$map(applyScaleFactors)$select('SR_B5')$map(clip_fun)$mean()$rename(as.character(dates_dynamics_select))
      LS_B6 = ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")$filterBounds(roi)$filterDate(as.character(dates_dynamics_select), as.character(dates_dynamics_select + 16))$map(maskL8sr)$map(applyScaleFactors)$select('SR_B6')$map(clip_fun)$mean()$rename(as.character(dates_dynamics_select))
      LS_B7 = ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")$filterBounds(roi)$filterDate(as.character(dates_dynamics_select), as.character(dates_dynamics_select + 16))$map(maskL8sr)$map(applyScaleFactors)$select('SR_B7')$map(clip_fun)$mean()$rename(as.character(dates_dynamics_select))
      LS_B10 = ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")$filterBounds(roi)$filterDate(as.character(dates_dynamics_select), as.character(dates_dynamics_select + 16))$map(maskL8sr)$map(applyScaleFactors)$select('ST_B10')$map(clip_fun)$mean()$rename(as.character(dates_dynamics_select))

      cat("\014")
      writeLines(c("Tasks done: 1/10","Setting bands..."))
      pb <- txtProgressBar(min=1, max=dates_length_LS, style=3)
      for (l in 2:dates_length_LS)
      {
        dates_dynamics_select <- dates_dynamics[1] - 192 +(l-1)*16

        ## Landsat-8 bands - image collection from GEE
        LS_B4_temp = ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")$filterBounds(roi)$filterDate(as.character(dates_dynamics_select), as.character(dates_dynamics_select + 16))$map(maskL8sr)$map(applyScaleFactors)$select('SR_B4')$map(clip_fun)$mean()
        LS_B5_temp = ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")$filterBounds(roi)$filterDate(as.character(dates_dynamics_select), as.character(dates_dynamics_select + 16))$map(maskL8sr)$map(applyScaleFactors)$select('SR_B5')$map(clip_fun)$mean()
        LS_B6_temp = ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")$filterBounds(roi)$filterDate(as.character(dates_dynamics_select), as.character(dates_dynamics_select + 16))$map(maskL8sr)$map(applyScaleFactors)$select('SR_B6')$map(clip_fun)$mean()
        LS_B7_temp = ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")$filterBounds(roi)$filterDate(as.character(dates_dynamics_select), as.character(dates_dynamics_select + 16))$map(maskL8sr)$map(applyScaleFactors)$select('SR_B7')$map(clip_fun)$mean()
        LS_B10_temp = ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")$filterBounds(roi)$filterDate(as.character(dates_dynamics_select), as.character(dates_dynamics_select + 16))$map(maskL8sr)$map(applyScaleFactors)$select('ST_B10')$map(clip_fun)$mean()


        skip_to_next <- FALSE

        tryCatch(ee_print(LS_B4_temp), error = function(e) { skip_to_next <<- TRUE})

        if(skip_to_next) {
          bunk(1,"Setting bands...",pb,l)
          next
        }

        if(!skip_to_next)
        {
          cat("\014")
          LS_B4 <- LS_B4$addBands(LS_B4_temp$rename(as.character(dates_dynamics_select)))
          LS_B5 <- LS_B5$addBands(LS_B5_temp$rename(as.character(dates_dynamics_select)))
          LS_B6 <- LS_B6$addBands(LS_B6_temp$rename(as.character(dates_dynamics_select)))
          LS_B7 <- LS_B7$addBands(LS_B7_temp$rename(as.character(dates_dynamics_select)))
          LS_B10 <- LS_B10$addBands(LS_B10_temp$rename(as.character(dates_dynamics_select)))
          bunk(1,"Setting bands...",pb,l)
        }
      }
      close(pb)


      ##################################################################################
      ##################################################################################
      ## Start downloading from GEE and Google Cloud Storage and load images to R memory
      ##################################################################################
      ##################################################################################

      ## SMAP
      cat("\014")
      writeLines(c("Tasks done: 2/10","Downloading map..."))
      pb <- txtProgressBar(min=0, max=3, style=3)
      ee_as_raster(image = ssm, region = roi, dsn = file.path(paste0(WD, "/covariates_temp"), "ssm_NASA_raw"), scale = resolution, quiet =T, maxPixels = 1e+13)
      if (length(list.files(paste0(WD, "/covariates_temp"),"ssm"))>1) return("File size is too large. Please split the region of interest.")
      ssm_raster <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "ssm_NASA_raw", full.names = T))
      bunk(2,"Downloading map...",pb,1)
      ee_as_raster(image = susm, region = roi, dsn = file.path(paste0(WD, "/covariates_temp"), "susm_NASA_raw"), scale = resolution, quiet =T, maxPixels = 1e+13)
      if (length(list.files(paste0(WD, "/covariates_temp"),"susm"))>1) return("File size is too large. Please split the region of interest.")
      susm_raster <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "susm_NASA_raw", full.names = T))
      bunk(2,"Downloading map...",pb,2)

      SMAP_dates <- substr(names(ssm_raster), 18, 25)
      SMAP_unique_dates <- sort(unique(SMAP_dates))
      SMAP_index <- NULL
      SMAP_list <- 1:length(SMAP_dates)
      for (i in 1:length(SMAP_unique_dates))
      {
        temp <-   SMAP_list[SMAP_dates %in%  SMAP_unique_dates[i]][1]
        SMAP_index <- c(SMAP_index, temp)
      }

      ## MODIS
      ee_as_raster(image = LST, region = roi, dsn = file.path(paste0(WD, "/covariates_temp"), "LST_NASA_raw"), scale = resolution, quiet =T, maxPixels = 1e+13)
      if (length(list.files(paste0(WD, "/covariates_temp"),"LST"))>1) return("File size is too large. Please split the region of interest.")
      LST_raster <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "LST_NASA_raw", full.names = T))
      bunk(2,"Downloading map...",pb,3)
      close(pb)

      # In case there are large gaps between subregion polygons
      # {
      #   test_days <- 30
      #   if (nlayers(LST_raster)<30) test_days <- nlayers(LST_raster)
      #   if (sum(sum(values(LST_raster[[1:test_days]])==0, na.rm=T)/length(LST_raster[[1]]) < 4/5 ) == 0) {
      #     for (i in seq(nlayers(LST_raster))) {
      #       cat("\014")
      #       writeLines(paste0("Masking LST raster ", i, "/", nlayers(LST_raster)))
      #       LST_raster[[i]] <- mask(LST_raster[[i]], boundary)
      #     }
      #     writeRaster(LST_raster, paste0(WD, "/covariates_temp/LST_NASA_raw.tif"), overwrite=T)
      #   }
      # }


      ## Remove MODIS images that have a large cloud cover
      {

        ## Get the index of MODIS images that are not masked by cloud mask (0 values) and  with less than 33% cloud cover
        mask_dates_MODIS <- NULL
        cat("\014")
        writeLines(c("Tasks done: 3/10","Removing cloud cover..."))
        pb <- txtProgressBar(min=0, max=dim(LST_raster)[3], style=3)
        for (i in 1:dim(LST_raster)[3])
        {
          if(summary(values(LST_raster[[i]]))[3]!=0 & (sum(values(LST_raster[[i]])==0)/length(LST_raster[[1]]) < 4/5 ) )
            mask_dates_MODIS <- c(mask_dates_MODIS, i)
          bunk(3,"Removing cloud cover...",pb,i)
        }
        close(pb)


      }

      values(LST_raster)[values(LST_raster) == 0] = NA
      LST_dates <- unlist(strsplit(ee_print(LST, quiet=T)$img_bands_names, " "))
      LST_dates <- gsub("(.*)(\\d{4}_\\d{2}_\\d{2})(.*)", "\\2", LST_dates)
      # LST_dates <- paste0(substr(names(LST_raster), 2, 5), substr(names(LST_raster), 7, 8), substr(names(LST_raster), 10, 11))


      ## Sentinel-1
      vv_na <- FALSE
      tryCatch(ee_print(vv), error = function(e) { vv_na <<- TRUE})
      cat("\014")
      writeLines(c("Tasks done: 4/10","Downloading map..."))
      pb <- txtProgressBar(min=0, max=8, style=3)
      if(!vv_na)
      {
        ee_as_raster(image = vv, region = roi, dsn = file.path(paste0(WD, "/covariates_temp"), "vv_S1_raw"), scale = resolution, quiet =T, maxPixels = 1e+13)
        if (length(list.files(paste0(WD, "/covariates_temp"),"vv"))>1) return("File size is too large. Please split the region of interest.")
        vv_raster <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "vv_S1_raw", full.names = T))
        bunk(4,"Downloading map...",pb,1)
        ee_as_raster(image = vh, region = roi, dsn = file.path(paste0(WD, "/covariates_temp"), "vh_S1_raw"), scale = resolution, quiet =T, maxPixels = 1e+13)
        if (length(list.files(paste0(WD, "/covariates_temp"),"vh"))>1) return("File size is too large. Please split the region of interest.")
        vh_raster <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "vh_S1_raw", full.names = T))
        bunk(4,"Downloading map...",pb,2)
        ee_as_raster(image = angle, region = roi, dsn = file.path(paste0(WD, "/covariates_temp"), "angle_S1_raw"), scale = resolution, quiet =T, maxPixels = 1e+13)
        if (length(list.files(paste0(WD, "/covariates_temp"),"angle"))>1) return("File size is too large. Please split the region of interest.")
        angle_raster <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "angle_S1_raw", full.names = T))
        bunk(4,"Downloading map...",pb,3)

        S1_dates <- paste0(substr(names(vv_raster), 2, 5), substr(names(vv_raster), 7, 8), substr(names(vv_raster), 10, 11))
      }

      ## Landsat-8
      ee_as_raster(image = LS_B4, region = roi, dsn = file.path(paste0(WD, "/covariates_temp"), "LS_B4_NASA_raw"), scale = resolution, quiet =T, maxPixels = 1e+13)
      if (length(list.files(paste0(WD, "/covariates_temp"),"LS_B4"))>1) return("File size is too large. Please split the region of interest.")
      bunk(4,"Downloading map...",pb,4)
      ee_as_raster(image = LS_B5, region = roi, dsn = file.path(paste0(WD, "/covariates_temp"), "LS_B5_NASA_raw"), scale = resolution, quiet =T, maxPixels = 1e+13)
      if (length(list.files(paste0(WD, "/covariates_temp"),"LS_B5"))>1) return("File size is too large. Please split the region of interest.")
      bunk(4,"Downloading map...",pb,5)
      ee_as_raster(image = LS_B6, region = roi, dsn = file.path(paste0(WD, "/covariates_temp"), "LS_B6_NASA_raw"), scale = resolution, quiet =T, maxPixels = 1e+13)
      if (length(list.files(paste0(WD, "/covariates_temp"),"LS_B6"))>1) return("File size is too large. Please split the region of interest.")
      bunk(4,"Downloading map...",pb,6)
      ee_as_raster(image = LS_B7, region = roi, dsn = file.path(paste0(WD, "/covariates_temp"), "LS_B7_NASA_raw"), scale = resolution, quiet =T, maxPixels = 1e+13)
      if (length(list.files(paste0(WD, "/covariates_temp"),"LS_B7"))>1) return("File size is too large. Please split the region of interest.")
      bunk(4,"Downloading map...",pb,7)
      ee_as_raster(image = LS_B10, region = roi, dsn = file.path(paste0(WD, "/covariates_temp"), "LS_B10_NASA_raw"), scale = resolution, quiet =T, maxPixels = 1e+13)
      if (length(list.files(paste0(WD, "/covariates_temp"),"LS_B10"))>1) return("File size is too large. Please split the region of interest.")
      bunk(4,"Downloading map...",pb,8)
      close(pb)

      ## Remove Landsat-8 images that have repeated dates and a large cloud cover
      {

        ## Load the images
        LS_B4_raster  <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "LS_B4_NASA_raw", full.names = T))
        LS_B5_raster  <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "LS_B5_NASA_raw", full.names = T))
        LS_B6_raster  <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "LS_B6_NASA_raw", full.names = T))
        LS_B7_raster  <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "LS_B7_NASA_raw", full.names = T))
        LS_B10_raster <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "LS_B10_NASA_raw", full.names = T))

        LS_dates <- paste0(substr(names(LS_B4_raster), 2, 5), substr(names(LS_B4_raster), 7, 8), substr(names(LS_B4_raster), 10, 11))


        names(LS_B4_raster) <- names(LS_B5_raster) <- names(LS_B6_raster) <- names(LS_B7_raster) <- names(LS_B10_raster) <- LS_dates


        ## Get the index of landsat-8 images with less than 80% cloud cover
        mask_dates_LS <- NULL
        cat("\014")
        writeLines(c("Tasks done: 5/10","Removibg cloud cover..."))
        pb <- txtProgressBar(min=0, max=dim(LS_B4_raster)[3], style=3)
        for (i in 1:dim(LS_B4_raster)[3])
        {
          if(sum(is.na(values(LS_B4_raster[[i]])))/length(LS_B4_raster[[1]]) < 0.8  & summary(values(LS_B4_raster[[i]]))[3] < 0.6 )
            mask_dates_LS <- c(mask_dates_LS, i)
          bunk(5,"Removing cloud cover...",pb,i)
        }
        close(pb)


      }

      values(LS_B4_raster)[values(LS_B4_raster) <= 0] = NA
      values(LS_B5_raster)[values(LS_B5_raster) <= 0] = NA
      values(LS_B6_raster)[values(LS_B6_raster) <= 0] = NA
      values(LS_B7_raster)[values(LS_B7_raster) <= 0] = NA
      values(LS_B10_raster)[values(LS_B10_raster) <= 0] = NA
    }

    ## Temporal gap-filling dynamic covariates using pixel-wise linear interpolation and load the gap-filled images into R
    {
      cat("\014")
      writeLines(c("Tasks done: 6/10","Temporal gap-filling..."))
      pb <- txtProgressBar(min=0, max=11, style=3)
      ## SMAP
      interpolateTemporal(s = ssm_raster[[SMAP_index]], xin = as.numeric(as.Date(SMAP_dates[SMAP_index], "%Y%m%d")), xout = as.numeric(seq.Date(doy_frame[1],doy_frame[2], 1)),
                          outdir = paste0(WD, "/covariates_temp"), prefix = "ssm_doy", progress=F, writechange=F, returnstack=F, overwrite=TRUE)
      bunk(6,"Temporal gap-filling...",pb,1)
      interpolateTemporal(s = susm_raster[[SMAP_index]], xin = as.numeric(as.Date(SMAP_dates[SMAP_index], "%Y%m%d")), xout = as.numeric(seq.Date(doy_frame[1],doy_frame[2], 1)),
                          outdir = paste0(WD, "/covariates_temp"), prefix = "susm_doy", progress=F, writechange=F, returnstack=F, overwrite=TRUE)
      bunk(6,"Temporal gap-filling...",pb,2)
      ## MODIS
      interpolateTemporal(s = LST_raster[[mask_dates_MODIS]], xin = as.numeric(as.Date(LST_dates[mask_dates_MODIS], "%Y_%m_%d")), xout = as.numeric(seq.Date(doy_frame[1],doy_frame[2], 1)),
                          outdir = paste0(WD, "/covariates_temp"), prefix = "LST_doy", progress=F, writechange=F, returnstack=F, overwrite=TRUE)
      bunk(6,"Temporal gap-filling...",pb,3)

      # For mac: deleting unwanted files
      xml_names <- list.files(paste0(WD, "/covariates_temp"), pattern = "xml", full.names = T)
      unlink(xml_names, recursive=TRUE)

      ssm_raster <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "ssm_doy", full.names = T))
      susm_raster <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "susm_doy", full.names = T))
      LST_raster <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "LST_doy", full.names = T))

      ## Sentinel-1
      if(!vv_na)
      {
        xin <- as.numeric(as.Date(S1_dates, "%Y%m%d"))
        xout <- as.numeric(seq.Date(doy_frame[1],doy_frame[2], 1))
        if (min(xin) > min(xout)) {
          vv_fill <- vv_raster[[1]] * 0 + mean(all$vv)
          vh_fill <- vh_raster[[1]] * 0 + mean(all$vh)
          angle_fill <- angle_raster[[1]] * 0 + mean(all$angle)
          for (i in seq(min(xout), min(xin)-1, 1)) {
            writeRaster(vv_fill, paste0(WD, "/covariates_temp/vv_doy_", i, ".tif"), overwrite=T)
            writeRaster(vh_fill, paste0(WD, "/covariates_temp/vh_doy_", i, ".tif"), overwrite=T)
            writeRaster(angle_fill, paste0(WD, "/covariates_temp/angle_doy_", i, ".tif"), overwrite=T)
          }
          xout <- seq(min(xin), max(xout), 1)
        }
        interpolateTemporal(s = vv_raster, xin = xin, xout = xout,
                            outdir = paste0(WD, "/covariates_temp"), prefix = "vv_doy", progress=F, writechange=F, returnstack=F, overwrite=TRUE)
        bunk(6,"Temporal gap-filling...",pb,4)
        interpolateTemporal(s = vh_raster, xin = xin, xout = xout,
                            outdir = paste0(WD, "/covariates_temp"), prefix = "vh_doy", progress=F, writechange=F, returnstack=F, overwrite=TRUE)
        bunk(6,"Temporal gap-filling...",pb,5)
        interpolateTemporal(s = angle_raster, xin = xin, xout = xout,
                            outdir = paste0(WD, "/covariates_temp"), prefix = "angle_doy", progress=F, writechange=F, returnstack=F, overwrite=TRUE)
        bunk(6,"Temporal gap-filling...",pb,6)

        xml_names <- list.files(paste0(WD, "/covariates_temp"), pattern = "xml", full.names = T)
        unlink(xml_names, recursive=TRUE)

        vv_raster <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "vv_doy", full.names = T))
        vh_raster <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "vh_doy", full.names = T))
        angle_raster <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "angle_doy", full.names = T))

        # ## If there are still NAs, fill NAs for each date for S1 images
        # for (i in 1:dim(vv_raster)[3])
        # {
        #   values(vv_raster[[i]])[is.na(values(vv_raster[[i]]))] = mean(values(vv_raster[[i]]), na.rm =T)
        #   values(vh_raster[[i]])[is.na(values(vh_raster[[i]]))] = mean(values(vh_raster[[i]]), na.rm =T)
        #   values(angle_raster[[i]])[is.na(values(angle_raster[[i]]))] = mean(values(angle_raster[[i]]), na.rm =T)
        #
        # }

      }

      ## if there is no Sentinel-1 data
      if(vv_na)
      {
        vv_raster <- LST_raster
        values(vv_raster) <- NA
        vh_raster <- angle_raster <- vv_raster

      }

      ## Landsat-8
      interpolateTemporal(s = LS_B4_raster[[mask_dates_LS]], xin = as.numeric(as.Date(LS_dates[mask_dates_LS], "%Y%m%d")), xout = as.numeric(seq.Date(doy_frame[1],doy_frame[2], 1)),
                          outdir = paste0(WD, "/covariates_temp"), prefix = "LS_B4_doy", progress=F, writechange=F, returnstack=F, overwrite=TRUE)
      bunk(6,"Temporal gap-filling...",pb,7)

      interpolateTemporal(s = LS_B5_raster[[mask_dates_LS]], xin = as.numeric(as.Date(LS_dates[mask_dates_LS], "%Y%m%d")), xout = as.numeric(seq.Date(doy_frame[1],doy_frame[2], 1)),
                          outdir = paste0(WD, "/covariates_temp"), prefix = "LS_B5_doy", progress=F, writechange=F, returnstack=F, overwrite=TRUE)
      bunk(6,"Temporal gap-filling...",pb,8)

      interpolateTemporal(s = LS_B6_raster[[mask_dates_LS]], xin = as.numeric(as.Date(LS_dates[mask_dates_LS], "%Y%m%d")), xout = as.numeric(seq.Date(doy_frame[1],doy_frame[2], 1)),
                          outdir = paste0(WD, "/covariates_temp"), prefix = "LS_B6_doy", progress=F, writechange=F, returnstack=F, overwrite=TRUE)
      bunk(6,"Temporal gap-filling...",pb,9)
      interpolateTemporal(s = LS_B7_raster[[mask_dates_LS]], xin = as.numeric(as.Date(LS_dates[mask_dates_LS], "%Y%m%d")), xout = as.numeric(seq.Date(doy_frame[1],doy_frame[2], 1)),
                          outdir = paste0(WD, "/covariates_temp"), prefix = "LS_B7_doy", progress=F, writechange=F, returnstack=F, overwrite=TRUE)
      bunk(6,"Temporal gap-filling...",pb,10)
      interpolateTemporal(s = LS_B10_raster[[mask_dates_LS]], xin = as.numeric(as.Date(LS_dates[mask_dates_LS], "%Y%m%d")), xout = as.numeric(seq.Date(doy_frame[1],doy_frame[2], 1)),
                          outdir = paste0(WD, "/covariates_temp"), prefix = "LS_B10_doy", progress=F, writechange=F, returnstack=F, overwrite=TRUE)
      bunk(6,"Temporal gap-filling...",pb,11)
      close(pb)

      xml_names <- list.files(paste0(WD, "/covariates_temp"), pattern = "xml", full.names = T)
      unlink(xml_names, recursive=TRUE)

      LS_B4_raster <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "LS_B4_doy", full.names = T))
      LS_B5_raster <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "LS_B5_doy", full.names = T))
      LS_B6_raster <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "LS_B6_doy", full.names = T))
      LS_B7_raster <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "LS_B7_doy", full.names = T))
      LS_B10_raster <- stack(list.files(paste0(WD, "/covariates_temp"), pattern = "LS_B10_doy", full.names = T))


      # ## If there are still NAs, fill NAs for each date for LS images
      # cat("\014")
      # writeLines(c("Tasks done: 7/10","Filling NAs..."))
      # pb <- txtProgressBar(min=0, max=dim(LS_B4_raster)[3], style=3)
      # for (i in 1:dim(LS_B4_raster)[3])
      # {
      #   values(LS_B4_raster[[i]])[is.na(values(LS_B4_raster[[i]]))] = mean(values(LS_B4_raster[[i]]), na.rm =T)
      #   values(LS_B5_raster[[i]])[is.na(values(LS_B5_raster[[i]]))] = mean(values(LS_B5_raster[[i]]), na.rm =T)
      #   values(LS_B6_raster[[i]])[is.na(values(LS_B6_raster[[i]]))] = mean(values(LS_B6_raster[[i]]), na.rm =T)
      #   values(LS_B7_raster[[i]])[is.na(values(LS_B7_raster[[i]]))] = mean(values(LS_B7_raster[[i]]), na.rm =T)
      #   values(LS_B10_raster[[i]])[is.na(values(LS_B10_raster[[i]]))] = mean(values(LS_B10_raster[[i]]), na.rm =T)
      #   bunk(7,"Filling NAs...",pb,i)
      # }
      # close(pb)

      cat("\014")
      writeLines(c("Tasks done: 7/10","Calculating NWDI..."))
      pb <- txtProgressBar(min=0, max=8, style=3)

      ## Calculate NDWI from Landsat-8
      LS_NDWI_raster <- (LS_B5_raster - LS_B7_raster)/(LS_B5_raster + LS_B7_raster)
      bunk(7, "Calculating NDVI...", pb, 1)
      ## Calculate NDVI from Landsat-8
      LS_NDVI_raster <- (LS_B5_raster - LS_B4_raster)/(LS_B5_raster + LS_B4_raster)
      bunk(7,"SSM aggregation...", pb, 2)


      ssm_criteria <- aggregate(ssm_raster[[1]], fact=round(10000/resolution), method='bilinear')

      if (dim(ssm_criteria[[1]])[1]>9){
        ## Resample SMAP
        ssm_raster2 <- aggregate(ssm_raster, fact=round(10000/resolution), method='bilinear')
        bunk(7,"SSM resampling...", pb, 3)
        ssm_raster3 <- resample(ssm_raster2, LS_B10_raster[[1]], method='bilinear')
        bunk(7,"SUSM aggregation...", pb, 4)
        ssm_raster <- ssm_raster3

        susm_raster2 <- aggregate(susm_raster, fact=round(10000/resolution), method='bilinear')
        bunk(7,"SUSM resampling...", pb, 5)
        susm_raster3 <- resample(susm_raster2, LS_B10_raster[[1]], method='bilinear')
        susm_raster <- susm_raster3
      }

      bunk(7,"LST aggregation...", pb, 6)

      ## Resample MODIS LST
      if (resolution<500) {
        LST_raster2 <- aggregate(LST_raster, fact=round(500/resolution), method='bilinear')
        bunk(7,"LST resampling...", pb, 7)
        LST_raster3 <- resample(LST_raster2, LS_B10_raster[[1]], method='bilinear')
        LST_raster <- LST_raster3
      }
      bunk(7,"LST resampling...", pb, 8)
      close(pb)

    }

  }

  ## Create folder to store downloaded covariates
  {
    dir.create(paste0(WD, "/covariates"))

    cat("\014")
    writeLines(c("Tasks done: 8/10","Saving covariates..."))
    pb <- txtProgressBar(min=0, max=17, style=3)
    writeRaster(constant_raster, file.path(paste0(WD, "/covariates"), "constant_raster.tif"), overwrite=TRUE)
    bunk(8,"Saving covariates...",pb,1)
    writeRaster(soil_5_raster, file.path(paste0(WD, "/covariates"), "soil_5_raster.tif"), overwrite=TRUE)
    bunk(8,"Saving covariates...",pb,2)
    writeRaster(soil_100_raster, file.path(paste0(WD, "/covariates"), "soil_100_raster.tif"), overwrite=TRUE)
    bunk(8,"Saving covariates...",pb,3)
    writeRaster(LC_raster, file.path(paste0(WD, "/covariates"), "LC_raster.tif"), overwrite=TRUE)
    bunk(8,"Saving covariates...",pb,4)

    writeRaster(vv_raster, file.path(paste0(WD, "/covariates"), "vv_raster.tif"), overwrite=TRUE)
    bunk(8,"Saving covariates...",pb,5)
    writeRaster(vh_raster, file.path(paste0(WD, "/covariates"), "vh_raster.tif"), overwrite=TRUE)
    bunk(8,"Saving covariates...",pb,6)
    writeRaster(angle_raster, file.path(paste0(WD, "/covariates"), "angle_raster.tif"), overwrite=TRUE)
    bunk(8,"Saving covariates...",pb,7)

    writeRaster(ssm_raster, file.path(paste0(WD, "/covariates"), "ssm_raster.tif"), overwrite=TRUE)
    bunk(8,"Saving covariates...",pb,8)
    writeRaster(susm_raster, file.path(paste0(WD, "/covariates"), "susm_raster.tif"), overwrite=TRUE)
    bunk(8,"Saving covariates...",pb,9)

    writeRaster(LST_raster, file.path(paste0(WD, "/covariates"), "LST_raster.tif"), overwrite=TRUE)
    bunk(8,"Saving covariates...",pb,10)

    writeRaster(LS_B4_raster, file.path(paste0(WD, "/covariates"), "LS_B4_raster.tif"), overwrite=TRUE)
    bunk(8,"Saving covariates...",pb,11)
    writeRaster(LS_B5_raster, file.path(paste0(WD, "/covariates"), "LS_B5_raster.tif"), overwrite=TRUE)
    bunk(8,"Saving covariates...",pb,12)
    writeRaster(LS_B6_raster, file.path(paste0(WD, "/covariates"), "LS_B6_raster.tif"), overwrite=TRUE)
    bunk(8,"Saving covariates...",pb,13)
    writeRaster(LS_B7_raster, file.path(paste0(WD, "/covariates"), "LS_B7_raster.tif"), overwrite=TRUE)
    bunk(8,"Saving covariates...",pb,14)
    writeRaster(LS_B10_raster, file.path(paste0(WD, "/covariates"), "LS_B10_raster.tif"), overwrite=TRUE)
    bunk(8,"Saving covariates...",pb,15)
    writeRaster(LS_NDVI_raster, file.path(paste0(WD, "/covariates"), "LS_NDVI_raster.tif"), overwrite=TRUE)
    bunk(8,"Saving covariates...",pb,16)
    writeRaster(LS_NDWI_raster, file.path(paste0(WD, "/covariates"), "LS_NDWI_raster.tif"), overwrite=TRUE)
    bunk(8,"Saving covariates...",pb,17)
    close(pb)
  }

  ### Read all covariate rasters
  {
    constant_raster <- stack(paste0(WD, "/covariates/constant_raster.tif"))
    soil_5_raster <- stack(paste0(WD, "/covariates/soil_5_raster.tif"))
    soil_100_raster <- stack(paste0(WD, "/covariates/soil_100_raster.tif"))
    LC_raster <- raster(paste0(WD, "/covariates/LC_raster.tif"))
    vv_raster <- stack(paste0(WD, "/covariates/vv_raster.tif"))
    vh_raster <- stack(paste0(WD, "/covariates/vh_raster.tif"))
    angle_raster <- stack(paste0(WD, "/covariates/angle_raster.tif"))
    ssm_raster <- stack(paste0(WD, "/covariates/ssm_raster.tif"))
    susm_raster <- stack(paste0(WD, "/covariates/susm_raster.tif"))
    LST_raster <- stack(paste0(WD, "/covariates/LST_raster.tif"))
    LS_B4_raster <- stack(paste0(WD, "/covariates/LS_B4_raster.tif"))
    LS_B5_raster <- stack(paste0(WD, "/covariates/LS_B5_raster.tif"))
    LS_B6_raster <- stack(paste0(WD, "/covariates/LS_B6_raster.tif"))
    LS_B7_raster <- stack(paste0(WD, "/covariates/LS_B7_raster.tif"))
    LS_B10_raster <- stack(paste0(WD, "/covariates/LS_B10_raster.tif"))
    LS_NDVI_raster <- stack(paste0(WD, "/covariates/LS_NDVI_raster.tif"))
    LS_NDWI_raster <- stack(paste0(WD, "/covariates/LS_NDWI_raster.tif"))


    #### Delete temporary datasets
    unlink(paste0(WD, "/covariates_temp"), recursive=TRUE)


  }

  ## Create output folders
  dir.create(paste0(WD, "/VWC"))

  ## Generate a regular grid from the ROI for mapping
  grid <- rasterToPoints(constant_raster, spatial = TRUE)
  grid_US <- as.data.frame(spTransform(grid, CRS("+init=epsg:5070")))
  grid_US.shape <- spTransform(grid, CRS("+init=epsg:5070"))
  grid_US_raster <- rasterFromXYZ(grid_US[,c("x","y","x")], res= resolution, crs =  CRS("+init=epsg:5070"), digits = 0)

  ## Get the list of targeted dates
  dates <- seq.Date(doy_frame[1],doy_frame[2], 1)

  ################## ML models for surface and rootzone soil moisture maps
  {
    ## For each date, apply the established ML model (reduced model) to the covariate rasters
    cat("\014")
    writeLines(c("Tasks done: 9/10","Applying ML model..."))
    pb <- txtProgressBar(min=0, max=length(dates), style=3)
    for (i in 1:length(dates))
    {
      print(paste0("now ", i, " day in ", length(dates), " days"))
      ## Combine all the covariates
      covariates_others <- list(constant_raster,
                                LC_raster,
                                vv_raster[[i]], vh_raster[[i]], angle_raster[[i]],
                                ssm_raster[[i]], susm_raster[[i]],
                                LST_raster[[i]],
                                LS_B5_raster[[i]], LS_B6_raster[[i]], LS_B7_raster[[i]], LS_B10_raster[[i]], LS_NDVI_raster[[i]], LS_NDWI_raster[[i]])

      ## Extract covariates to the mapping grid
      map_covariate <- as.data.frame(do.call("cbind", lapply(covariates_others, raster::extract, grid)))
      # cat("\014")
      # setTxtProgressBar(pb, i-0.5)

      colnames(map_covariate) <- c("elevation", "slope", "aspect", "hillshade",
                                   "clay_5", "sand_5", "bd_5",
                                   "clay_100", "sand_100", "bd_100",
                                   "LC",
                                   "vv", "vh", "angle",
                                   "ssm", "susm",
                                   "LST",
                                   "LS_B5", "LS_B6", "LS_B7", "LS_B10", "LS_NDVI", "LS_NDWI")
      map_covariate <- cbind(coordinates(grid), map_covariate)

      ## Convert land cover to categorical data
      map_covariate$landcover <- "Others"

      if (sum(map_covariate$LC==21|map_covariate$LC==22|map_covariate$LC==23|map_covariate$LC==24)>0)
        map_covariate[map_covariate$LC==21|map_covariate$LC==22|map_covariate$LC==23|map_covariate$LC==24,]$landcover <- "Developed"
      if (sum(map_covariate$LC==31)>0)
        map_covariate[map_covariate$LC==31,]$landcover <- "Barren"

      if (sum(map_covariate$LC==41|map_covariate$LC==42|map_covariate$LC==43)>0)
        map_covariate[map_covariate$LC==41|map_covariate$LC==42|map_covariate$LC==43,]$landcover <- "Forest"

      if (sum(map_covariate$LC==51|map_covariate$LC==52)>0)
        map_covariate[map_covariate$LC==51|map_covariate$LC==52,]$landcover <- "Shrub"

      if (sum(map_covariate$LC==71|map_covariate$LC==72)>0)
        map_covariate[map_covariate$LC==71|map_covariate$LC==72,]$landcover <- "Grassland"

      if (sum(map_covariate$LC==81)>0)
        map_covariate[map_covariate$LC==81,]$landcover <- "Pasture"

      if (sum(map_covariate$LC==82)>0)
        map_covariate[map_covariate$LC==82,]$landcover <- "Crop"

      if (sum(map_covariate$LC==90|map_covariate$LC==95)>0)
        map_covariate[map_covariate$LC==90|map_covariate$LC==95,]$landcover <- "Wetland"


      map_covariate <- map_covariate[!map_covariate$landcover=="Others",]


      index_grid <- names(map_covariate) %in% c( "x","y", "landcover",
                                                 "elevation", "slope", "aspect", "hillshade",
                                                 "clay_5", "sand_5", "bd_5",
                                                 "clay_100", "sand_100", "bd_100",
                                                 "ssm", "susm",
                                                 "vv", "vh", "angle",
                                                 "LST",
                                                 "LS_B5", "LS_B6", "LS_B7", "LS_B10", "LS_NDVI", "LS_NDWI")

      brige <- all[,c( "Longitude","Latitude","landcover",
                        "elevation", "slope", "aspect", "hillshade",
                        "clay_5", "sand_5", "bd_5",
                        "clay_100", "sand_100", "bd_100",
                        "ssm","susm",
                        "vv", "vh", "angle",
                        "LST",
                        "LS_B5", "LS_B6", "LS_B7", "LS_B10", "LS_NDVI", "LS_NDWI")]

      colnames(brige) <- c( "x","y", "landcover",
                            "elevation", "slope", "aspect", "hillshade",
                            "clay_5", "sand_5", "bd_5",
                            "clay_100", "sand_100", "bd_100",
                            "ssm","susm",
                            "vv", "vh", "angle",
                            "LST",
                            "LS_B5", "LS_B6", "LS_B7", "LS_B10", "LS_NDVI", "LS_NDWI")

      map_covariate_temp <- rbind(brige, map_covariate[,index_grid])

      map_covariate_temp$landcover <- as.factor(map_covariate_temp$landcover)
      map_covariate_grid <- map_covariate_temp[-c(1:nrow(brige)),]


      for (w in 1:ncol(map_covariate_grid))
      {
        if (!is.factor(map_covariate_grid[,w] ))
          map_covariate_grid[is.na(map_covariate_grid[,w]),w]  <- mean(map_covariate_grid[,w], na.rm =T)

      }


      if(sum(!is.na(map_covariate_grid$vv))<5)
      {
        map_covariate_grid$vv  <- mean(all$vv)
        map_covariate_grid$vh  <- mean(all$vh)
        map_covariate_grid$angle  <- mean(all$angle)
      }



      ## Apply ML models to the grid
      ## Reduced model
      grid_VWC_5_mean   <-  predict(model_surface, newdata = map_covariate_grid, what = mean)
      grid_VWC_5_sd     <-  predict(model_surface, newdata = map_covariate_grid, what = sd)
      grid_VWC_100_mean <-  predict(model_rz, newdata = map_covariate_grid, what = mean)
      grid_VWC_100_sd   <-  predict(model_rz, newdata = map_covariate_grid, what = sd)

      ## Output predictions as rasters
      VWC_5_mean_raster <- rasterFromXYZ(cbind(map_covariate_grid[,c("x","y")], grid_VWC_5_mean), res=res(constant_raster[[1]]), crs =   "+proj=longlat +datum=WGS84 +no_defs ")
      VWC_5_mean_raster <- projectRaster(mask(VWC_5_mean_raster, boundary), grid_US_raster, method='bilinear')
      writeRaster(round(VWC_5_mean_raster, digits=3), file.path(paste0(WD, "/VWC"), paste("VWC_5_mean_", dates[i], ".tif", sep = "")), overwrite=TRUE)
      VWC_5_sd_raster <- rasterFromXYZ(cbind(map_covariate_grid[,c("x","y")], grid_VWC_5_sd), res=res(constant_raster[[1]]), crs =   "+proj=longlat +datum=WGS84 +no_defs ")
      VWC_5_sd_raster <- projectRaster(mask(VWC_5_sd_raster, boundary), grid_US_raster, method='bilinear')
      writeRaster(round(VWC_5_sd_raster, digits=3), file.path(paste0(WD, "/VWC"), paste("VWC_5_sd_", dates[i], ".tif", sep = "")), overwrite=TRUE)

      VWC_100_mean_raster <- rasterFromXYZ(cbind(map_covariate_grid[,c("x","y")], grid_VWC_100_mean), res=res(constant_raster[[1]]), crs =   "+proj=longlat +datum=WGS84 +no_defs ")
      VWC_100_mean_raster <- projectRaster(mask(VWC_100_mean_raster, boundary), grid_US_raster, method='bilinear')
      writeRaster(round(VWC_100_mean_raster, digits=3), file.path(paste0(WD, "/VWC"), paste("VWC_100_mean_", dates[i], ".tif", sep = "")), overwrite=TRUE)
      VWC_100_sd_raster <- rasterFromXYZ(cbind(map_covariate_grid[,c("x","y")], grid_VWC_100_sd), res=res(constant_raster[[1]]), crs =   "+proj=longlat +datum=WGS84 +no_defs ")
      VWC_100_sd_raster <- projectRaster(mask(VWC_100_sd_raster, boundary), grid_US_raster, method='bilinear')
      writeRaster(round(VWC_100_sd_raster, digits=3), file.path(paste0(WD, "/VWC"), paste("VWC_100_sd_", dates[i], ".tif", sep = "")), overwrite=TRUE)

      # Percentile
      if(percentile){
        grid_VWC_5_lower   <-  predict(model_surface, newdata = map_covariate_grid, what = 0.05)
        VWC_5_lower_raster <- rasterFromXYZ(cbind(map_covariate_grid[,c("x","y")], grid_VWC_5_lower), res=res(constant_raster[[1]]), crs =   "+proj=longlat +datum=WGS84 +no_defs ")
        VWC_5_lower_raster <- projectRaster(mask(VWC_5_lower_raster, boundary), grid_US_raster, method='bilinear')
        writeRaster(round(VWC_5_lower_raster, digits=3), file.path(paste0(WD, "/VWC"), paste("VWC_5_lower_", dates[i], ".tif", sep = "")), overwrite=TRUE)

        grid_VWC_5_upper   <-  predict(model_surface, newdata = map_covariate_grid, what = 0.95)
        VWC_5_upper_raster <- rasterFromXYZ(cbind(map_covariate_grid[,c("x","y")], grid_VWC_5_upper), res=res(constant_raster[[1]]), crs =   "+proj=longlat +datum=WGS84 +no_defs ")
        VWC_5_upper_raster <- projectRaster(mask(VWC_5_upper_raster, boundary), grid_US_raster, method='bilinear')
        writeRaster(round(VWC_5_upper_raster, digits=3), file.path(paste0(WD, "/VWC"), paste("VWC_5_upper_", dates[i], ".tif", sep = "")), overwrite=TRUE)

        grid_VWC_100_lower   <-  predict(model_rz, newdata = map_covariate_grid, what = 0.05)
        VWC_100_lower_raster <- rasterFromXYZ(cbind(map_covariate_grid[,c("x","y")], grid_VWC_100_lower), res=res(constant_raster[[1]]), crs =   "+proj=longlat +datum=WGS84 +no_defs ")
        VWC_100_lower_raster <- projectRaster(mask(VWC_100_lower_raster, boundary), grid_US_raster, method='bilinear')
        writeRaster(round(VWC_100_lower_raster, digits=3), file.path(paste0(WD, "/VWC"), paste("VWC_100_lower_", dates[i], ".tif", sep = "")), overwrite=TRUE)

        grid_VWC_100_upper   <-  predict(model_rz, newdata = map_covariate_grid, what = 0.95)
        VWC_100_upper_raster <- rasterFromXYZ(cbind(map_covariate_grid[,c("x","y")], grid_VWC_100_upper), res=res(constant_raster[[1]]), crs =   "+proj=longlat +datum=WGS84 +no_defs ")
        VWC_100_upper_raster <- projectRaster(mask(VWC_100_upper_raster, boundary), grid_US_raster, method='bilinear')
        writeRaster(round(VWC_100_upper_raster, digits=3), file.path(paste0(WD, "/VWC"), paste("VWC_100_upper_", dates[i], ".tif", sep = "")), overwrite=TRUE)
      }

      if(i<length(dates)) bunk(9,"Applying ML model...",pb,i)
      if(i==length(dates)) bunk(10,"Done.",pb,i)
    }
    close(pb)
  }

  xml_names <- list.files(paste0(WD, "/VWC"), pattern = "xml", full.names = T)
  unlink(xml_names, recursive=TRUE)

  xml_names <- list.files(paste0(WD, "/covariates"), pattern = "xml", full.names = T)
  unlink(xml_names, recursive=TRUE)

  # Clean memory
  gc()
  unlink(file.path(tempdir(), "raster"), recursive=T)
  ee_clean_container(name = "rgee_backup", type = "drive", quiet = FALSE)
}
