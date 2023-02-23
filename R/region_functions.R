#' Daily map visualization
#'
#' Visualizes the maps created by the regional prediction function (one day, one depth at a time).
#'
#' @param date  character. Defines the date of the visualized maps
#' @param depth numeric. Defines the depth of the VWC values to visualize.
#' @param percentile logical. If TRUE, the upper and lower bounds of the 90\% Confidence Interval of the prediction will be included.
#' @param project character. Defines the name of the folder where the maps are saved. If not defined, the function will look for maps in the working directory.
#'
#' @return a leaflet map containing the VWC maps.
#'
#' @import leaflet raster RColorBrewer
#'
#' @export
#'
#' @examples
#' \dontrun{ plot_map("2020-06-15", 5, TRUE, project="WI_region") }
plot_map <- function(date, depth, percentile = FALSE, project=NULL){
  setwd(getwd())
  WD <- getwd()
  if(!is.null(project)){
    if (!file.exists(project)){
      return("Project folder doesn't exist.")
    } else {
      WD <- paste0(WD, "/", project)
    }
  }
  WD <- paste0(WD,"/VWC")

  xml_names <- list.files(paste0(WD, "/VWC"), pattern = "xml", full.names = T)
  unlink(xml_names, recursive=TRUE)

  mode <- c("mean", "sd")
  if (percentile==TRUE) mode <- c(mode, "lower", "upper")

  maps <- stack(paste0(WD, "/", paste0("VWC_", depth, "_", mode, "_", date, ".tif")))
  mode <- c("Mean", "SD")
  if (percentile==TRUE) mode <- c(mode, "0.05", "0.95")

  m <- leaflet() %>%
    addProviderTiles("Esri.WorldImagery")

  for (i in seq(length(mode))){
    vwc_val <- maps[[i]]
    min_val <- min(values(vwc_val),na.rm=T)
    max_val <- max(values(vwc_val),na.rm=T)
    pal1 <- colorNumeric(
      brewer.pal(n = 5, name = "RdYlBu"),
      domain = c(min_val, max_val),
      na.color = "transparent")
    m <- m %>%
      addRasterImage(vwc_val, colors=pal1, opacity = 0.7, group=mode[i], project=FALSE) %>%
      addLegend("bottomright", opacity=0.8, pal = pal1, values = values(vwc_val),
                bins = seq(min_val, max_val, (max_val-min_val)/6),
                title=mode[i], group=mode[i])
  }

  m <- m %>%
    addLayersControl(
      baseGroups = mode,
      options = layersControlOptions(collapsed = FALSE)
    )%>%
    htmlwidgets::onRender("
        function(el, x) {
          var updateLegend = function () {
              var selectedGroup = document.querySelectorAll('input:checked')[0].nextSibling.innerText.substr(1);

              document.querySelectorAll('.legend').forEach(a => a.hidden=true);
              document.querySelectorAll('.legend').forEach(l => {
                if (l.children[0].children[0].innerText == selectedGroup) l.hidden=false;
              });
          };
          updateLegend();
          this.on('baselayerchange', e => updateLegend());
        }") # https://stackoverflow.com/questions/52393310/is-it-possible-to-switch-between-multiple-legends-when-switching-between-base-gr

  return(m)
}



#' Temporal variation for the entire area
#'
#' Visualizes the change in basic statistics (mean, min, median, max) of the regional VWC values across time.
#'
#' @param depth numeric. Defines the depth of the VWC values to visualize.
#' @param project character. Defines the name of the folder where the maps are saved. If not defined, the function will look for maps in the working directory.
#'
#' @return a time series plot visualizing the basic statistics of the region.
#'
#' @import tidyverse raster viridis
#'
#' @export
#'
#' @examples
#' \dontrun{ area_sum(5, project="WI_region") }
area_sum <- function(depth, project=NULL){
  setwd(getwd())
  WD <- getwd()
  if(!is.null(project)){
    if (!file.exists(project)){
      return("Project folder doesn't exist.")
    } else {
      WD <- paste0(WD, "/", project)
    }
  }

  if(!file.exists(paste0(WD, paste0("/VWC_", depth, "_ts_data.csv")))) {
    xml_names <- list.files(paste0(WD, "/VWC"), pattern = "xml", full.names = T)
    unlink(xml_names, recursive=TRUE)

    files <- list.files(paste0(WD,"/VWC"), paste0("VWC_",as.character(depth),"_mean"))
    rasters <- stack(paste0(WD,"/VWC/",files))

    dates <- gsub(".*_(\\d*-\\d*-\\d*).tif", "\\1", files)

    ts_table <- data.frame(Date = rep(dates,each=4),
                           Summary = rep(c("Mean", "Median", "Min", "Max"), times=length(dates)))

    print("Calculating time series variables...")
    pb <- txtProgressBar(min = 0, max = length(dates), style = 3)
    for (i in seq(length(dates))){
      ts_table$value[i*4-3] <- round(mean(values(rasters[[i]]), na.rm=T), digits=3)
      ts_table$value[i*4-2] <- round(median(values(rasters[[i]]), na.rm=T), digits=3)
      ts_table$value[i*4-1] <- round(min(values(rasters[[i]]), na.rm=T), digits=3)
      ts_table$value[i*4] <- round(max(values(rasters[[i]]), na.rm=T), digits=3)
      setTxtProgressBar(pb, i)
    }
    close(pb)

    write.csv(ts_table, paste0(WD, paste0("/VWC_", depth, "_ts_data.csv")), row.names=F)
  } else {
    ts_table <- read.csv(paste0(WD, paste0("/VWC_", depth, "_ts_data.csv")))
  }

  p <- ggplot(ts_table)+
    geom_line(aes(as.Date(Date),value,col=Summary))+
    labs(y="VWC", x="Date", color="Summary", title="Area Aggregated Statistics")+
    theme_bw()+
    scale_x_date(date_labels =  "%y-%m-%d")+
    theme(axis.text.x = element_text(angle = -10, vjust = 0.1, hjust=0.5,size=8))

  return(p)


}



#' Pixel-wise variation over time
#'
#' Visualizes the pixel-wise variations in the regional VWC values across time.
#' NOTE: very time consuming for maps with high resolution and large area.
#'
#' @param depth numeric. Defines the depth of the VWC values to visualize.
#' @param project character. Defines the name of the folder where the maps are saved. If not defined, the function will look for maps in the working directory.
#'
#' @return a leaflet plot visualizing the basic statistics (mean, min, median, max) of the pixel-wise variation in the region.
#'
#' @import leaflet raster
#'
#' @export
#'
#' @examples
#' \dontrun{ pixelwise_sum(5, project="WI_region") }
pixelwise_sum <- function(depth, project=NULL){
  setwd(getwd())
  WD <- getwd()
  if(!is.null(project)){
    if (!file.exists(project)){
      return("Project folder doesn't exist.")
    } else {
      WD <- paste0(WD, "/", project)
      if(!dir.exists(paste0(WD,"/temporal_VWC"))){
        dir.create(paste0(WD, "/temporal_VWC"))
      }
    }
  }

  {
    if (!dir.exists(paste0(WD,"/temporal_VWC/VWC_", depth, "cm"))) {
      dir.create(paste0(WD, "/temporal_VWC/VWC_", depth, "cm"))

      xml_names <- list.files(paste0(WD, "/VWC"), pattern = "xml", full.names = T)
      unlink(xml_names, recursive=TRUE)

      files <- list.files(paste0(WD,"/VWC"), paste0("VWC_",as.character(depth),"_mean"))
      rasters <- stack(paste0(WD,"/VWC/",files))

      print("Calculating temporal VWC... (NOTE: this may take a long time.)")
      pb <- txtProgressBar(min = 0, max = 6, style = 3)
      raster_min <- calc(rasters, min)
      setTxtProgressBar(pb, 1)
      raster_median <- calc(rasters, median)
      setTxtProgressBar(pb, 2)
      raster_max <- calc(rasters, max)
      setTxtProgressBar(pb, 3)
      raster_mean <- calc(rasters, mean)
      setTxtProgressBar(pb, 4)
      raster_sd <- calc(rasters, sd)
      setTxtProgressBar(pb, 5)
      raster_range <- calc(rasters, function(i){max(i) - min(i)})
      setTxtProgressBar(pb, 6)
      close(pb)

      writeRaster(round(raster_min, digits=3), file.path(paste0(WD, "/temporal_VWC/VWC_", depth, "cm"), paste0("VWC_", depth, "_min.tif")), overwrite=TRUE)
      writeRaster(round(raster_median, digits=3), file.path(paste0(WD, "/temporal_VWC/VWC_", depth, "cm"), paste0("VWC_", depth, "_median.tif")), overwrite=TRUE)
      writeRaster(round(raster_max, digits=3), file.path(paste0(WD, "/temporal_VWC/VWC_", depth, "cm"), paste0("VWC_", depth, "_max.tif")), overwrite=TRUE)
      writeRaster(round(raster_mean, digits=3), file.path(paste0(WD, "/temporal_VWC/VWC_", depth, "cm"), paste0("VWC_", depth, "_mean.tif")), overwrite=TRUE)
      writeRaster(round(raster_sd, digits=3), file.path(paste0(WD, "/temporal_VWC/VWC_", depth, "cm"), paste0("VWC_", depth, "_sd.tif")), overwrite=TRUE)
      writeRaster(round(raster_range, digits=3), file.path(paste0(WD, "/temporal_VWC/VWC_", depth, "cm"), paste0("VWC_", depth, "_range.tif")), overwrite=TRUE)

      xml_names <- list.files(paste0(WD, "/temporal_VWC/VWC_", depth, "cm"), pattern = "xml", full.names = T)
      unlink(xml_names, recursive=TRUE)
    } else {
      raster_min <- raster(paste0(WD, "/temporal_VWC/VWC_", depth, "cm/VWC_",depth, "_min.tif"))
      raster_median <- raster(paste0(WD, "/temporal_VWC/VWC_", depth, "cm/VWC_",depth, "_median.tif"))
      raster_max <- raster(paste0(WD, "/temporal_VWC/VWC_", depth, "cm/VWC_",depth, "_max.tif"))
      raster_mean <- raster(paste0(WD, "/temporal_VWC/VWC_", depth, "cm/VWC_",depth, "_mean.tif"))
      raster_sd <- raster(paste0(WD, "/temporal_VWC/VWC_", depth, "cm/VWC_",depth, "_sd.tif"))
      raster_range <- raster(paste0(WD, "/temporal_VWC/VWC_", depth, "cm/VWC_",depth, "_range.tif"))
    }

    r <- c(raster_mean, raster_min, raster_median, raster_max, raster_sd, raster_range)
    m <- c("Mean", "Min", "Median", "Max", "SD", "Range")

    map <- leaflet() %>%
      addProviderTiles("Esri.WorldImagery")

    print("Plotting temporal VWC...")
    pb <- txtProgressBar(min=0, max=length(r), style=3)
    for (i in seq(length(r))){
      vwc_val <- r[[i]]
      min_val <- min(values(vwc_val),na.rm=T)
      max_val <- max(values(vwc_val),na.rm=T)
      pal1 <- colorNumeric(
        brewer.pal(n = 5, name = "RdYlBu"),
        domain = c(min_val, max_val),
        na.color = "transparent")
      map <- map %>%
        addRasterImage(vwc_val, colors=pal1, opacity = 0.7, group=m[i], project=FALSE) %>%
        addLegend("bottomright", opacity=0.8, pal = pal1, values = values(vwc_val),
                  bins = seq(min_val, max_val, (max_val-min_val)/6),
                  title=m[i], group=m[i])
      setTxtProgressBar(pb,i)
    }
    close(pb)

    map <- map %>%
      addLayersControl(
        baseGroups = m,
        options = layersControlOptions(collapsed = FALSE)
      )%>%
      htmlwidgets::onRender("
      function(el, x) {
        var updateLegend = function () {
            var selectedGroup = document.querySelectorAll('input:checked')[0].nextSibling.innerText.substr(1);

            document.querySelectorAll('.legend').forEach(a => a.hidden=true);
            document.querySelectorAll('.legend').forEach(l => {
              if (l.children[0].children[0].innerText == selectedGroup) l.hidden=false;
            });
        };
        updateLegend();
        this.on('baselayerchange', e => updateLegend());
      }") # https://stackoverflow.com/questions/52393310/is-it-possible-to-switch-between-multiple-legends-when-switching-between-base-gr

    return(map)
  }
}



#' Aggregation by interval
#'
#' Aggregates VWC values according to starting date & time interval and computes the basic statistics (mean, min, max, sd) of each time interval.
#'
#' @param depth numeric. Defines the depth of the VWC values to aggregate.
#' @param start_date character. Defines the start date of the aggregation period.
#' @param end_date character. Defines the end date of the aggregation period.
#' @param frequency numeric. Defines the number of days of the aggregation time interval.
#' @param project character. Defines the name of the folder where the maps are saved. If not defined, the function will look for maps in the working directory.
#'
#' @return None
#'
#' @import raster tidyverse
#'
#' @export
#'
#' @examples
#' \dontrun{ aggregate_interval(5, "2020-06-15", "2020-08-15", 7, "WI_region") }
aggregate_interval <- function(depth, start_date, end_date, frequency, project=NULL){
  setwd(getwd())
  WD <- getwd()
  if(!is.null(project)){
    if (!file.exists(project)){
      return("Project folder doesn't exist.")
    } else {
      WD <- paste0(WD, "/", project)
    }
  }

  xml_names <- list.files(paste0(WD, "/VWC"), pattern = "xml", full.names = T)
  unlink(xml_names, recursive=TRUE)

  files <- list.files(paste0(WD,"/VWC"), paste0("VWC_",as.character(depth),"_mean"))
  start_idx <- match(paste0("VWC_",as.character(depth),"_mean_",start_date,".tif"), files)
  end_idx <- match(paste0("VWC_",as.character(depth),"_mean_",end_date,".tif"), files)
  files <- files[start_idx :end_idx]
  groups <- split(files, ceiling(seq_along(files)/frequency))
  dates <- gsub(".*_(\\d*-\\d*-\\d*).tif", "\\1",sapply(groups,"[[",1))

  if(!dir.exists(paste0(WD, "/VWC_aggregation"))) dir.create(paste0(WD, "/VWC_aggregation"))
  if(!dir.exists(paste0(WD, "/VWC_aggregation/", start_date, "_", end_date))) dir.create(paste0(WD, "/VWC_aggregation/", start_date, "_", end_date))

  if(!dir.exists(paste0(WD, paste0("/VWC_aggregation/", start_date, "_", end_date, "/", frequency,"days_", depth, "cm")))){
    dir.create(paste0(WD, paste0("/VWC_aggregation/", start_date, "_", end_date, "/", frequency,"days_", depth, "cm")))
    save_dir <- paste0(WD, paste0("/VWC_aggregation/", start_date, "_", end_date, "/", frequency,"days_", depth, "cm"))
    dir.create(paste0(save_dir, "/aggregated_VWC"))

    sub_mean <- NULL
    sub_min <- NULL
    sub_max <- NULL
    sub_sd <- NULL
    sub_median <- NULL

    print("Calculating aggregated values...")
    pb <- txtProgressBar(min=0, max=length(groups), style=3)
    for (i in seq(length(groups))){
      if (length(groups[[i]])<2) {
        setTxtProgressBar(pb, i)
        print(paste0("Eliminate ", groups[[i]], " - interval only has one day."))
        break
      }

      rasters <- stack(paste0(WD,"/VWC/",groups[[i]]))

      raster_min <- calc(rasters, min)
      raster_max <- calc(rasters, max)
      raster_mean <- calc(rasters, mean)
      raster_sd <- calc(rasters, sd)
      raster_median <- calc(rasters, median)

      sub_min <-rbind(sub_min, as.data.frame(round(raster_min, digits=3), xy=TRUE, na.rm=T) %>%
                        mutate(date=dates[[i]]))
      sub_max <-rbind(sub_max, as.data.frame(round(raster_max, digits=3), xy=TRUE, na.rm=T) %>%
                        mutate(date=dates[[i]]))
      sub_mean <-rbind(sub_mean, as.data.frame(round(raster_mean, digits=3), xy=TRUE, na.rm=T) %>%
                         mutate(date=dates[[i]]))
      sub_sd <-rbind(sub_sd, as.data.frame(round(raster_sd, digits=3), xy=TRUE, na.rm=T) %>%
                       mutate(date=dates[[i]]))
      sub_median <-rbind(sub_median, as.data.frame(round(raster_median, digits=3), xy=TRUE, na.rm=T) %>%
                       mutate(date=dates[[i]]))

      writeRaster(round(raster_min, digits=3), file.path(paste0(save_dir, "/aggregated_VWC/min_", dates[[i]], ".tif")), overwrite=TRUE)
      writeRaster(round(raster_max, digits=3), file.path(paste0(save_dir, "/aggregated_VWC/max_", dates[[i]], ".tif")), overwrite=TRUE)
      writeRaster(round(raster_mean, digits=3), file.path(paste0(save_dir, "/aggregated_VWC/mean_", dates[[i]], ".tif")), overwrite=TRUE)
      writeRaster(round(raster_sd, digits=3), file.path(paste0(save_dir, "/aggregated_VWC/sd_", dates[[i]], ".tif")), overwrite=TRUE)
      writeRaster(round(raster_median, digits=3), file.path(paste0(save_dir, "/aggregated_VWC/median_", dates[[i]], ".tif")), overwrite=TRUE)

      setTxtProgressBar(pb, i)
    }
    close(pb)

    xml_names <- list.files(paste0(save_dir, "/aggregated_VWC"), pattern = "xml", full.names = T)
    unlink(xml_names, recursive=TRUE)

    write.csv(sub_min, paste0(save_dir, "/", frequency, "_day_min.csv"), row.names=F)
    write.csv(sub_max, paste0(save_dir, "/", frequency, "_day_max.csv"), row.names=F)
    write.csv(sub_mean, paste0(save_dir, "/", frequency, "_day_mean.csv"), row.names=F)
    write.csv(sub_sd, paste0(save_dir, "/", frequency, "_day_sd.csv"), row.names=F)
    write.csv(sub_median, paste0(save_dir, "/", frequency, "_day_median.csv"), row.names=F)

    return(dates)
  } else {
    print("Files already exist.")
    return(dates)
  }
}



#' Aggregation by dates
#'
#' Summarizes the basic statistics (mean, min, max, sd) of VWC values within one time interval.
#'
#' @param depth numeric. Defines the depth of the VWC values to visualize.
#' @param start_date character. Defines the start date of the aggregation period.
#' @param end_date character. Defines the end date of the aggregation period.
#' @param project character. Defines the name of the folder where the maps are saved. If not defined, the function will look for maps in the working directory.
#'
#' @return None
#'
#' @import raster
#'
#' @export
#'
#' @examples
#' \dontrun{ aggregate_sum(5, "2020-06-15", "2020-08-15", project="WI_region") }
aggregate_sum <- function(depth, start_date, end_date, project=NULL) {
  setwd(getwd())
  WD <- getwd()
  if(!is.null(project)){
    if (!file.exists(project)){
      return("Project folder doesn't exist.")
    } else {
      WD <- paste0(WD, "/", project)
    }
  }

  if (!dir.exists(paste0(WD, "/VWC_aggregation"))) dir.create(paste0(WD, "/VWC_aggregation"))
  if(!dir.exists(paste0(WD, "/VWC_aggregation/", start_date, "_", end_date))) dir.create(paste0(WD, "/VWC_aggregation/", start_date, "_", end_date))

  if (!dir.exists(paste0(WD, paste0("/VWC_aggregation/", start_date,"_", end_date, "/",depth,"cm")))){
    dir.create(paste0(WD, paste0("/VWC_aggregation/", start_date,"_", end_date, "/",depth,"cm")))
    save_dir <- paste0(WD, paste0("/VWC_aggregation/", start_date,"_", end_date, "/",depth,"cm"))
    dir.create(paste0(save_dir, "/aggregated_VWC"))

    xml_names <- list.files(paste0(WD, "/VWC"), pattern = "xml", full.names = T)
    unlink(xml_names, recursive=TRUE)

    files <- list.files(paste0(WD,"/VWC"), paste0("VWC_",as.character(depth),"_mean"))
    start_idx <- match(paste0("VWC_",as.character(depth),"_mean_",start_date,".tif"), files)
    end_idx <- match(paste0("VWC_",as.character(depth),"_mean_",end_date,".tif"), files)
    files <- files[start_idx :end_idx]

    rasters <- stack(paste0(WD,"/VWC/",files))

    print("Calculating aggregated values...")
    pb <- txtProgressBar(min=0, max=5, style=3)
    raster_min <- calc(rasters, min)
    setTxtProgressBar(pb, 1)
    raster_max <- calc(rasters, max)
    setTxtProgressBar(pb, 2)
    raster_mean <- calc(rasters, mean)
    setTxtProgressBar(pb, 3)
    raster_sd <- calc(rasters, sd)
    setTxtProgressBar(pb, 4)
    raster_median <- calc(rasters, median)
    setTxtProgressBar(pb, 5)

    writeRaster(round(raster_min, digits=3), file.path(paste0(save_dir, "/aggregated_VWC/min.tif")), overwrite=TRUE)
    writeRaster(round(raster_max, digits=3), file.path(paste0(save_dir, "/aggregated_VWC/max.tif")), overwrite=TRUE)
    writeRaster(round(raster_mean, digits=3), file.path(paste0(save_dir, "/aggregated_VWC/mean.tif")), overwrite=TRUE)
    writeRaster(round(raster_sd, digits=3), file.path(paste0(save_dir, "/aggregated_VWC/sd.tif")), overwrite=TRUE)
    writeRaster(round(raster_median, digits=3), file.path(paste0(save_dir, "/aggregated_VWC/median.tif")), overwrite=TRUE)
    xml_names <- list.files(paste0(save_dir, "/aggregated_VWC"), pattern = "xml", full.names = T)
    unlink(xml_names, recursive=TRUE)

    write.csv(as.data.frame(round(raster_min, digits=3), xy=T, na.rm=T), paste0(save_dir, "/", as.numeric(as.Date(end_date)-as.Date(start_date))+1, "_day_min.csv"))
    write.csv(as.data.frame(round(raster_max, digits=3), xy=T, na.rm=T), paste0(save_dir, "/", as.numeric(as.Date(end_date)-as.Date(start_date))+1, "_day_max.csv"))
    write.csv(as.data.frame(round(raster_mean, digits=3), xy=T, na.rm=T), paste0(save_dir, "/", as.numeric(as.Date(end_date)-as.Date(start_date))+1, "_day_mean.csv"))
    write.csv(as.data.frame(round(raster_sd, digits=3),  xy=T, na.rm=T), paste0(save_dir, "/", as.numeric(as.Date(end_date)-as.Date(start_date))+1, "_day_sd.csv"))
    write.csv(as.data.frame(round(raster_median, digits=3),  xy=T, na.rm=T), paste0(save_dir, "/", as.numeric(as.Date(end_date)-as.Date(start_date))+1, "_day_median.csv"))
  } else {
    return("Files already exist.")
  }
}



#' Aggregation visualization
#'
#' Visualizes the output maps saved by the two aggregation functions.
#'
#' @param depth numeric. Defines the depth of the VWC values to visualize.
#' @param start_date character. Defines the start date of the aggregation period.
#' @param date character. Defines the starting date of an interval saved by aggregation_summary().
#' @param frequency numeric. Defines the length of the aggregation time interval (in days).
#' @param end_date character. Defines the end date of an aggregation period saved by aggregation_by_date().
#' @param project character. Defines the name of the folder where the maps are saved. If not defined, the function will look for maps in the working directory.
#'
#' @return a leaflet map with layers of the basic statistics of the interval.
#'
#' @import leaflet raster
#'
#' @export
#'
#' @examples
#' \dontrun{ plot_aggregated_VWC(5, "2020-06-15", "2020-08-15", date="2020-06-29", frequency=14, project="WI_region")
#'           plot_aggregated_VWC(5, "2020-06-15", "2020-08-15", project="WI_region") }
plot_aggregated_VWC <- function(depth, start_date, end_date, date=NULL, frequency=NULL, project=NULL){
  setwd(getwd())
  WD <- getwd()
  if(!is.null(project)){
    if (!file.exists(project)){
      return("Project folder doesn't exist.")
    } else {
      WD <- paste0(WD, "/", project)
    }
  }

  if (!is.null(frequency)){
    if(!dir.exists(paste0(WD, paste0("/VWC_aggregation/", start_date, "_", end_date, "/", frequency,"days_", depth, "cm")))) return("Files don't exist.")
    WD <- paste0(WD, paste0("/VWC_aggregation/", start_date, "_", end_date, "/", frequency,"days_", depth, "cm/aggregated_VWC"))
    if(is.null(date)) return("No date specified.")
    xml_names <- list.files(WD, pattern = "xml", full.names = T)
    unlink(xml_names, recursive=TRUE)
    rasters <- stack(paste0(WD, "/", paste0(c("mean_", "min_", "median_", "max_","sd_"), date, ".tif")))
  } else if (is.null(frequency) && !is.null(end_date)){
    if(!dir.exists(paste0(WD, paste0("/VWC_aggregation/", start_date,"_", end_date,"/", depth,"cm")))) return("Files don't exist.")
    WD <- paste0(WD, paste0("/VWC_aggregation/", start_date,"_", end_date, "/", depth, "cm/aggregated_VWC"))
    xml_names <- list.files(WD, pattern = "xml", full.names = T)
    unlink(xml_names, recursive=TRUE)
    rasters <- stack(paste0(WD, "/", paste0(c("mean", "min", "median", "max","sd"), ".tif")))
  } else if (is.null(frequency) && is.null(end_date)){
    return("Please specify the files to visualize.")
  }

  m <- c("Mean", "Min", "Median", "Max","SD")

  map <- leaflet() %>%
    addProviderTiles("Esri.WorldImagery")

  print("Plotting temporal VWC...")
  pb <- txtProgressBar(min=0, max=length(m), style=3)
  for (i in seq(length(m))){
    vwc_val <- rasters[[i]]
    min_val <- min(values(vwc_val),na.rm=T)
    max_val <- max(values(vwc_val),na.rm=T)
    pal1 <- colorNumeric(
      brewer.pal(n = 5, name = "RdYlBu"),
      domain = c(min_val, max_val),
      na.color = "transparent")
    map <- map %>%
      addRasterImage(vwc_val, colors=pal1, opacity = 0.7, group=m[i], project=FALSE) %>%
      addLegend("bottomright", opacity=0.8, pal = pal1, values = values(vwc_val),
                bins = seq(min_val, max_val, (max_val-min_val)/6),
                title=m[i], group=m[i])
    setTxtProgressBar(pb,i)
  }
  close(pb)

  map <- map %>%
    addLayersControl(
      baseGroups = m,
      options = layersControlOptions(collapsed = FALSE)
    )%>%
    htmlwidgets::onRender("
        function(el, x) {
          var updateLegend = function () {
              var selectedGroup = document.querySelectorAll('input:checked')[0].nextSibling.innerText.substr(1);

              document.querySelectorAll('.legend').forEach(a => a.hidden=true);
              document.querySelectorAll('.legend').forEach(l => {
                if (l.children[0].children[0].innerText == selectedGroup) l.hidden=false;
              });
          };
          updateLegend();
          this.on('baselayerchange', e => updateLegend());
        }") # https://stackoverflow.com/questions/52393310/is-it-possible-to-switch-between-multiple-legends-when-switching-between-base-gr

  return(map)
}
