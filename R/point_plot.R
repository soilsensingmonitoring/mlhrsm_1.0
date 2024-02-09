#' Point Extraction
#'
#' Extracts VWC values of specified points from the predicted regional maps.
#' The point map should be saved within project folder for function access.
#'
#' @param site_  character. Defines the file name that stores the location of the points (in shapefile format).
#' @param depth numeric. Indicating the depth of the map to extract.
#' @param percentile logical. If TRUE, the upper and lower bounds of the 90\% Confidence Interval of the prediction will be included.
#' @param project character. Defines the name of the folder where the maps are saved. If not defined, the function will look for maps in the working directory.
#' @param filename character. Defines the file name where the predicted data is stored (in csv format). Default as VWC_point_data.csv.
#'
#' @return None
#'
#' @import rgdal tidyverse raster
#'
#' @export
#'
#' @examples
#' \dontrun{ point_extraction("WI_SM.shp", TRUE, project="WI_region") }
point_extraction <- function(site_, percentile=FALSE, project=NULL,
                             filename="VWC_point_data.csv") {
  setwd(getwd())
  WD <- getwd()
  if(!is.null(project)){
    if (!file.exists(project)){
      print("Project folder doesn't exist.")
    } else {
      WD <- paste0(WD, "/", project)
    }
  }

  if (file.exists(paste0(WD, "/", filename))) writeLines("***Extracted point file already exist.***")
  Sys.sleep(2)

  xml_names <- list.files(paste0(WD, "/VWC"), pattern = "xml", full.names = T)
  unlink(xml_names, recursive=TRUE)

  validation.pt <- readOGR(WD, unlist(str_split(site_,"\\."))[1]) # `point` might include ".shp"

  vwc_5_mean <- list.files(paste0(WD,"/VWC"), "VWC_5_mean", full.names = T )
  vwc_5_sd <- list.files(paste0(WD,"/VWC"), "VWC_5_sd", full.names = T )
  map_5_mean <- lapply(vwc_5_mean, raster)
  map_5_sd <- lapply(vwc_5_sd, raster)
  vwc_100_mean <- list.files(paste0(WD,"/VWC"), "VWC_100_mean", full.names = T )
  vwc_100_sd <- list.files(paste0(WD,"/VWC"), "VWC_100_sd", full.names = T )
  map_100_mean <- lapply(vwc_100_mean, raster)
  map_100_sd <- lapply(vwc_100_sd, raster)

  percentile <- FALSE
  if (length(list.files(paste0(WD,"/VWC"), "VWC_5_lower"))>0) percentile <- TRUE

  if (percentile){
    vwc_5_lower <- list.files(paste0(WD,"/VWC"), "VWC_5_lower", full.names = T )
    vwc_5_upper <- list.files(paste0(WD,"/VWC"), "VWC_5_upper", full.names = T )
    map_5_lower <- lapply(vwc_5_lower, raster)
    map_5_upper <- lapply(vwc_5_upper, raster)
    vwc_100_lower <- list.files(paste0(WD,"/VWC"), "VWC_100_lower", full.names = T )
    vwc_100_upper <- list.files(paste0(WD,"/VWC"), "VWC_100_upper", full.names = T )
    map_100_lower <- lapply(vwc_100_lower, raster)
    map_100_upper <- lapply(vwc_100_upper, raster)
  }

  dates <- gsub(".*_(\\d*-\\d*-\\d*).tif", "\\1", vwc_5_mean)

  ## Extract mapped soil moisture to the validation points
  VWC_value <- as.data.frame(validation.pt)
  print("Extracting site VWC...")
  pb <- txtProgressBar(0, length(dates), style=3)
  for (i in 1:length(dates))
  {
    temp_5_mean <- raster::extract(map_5_mean[[i]], validation.pt)
    temp_5_sd <- raster::extract(map_5_sd[[i]], validation.pt)
    temp_100_mean <- raster::extract(map_100_mean[[i]], validation.pt)
    temp_100_sd <- raster::extract(map_100_sd[[i]], validation.pt)
    VWC_value <- cbind(VWC_value, as.data.frame(round(temp_5_mean, digits=3)), as.data.frame(round(temp_5_sd, digits=3)), as.data.frame(round(temp_100_mean, digits=3)), as.data.frame(round(temp_100_sd, digits=3)))

    if (percentile){
      temp_5_lower <- raster::extract(map_5_lower[[i]], validation.pt)
      temp_5_upper <- raster::extract(map_5_upper[[i]], validation.pt)
      temp_100_lower <- raster::extract(map_100_lower[[i]], validation.pt)
      temp_100_upper <- raster::extract(map_100_upper[[i]], validation.pt)
      VWC_value <- cbind(VWC_value, as.data.frame(round(temp_5_lower, digits=3)), as.data.frame(round(temp_5_upper, digits=3)), as.data.frame(round(temp_100_lower, digits=3)), as.data.frame(round(temp_100_upper, digits=3)))
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)

  if(!percentile){
    colnames(VWC_value) <- c("ID", "Longitude", "Latitude", paste0(rep(dates,each=4),rep(c(".5Mean",".5SD",".100Mean",".100SD"),times=length(dates))))
  } else {
    colnames(VWC_value) <- c("ID", "Longitude", "Latitude", paste0(rep(dates,each=8),
                                                                   rep(c(".5Mean",".5SD", ".100Mean",".100SD", ".5Lower", ".5Upper", ".100Lower", ".100Upper"),times=length(dates))))
  }

  VWC_value <- VWC_value %>%
    pivot_longer(4:length(VWC_value),names_to="Date",values_to="value") %>%
    mutate(
      mode = gsub("(.*)\\.(.*)","\\2",Date),
      Date = as.Date(gsub("(.*)\\.(.*)","\\1",Date))
    ) %>%
    pivot_wider(names_from="mode",values_from=value)

  names(VWC_value)[names(VWC_value) == "5Mean"] <- "VWC_5_mean_pts"
  names(VWC_value)[names(VWC_value) == "5SD"] <- "VWC_5_sd_pts"
  names(VWC_value)[names(VWC_value) == "100Mean"] <- "VWC_100_mean_pts"
  names(VWC_value)[names(VWC_value) == "100SD"] <- "VWC_100_sd_pts"

  if (percentile) {
    names(VWC_value)[names(VWC_value) == "5Lower"] <- "VWC_5_lower_pts"
    names(VWC_value)[names(VWC_value) == "5Upper"] <- "VWC_5_upper_pts"
    names(VWC_value)[names(VWC_value) == "100Lower"] <- "VWC_100_lower_pts"
    names(VWC_value)[names(VWC_value) == "100Upper"] <- "VWC_100_upper_pts"
  }

  write.csv(VWC_value, paste0(WD, "/", filename), row.names=F)

}



#' Point Variation
#'
#' Plotting function for the change in VWC at specific points over time.
#'
#' @param depth numeric. Indicating the depth of the plotted VWC.
#' @param percentile logical. If TRUE, the upper and lower bounds of the 90\% Confidence Interval of the prediction will be included.
#' @param project character. Defines the name of the folder where the maps are saved. If not defined, the function will look for maps in the working directory.
#' @param site_ character or a vector of characters. Defines the specific points to include in the plot. If not defined, all points will be plotted (for site number less than or equal to 12). If site number exceeds 12, 12 points will be chosen at random.
#' @param filename character. Defines the file name where the predicted data is stored (in csv format). Default as VWC_point_data.csv.
#'
#' @return a list of plots. The first plot is a leaflet plot showing the locations of the points; the second plot is a time series plot showing the variation of each point’s predicted VWC overtime
#'
#' @import tidyverse leaflet latex2exp
#'
#' @export
#'
#' @examples
#' \dontrun{ plots <- site_variation(5, TRUE, project="WI_point", c("S13", "S49"))
#'  plots[[1]]
#'  plots[[2]]}
site_variation <- function(depth, percentile=FALSE, project=NULL, site_=NULL,
                           filename="VWC_point_data.csv") {
  setwd(getwd())
  WD <- getwd()
  if(!is.null(project)){
    if (!file.exists(project)){
      print("Project folder doesn't exist.")
    } else {
      WD <- paste0(WD, "/", project)
    }
  }

  data <- read_csv(paste0(WD, "/",filename))%>%
    dplyr::select(Date, ID, Longitude, Latitude, contains(paste0("VWC_", depth, "_")))%>%
    rename(
      "Mean"=paste0("VWC_", depth, "_mean_pts"),
      "SD"=paste0("VWC_", depth, "_sd_pts")
    )

  if (percentile) {
    data <- data %>%
      rename(
        "Lower"=paste0("VWC_", depth, "_lower_pts"),
        "Upper"=paste0("VWC_", depth, "_upper_pts"),
      )
  }

  if (!is.null(site_)){
    data <- data %>%
      filter(ID %in% site_)
  }

  p1 <- leaflet() %>%
    addProviderTiles("Esri.WorldImagery") %>%
    addCircleMarkers(lat = unique(data$Latitude),
                     lng = unique(data$Longitude),
                     label = unique(data$ID),
                     labelOptions = labelOptions(noHide = T))

  sites <- unique(data$ID)
  if(length(sites)>12){
    sites <- sample(sites, 12)
  }

  p2 <- ggplot(data %>% filter(ID %in% sites))+
    geom_line(aes(Date, Mean, linetype="a", col="a"))+
    geom_ribbon(aes(Date, ymin = Mean - SD,
                    ymax = Mean + SD, fill="sd"), alpha=0.3)+
    facet_wrap(~~factor(ID, stringr::str_sort(unique(ID), numeric = TRUE)))+
    labs(y=TeX("VWC (m^3m^{-3})"), x="Date", title="Site-specific VWC variations over time")+
    theme_bw()+
    # scale_linetype_manual("", "mean", labels = "Mean", values=1)+
    scale_color_manual(values="grey16")+
    scale_fill_manual("","sd",labels="Standard Deviation",values="grey70")+
    scale_x_date(date_labels =  "%y-%m-%d")+
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1,size=8))+
    labs(linetype=NULL, col=NULL)

  if(percentile==TRUE){
    p2 <- p2 +
      geom_line(aes(Date, Lower,linetype="b", col="b"))+
      geom_line(aes(Date, Upper,linetype="b", col="b"))+
      scale_linetype_manual(labels = c("Mean", "90% Confidence Interval"),values=c("solid", "dashed"))+
      scale_color_manual(labels = c("Mean", "90% Confidence Interval"),values=c("grey16", "grey52"))
  }

  return(list(p1, p2))
}
