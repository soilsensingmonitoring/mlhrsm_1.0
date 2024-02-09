#' Validation time series
#'
#' Compares the predicted and measured data via time series plot with 90 percent Confidence Interval.
#'
#' @param test  character. Defines the file name of the validation (measured) data (in csv format).
#' @param depth numeric. Indicating the depth of the predicted and measured data.
#' @param project character. Defines the name of the folder where the maps are saved. If not defined, the function will look for maps in the working directory.
#' @param site_ character or a vector of characters. Defines the specific points to include in the plot. If not defined, all points will be plotted (for site number less than or equal to 12). If site number exceeds 12, 12 points will be chosen at random.
#' @param filename character. Defines the file name where the predicted data is stored (in csv format). Default as VWC_point_data.csv.
#'
#' @return a time series plot with predicted and measured data, as well as the 90 percent CI of the predicted data.
#'
#' @import tidyverse latex2exp
#' @importFrom lubridate ymd
#'
#' @export
#'
#' @examples
#' \dontrun{ point_CI("Huges_VWC.csv", 5, project="WI_point") }
point_CI <- function(test, depth, project=NULL, site_=NULL,
                     filename="VWC_point_data.csv"){
  setwd(getwd())
  WD <- getwd()
  if(!is.null(project)){
    if (!file.exists(project)){
      print("Project folder doesn't exist.")
    } else {
      WD <- paste0(WD, "/", project)
    }
  }

  if (!file.exists(paste0(WD,"/validation_", depth, ".csv"))) {
    # where to store test csv
    test_data <- read.csv(test) %>%
      mutate(Date = as.Date(ymd(Date)))
    ## Model performance
    predicted_data  = read.csv(paste0(WD, "/", filename)) %>%
      mutate(Date=as.Date(Date))
    compare <- predicted_data %>%
      dplyr::select(ID, Date, contains(paste0("VWC_", depth)))%>%
      full_join(test_data, by=c("ID", "Date")) %>%
      rename("measured" = paste0("VWC_", depth), "predicted" = paste0("VWC_", depth, "_mean_pts"),
             "sd" = paste0("VWC_", depth, "_sd_pts"), "lower" = paste0("VWC_", depth, "_lower_pts"), "upper" = paste0("VWC_", depth, "_upper_pts")) %>%
      drop_na()
    write.csv(compare, paste0(WD, "/validation_", depth, ".csv"), row.names=F)
  } else {
    compare <- read.csv(paste0(WD, "/validation_", depth, ".csv")) %>%
      mutate(Date=as.Date(Date))
  }

  if (!is.null(site_)) {
    compare <- compare %>%
      filter(ID==site_)
  }

  p <- ggplot(compare)+
    geom_ribbon(aes(Date, ymin = lower,
                    ymax = upper, fill="ci"),alpha=0.5)+
    geom_line(aes(Date, predicted, linetype="a"), col="grey16")+
    geom_line(aes(Date, measured, linetype="b"), col="grey16")+
    facet_wrap(~~factor(ID, stringr::str_sort(unique(ID), numeric = TRUE)))+
    theme_bw()+
    scale_fill_manual("","ci",labels="90% Confidence\n Interval",values="grey76")+
    scale_x_date(date_labels =  "%y-%m-%d")+
    scale_linetype_manual(values=c(1,2), labels=c("Predicted", "Measured"))+
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1,size=8))+
    labs(x="Date", y=TeX("VWC (m^3m^{-3})"),linetype="")

  return(p)
}



#' Validation
#'
#' Plots the predicted v. measured scatter plot.
#'
#' @param test  character. Defines the file name of the validation (measured) data (in csv format).
#' @param depth numeric. Indicating the depth of the predicted and measured data.
#' @param stats logical. If TRUE, the values of r-squared, RMSE, bias and KGE will be included in the diagnosis plot.
#' @param project character. Defines the name of the folder where the maps are saved. If not defined, the function will look for maps in the working directory.
#' @param filename character. Defines the file name where the predicted data is stored (in csv format). Default as VWC_point_data.csv.
#'
#' @return a list of validation values including r-squared, RMSE, NSE, KGE, bias, RPD, and RPIQ.
#'
#' @import tidyverse hydroGOF chillR
#' @importFrom lubridate ymd
#'
#' @export
#'
#' @examples
#' \dontrun{ point_performance("Huges_VWC.csv", 5, TRUE, project="WI_point") }
point_performance <- function(test, depth, stats=FALSE, project=NULL,
                              filename="VWC_point_data.csv"){
  setwd(getwd())
  WD <- getwd()
  if(!is.null(project)){
    if (!file.exists(project)){
      print("Project folder doesn't exist.")
    } else {
      WD <- paste0(WD, "/", project)
    }
  }

  # need to specify which variables are needed & their names
  if (!file.exists(paste0(WD,"/validation_", depth, ".csv"))) {
    # where to store test csv
    test_data <- read.csv(test) %>%
      mutate(Date = ymd(Date))
    ## Model performance
    predicted_data  = read.csv(paste0(WD, "/", filename)) %>%
      mutate(Date=as.Date(Date))
    compare <- predicted_data %>%
      dplyr::select(ID, Date, contains(paste0("VWC_", depth)))%>%
      full_join(test_data, by=c("ID", "Date")) %>%
      rename("measured" = paste0("VWC_", depth), "predicted" = paste0("VWC_", depth, "_mean_pts"),
             "sd" = paste0("VWC_", depth, "_sd_pts"), "lower" = paste0("VWC_", depth, "_lower_pts"), "upper" = paste0("VWC_", depth, "_upper_pts")) %>%
      drop_na()
    write.csv(compare, paste0(WD, "/validation_", depth, ".csv"), row.names=F)
  } else {
    compare <- read.csv(paste0(WD, "/validation_", depth, ".csv")) %>%
      mutate(Date=as.Date(Date))
  }

  # R-squared: https://stackoverflow.com/questions/40901445/function-to-calculate-r2-r-squared-in-r
  rsq <- cor(compare$measured, compare$predicted) ^ 2

  # RMSE: https://stackoverflow.com/questions/43123462/how-to-obtain-rmse-out-of-lm-result
  rmse <- sqrt(mean((compare$measured-compare$predicted)^2))

  # Bias
  bias <- mean(compare$predicted-compare$measured)

  # RPD & RPIQ: https://search.r-project.org/CRAN/refmans/chillR/html/RPD.html
  #             https://search.r-project.org/CRAN/refmans/chillR/html/RPIQ.html
  rpd <- RPD(compare$predicted, compare$measured)
  rpiq <- RPIQ(compare$predicted, compare$measured)

  # NSE: https://www.rdocumentation.org/packages/hydroGOF/versions/0.4-0/topics/NSE
  nse <- NSE(compare$predicted, compare$measured)

  # KGE: https://www.rdocumentation.org/packages/hydroGOF/versions/0.3-0/topics/KGE
  kge <- KGE(compare$predicted, compare$measured)

  plot(compare$measured,  compare$predicted, xlim = c(0, 0.4), ylim = c(0, 0.4), main= "Reduced model - Validation",
       xlab="Measured VWC", ylab="Predicted VWC")
  abline(0,1)

  if(stats){
    legend("topleft", cex=0.6, legend=c(paste("R-squared = ", round(rsq, 4)),
                                        paste("RMSE = ", round(rmse, 4)),
                                        paste("Bias = ", round(bias, 4)),
                                        paste("KGE = ", round(kge, 4))))
  }

  return(list("r_squared"=rsq, "RMSE"=rmse, "NSE"=nse, "KGE"=kge,
              "Bias"=bias, "RPD"=rpd, "RPIQ"=rpiq))
}
