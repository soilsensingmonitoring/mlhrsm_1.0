# mlhrsm_1.0

Machine-Learning based High-Resolution Soil Moisture mapping and spatial analysis tool - Version 1.0

Authors: Yuliang Peng (peng68@wisc.edu), Jingyi Huang (jhuang426@wisc.edu)

For package support or future collaboration, please contact Dr. Jingyi Huang by email or join our Google Group (mlhrsm@googlegroups.com).

Package installation instructions 

1. Install the latest version RTools (RTools 4.2 or the version that is compatible with the user’s R console, https://cran.r-project.org/bin/windows/Rtools/).

2. Install the following dependency R packages.
> install.packages(c('raster', 'rgee', 'sf', 'tidyverse', 'viridis', 'FedData', 'RColorBrewer', 'caret', 'chillR', 'leaflet', 'hydroGOF', 'quantregForest', 'randomForest', 'reshape2', 'sp', 'lubridate', 'geojsonio', 'stars', 'Rcpp', 'fastmap', 'digest', 'fs', 'stringi', 'cachem', 'htmltools', 'curl', 'ps', 'processx'))

Since `rgdal` is no longer available on CRAN, we need to install it through `devtools`.
> install.packages("devtools")
> library(devtools)
> install_version("rgdal", version = "1.6-7", repos = "http://cran.us.r-project.org") 
> library(rgdal)
> options("rgdal_show_exportToProj4_warnings"="none")   ## Please run this line so that rgdal package can still be used without confusion with sf and terra packages

3. Install R package mlhrsm. The users can install it from GitHub.
> install.packages("R.rsp")
> devtools::install_github("soilsensingmonitoring/mlhrsm_1.0", build_vignettes=T)

4. Set up Google Earth Engine account, project, and API. 
First, all users need to create a free Google Earth Engine account (https://earthengine.google.com/signup/).
Second, install gcloud CLI before downloading maps from Google Earth Engine (https://dl.google.com/dl/cloudsdk/channels/rapid/GoogleCloudSDKInstaller.exe). 
Third, create a project on the Google Earth Account for future use. After installing the gcloud CLI, if a CMD window pops out (when the user enables configuration of gcloud) to ask the user to connect gcloud CLI, select “Y” to log in. Then a web page will appear with a message saying “Google Cloud SDK wants to access your Google Account”; select Allow and go back to the CMD where the system asks the user to “Pick cloud project to use.” Select the project the user wants to use, and close CMD. Lastly, relaunch R software and install Google Earth Engine API in the R environment.
In the future, if the user wants to reset the gcloud, please follow this page for detailed instructions. https://cloud.google.com/sdk/docs/configurations

> library(mlhrsm)

> ee_Initialize("Your email address", drive=T) # insert your email address

If it’s the first time the user uses ee_Initialize() on the computer, R will print downloading and installation messages when preparing for the initialization. Select “Y” when R asks to install Miniconda. If the computer does not have the Python package "earthengine-api" installed, an error message will appear and the user should run the following command line to install it.  
	
> .rs.restartR() ## If this does not work, please restart the R session manually

> ee_install()

Then R will ask the user to store environment variables EARTHENGINE_PYTHON and EARTHENGINE_ENV in the .Renviron file to use Python path in future sessions. Type “Y” to continue, and restart R session when prompted to do so after installation is completed. Run ee_Initialize("Your email address", drive=T) again. A new window will pop up in the browser saying “Google Earth Engine Authenticator wants to access your Google Account”, then select Allow to allow the local R environment to connect to the user’s Google Earth Engine. If successful, the user will see the following messages in the R console. The user can now access the maps in Earth Engine from the local R environment and download them to the user’s Google Drive. 
	
Fetching credentials using gcloud

Successfully saved authorization token 
	
###################### Demo ####################
To see some examples, please download the files in the "Demo" folder and run "Demo.R" file.
