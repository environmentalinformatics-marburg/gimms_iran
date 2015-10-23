################################################################################
### global settings
################################################################################

## required packages
library(gimms)
library(rworldmap)
library(doParallel)
library(remote)
library(Kendall)

## initialize parallel backend
cores <- detectCores() - 1
cl <- makeCluster(cores)
registerDoParallel(cl)


################################################################################
### data processing
################################################################################

## location of gimms raw data (jan 1982 to dec 2013)
ch_dir_extdata <- "/media/fdetsch/XChange/gimms_iran/data/"

gimms_files <- rearrangeFiles(dsn = ch_dir_extdata, full.names = TRUE)[-(1:12)]
gimms_rasters <- stack(gimms_files)

## crop global gimms images
data("countriesCoarse")
spy_iran <- subset(countriesCoarse, ADMIN == "Iran")

gimms_files_crop <- paste0(ch_dir_extdata, "crp/CRP_")
gimms_list_crop <- foreach(i = 1:nlayers(gimms_rasters), 
                              overwrite = rep(FALSE, nlayers(gimms_rasters)), 
                             .packages = c("raster", "rgdal")) %dopar% {
                               
  gimms_file_crop <- paste0(gimms_files_crop, names(gimms_rasters)[i], ".tif")                             
             
  if (!file.exists(gimms_file_crop) | overwrite) {
    crop(gimms_rasters[[i]], spy_iran, filename = gimms_file_crop, 
         overwrite = TRUE, format = "GTiff")
  } else {
    raster(gimms_file_crop)
  }
}

gimms_rasters_crop <- stack(gimms_list_crop)

## remove seasonal signal
gimms_rasters_deseason <- deseason(gimms_rasters_crop, 
                                   cycle.window = 24, use.cpp = TRUE)

gimms_files_deseason <- paste0(ch_dir_extdata, "dsn/DSN_", 
                               names(gimms_rasters_crop))
gimms_list_deseason <- foreach(i = 1:nlayers(gimms_rasters_deseason), 
                               .packages = c("raster", "rgdal")) %dopar%
  writeRaster(gimms_rasters_deseason[[i]], filename = gimms_files_deseason[i], 
              overwrite = TRUE, format = "GTiff")
gimms_rasters_deseason <- stack(gimms_list_deseason)

## apply mann-kendall trend test
kendallsTau <- function(x, p = 0.001) {
  
  # if not specified, set p to 1
  if (missing(p)) p <- 1
  
  mk <- MannKendall(x)

  if (mk$sl >= p) {
    return(NA) 
  } else {
    return(mk$tau)
  }
}

gimms_raster_trend <- overlay(gimms_rasters_deseason, fun = function(x) {
  cat("Processing cell no.", n, "of", ncell(gimms_rasters_deseason), ".\n")
  n <- n + 1
  kendallsTau(x, p = .001)
}, filename = paste0(ch_dir_extdata, "out/gimms_mk001_8213"), 
format = "GTiff", overwrite = TRUE)



## deregister parallel backend
stopCluster(cl)
