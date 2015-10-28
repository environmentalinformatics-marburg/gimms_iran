################################################################################
### global settings
################################################################################

## required packages
library(gimms)
library(rworldmap)
library(doParallel)
library(remote)
library(Kendall)

# library(devtools)
# install_github('dutri001/bfastSpatial')
library(bfastSpatial)

## initialize parallel backend
cores <- detectCores() - 1
cl <- makeCluster(cores)
registerDoParallel(cl)


################################################################################
### data processing
################################################################################

## location of gimms raw data (jan 1982 to dec 2013)
ch_dir_extdata <- "/media/fdetsch/dev/gimms_iran/data/"

gimms_files <- rearrangeFiles(dsn = ch_dir_extdata, full.names = TRUE)
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


################################################################################
## remove seasonal signal
################################################################################

## calculate long-term bi-monthly means
lst_gimms_means <- 
  foreach(i = 1:24, .packages = c("raster", "rgdal")) %dopar% {
    
    # layers corresponding to current period (e.g. '81jul15a')
    id <- seq(i, nlayers(gimms_rasters_crop), 24)
    rst_gimms_tmp <- gimms_rasters_crop[[id]]
    
    # calculate long-term mean of current period (e.g. for 1981-2013 'jul15a')
    calc(rst_gimms_tmp, fun = mean, na.rm = TRUE,
         filename = paste0(ch_dir_extdata, "longterm_means/mean_", formatC(i, width = 2, flag = "0")),
         format = "GTiff", overwrite = TRUE)
  }

rst_gimms_means <- stack(lst_gimms_means)

## replicate bi-monthly 'gimms_raster_means' to match up with number of layers of
## initial 'gimms_raster_agg' (as `foreach` does not support recycling!)
lst_gimms_means <- replicate(nlayers(gimms_rasters_crop) / nlayers(rst_gimms_means),
                              rst_gimms_means)
rst_gimms_means <- stack(lst_gimms_means)
rst_gimms_means <- stack(rst_gimms_means, rst_gimms_means[[1:12]])

## subtract long-term mean from bi-monthly values
lst_gimms_dsn <-
  foreach(i = 1:nlayers(gimms_rasters_crop),
          .packages = c("raster", "rgdal")) %dopar% {
            
            overlay(gimms_rasters_crop[[i]], rst_gimms_means[[i]],
                    fun = function(x, y) {x - y},
                    filename = paste0(ch_dir_extdata, "dsn/DSN_",
                                      names(gimms_rasters_crop[[i]]), ".tif"),
                    format = "GTiff", overwrite = TRUE)
            
          }

rst_gimms_dsn <- stack(lst_gimms_dsn)
fls_gimms_dsn <- list.files(paste0(ch_dir_extdata, "dsn"), 
                            pattern = "^DSN_.*VI3g.tif", full.names = TRUE)

mat_gimms_dsn <- as.matrix(rst_gimms_dsn)
plot(ts(mat_gimms_dsn[1293, ], start = c(1981, 12), end = c(2013, 11), frequency = 24), type = "l")

plot(mat_gimms_dsn[1293, ], type = "l")

################################################################################
### apply mann-kendall trend test (Mann, 1945) on a pixel basis
################################################################################

## custom function that returns significant values of tau only
significantTau <- function(x, p = 0.001) {
  mk <- Kendall::MannKendall(x)
  # reject value of tau if p >= 0.001
  if (mk$sl >= p) {
    return(NA)
    # keep value of tau if p < 0.001
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

################################################################################
## mann-kendall trend test (p < 0.001) per season
################################################################################

fls_gimms_dsn <- rearrangeFiles(dsn = paste0(ch_dir_extdata, "/dsn"), 
                                pattern = "^DSN_.*VI3g.tif", full.names = TRUE, 
                                pos = c(4+8, 6+8, 11+8))

id_82 <- grep("81dec15a", fls_gimms_dsn)
id_13 <- grep("13nov15b", fls_gimms_dsn)
fls_gimms_dsn <- fls_gimms_dsn[id_82:id_13]
rst_gimms_dsn <- stack(fls_gimms_dsn)

systime_locale <- Sys.getlocale(category = "LC_TIME")
if (Sys.info()[["sysname"]] == "Windows") {
  invisible(Sys.setlocale(category = "LC_TIME", locale = "C"))
} else {
  invisible(Sys.setlocale(category = "LC_TIME", locale = "en_US.UTF-8"))
}

yrmn <- zoo::as.yearmon(substr(basename(fls_gimms_dsn), 12, 16), "%y%b")
indices <- rep(1:(length(yrmn) / 6), each = 6)
rst_gimms_dsn <- stackApply(rst_gimms_dsn, indices = indices, fun = mean)

lst_ssn_trends <- 
  foreach(i = list(1, 2, 3, 4), j = list("DJF", "MAM", "JJA", "SON")) %dopar% {
    rst_ssn <- rst_gimms_dsn[[seq(i, nlayers(rst_gimms_dsn), 4)]]
    
    calc(rst_ssn, fun = function(x) significantTau(x, p = 0.001),
         filename = paste0(ch_dir_extdata, "out/gimms_mk001_", j, "_8213"),
         format = "GTiff", overwrite = TRUE)
}

rst_ssn_trends <- stack(lst_ssn_trends)

library(RColorBrewer)
cols <- colorRampPalette(brewer.pal(5, "BrBG"))
spplot(mask(rst_ssn_trends, spy_iran), col.regions = cols(1000), at = seq(-.85, .85, .01)) +
  latticeExtra::layer(sp.polygons(spy_iran, col = "black"))

################################################################################
### breakpoint detection (see http://www.loicdutrieux.com/bfastSpatial/)
################################################################################

gimms_dates <- monthlyIndices(gimms_files, to_date = TRUE, 
                              format = "%Y-%m-%d")
gimms_dates <- as.Date(gimms_dates)
gimms_rasters_bfm <- bfmSpatial(gimms_rasters_crop, dates = gimms_dates, 
                                mc.cores = cores)

## deregister parallel backend
stopCluster(cl)
