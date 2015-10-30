################################################################################
### global settings
################################################################################

## required packages
library(gimms)
library(rworldmap)
library(doParallel)
library(remote)
library(Kendall)
library(RColorBrewer)

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
ch_dir_extdata <- "/media/fdetsch/XChange/gimms_iran/data/"

gimms_files <- rearrangeFiles(dsn = ch_dir_extdata, full.names = TRUE)
gimms_rasters <- stack(gimms_files)

## crop global gimms images
if (!file.exists("data/IRN_adm0.rds"))
  download.file("http://biogeo.ucdavis.edu/data/gadm2.7/rds/IRN_adm0.rds",
                destfile = "data/IRN_adm0.rds")
spy_iran <- readRDS("data/IRN_adm0.rds")
suppressWarnings(proj4string(spy_iran) <- "+init=epsg:4326")

gimms_files_crop <- paste0(ch_dir_extdata, "crp/CRP_")
gimms_list_crop <- foreach(i = 1:nlayers(gimms_rasters), 
                              overwrite = rep(TRUE, nlayers(gimms_rasters)), 
                             .packages = c("raster", "rgdal")) %dopar% {
                               
  gimms_file_crop <- paste0(gimms_files_crop, names(gimms_rasters)[i], ".tif")                             
             
  if (!file.exists(gimms_file_crop) | overwrite) {
    crop(gimms_rasters[[i]], spy_iran, snap = "out", filename = gimms_file_crop, 
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
  foreach(i = 1:24, j = month.abb[rep(c(7:12, 1:6), each = 2)], 
          k = rep(c("15a", "15b"), 12), .packages = c("raster", "rgdal")) %dopar% {
    
    # layers corresponding to current period (e.g. '81jul15a')
    id <- seq(i, nlayers(gimms_rasters_crop), 24)
    rst_gimms_tmp <- gimms_rasters_crop[[id]]
    
    # calculate long-term mean of current period (e.g. for 1981-2013 'jul15a')
    calc(rst_gimms_tmp, fun = mean, na.rm = TRUE,
         filename = paste0(ch_dir_extdata, "longterm_means/mean_", j, k),
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


################################################################################
### apply mann-kendall trend test (Mann, 1945) on a pixel basis
################################################################################

## 1982-2013
fls_gimms_dsn <- rearrangeFiles(dsn = paste0(ch_dir_extdata, "/dsn"), 
                                pattern = "^DSN_.*VI3g.tif", full.names = TRUE, 
                                pos = c(4, 6, 11) + 8)

id_82 <- grep("82jan15a", fls_gimms_dsn)
id_13 <- grep("13dec15b", fls_gimms_dsn)
fls_gimms_dsn <- fls_gimms_dsn[id_82:id_13]
rst_gimms_dsn <- stack(fls_gimms_dsn)

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

cols <- colorRampPalette(brewer.pal(5, "BrBG"))

lst_gimms_trend <- foreach(p = c(0.05, 0.01, 0.001), 
                           .packages = c("raster", "rgdal", "Kendall", "latticeExtra")) %do% {
  
#   rst_gimms_trend <- calc(rst_gimms_dsn,
#                           fun = function(x) significantTau(x, p = p),
#                           filename = paste0(ch_dir_extdata, "out/gimms_mk", gsub("0\\.", "", p), "_8213"),
#                           format = "GTiff", overwrite = TRUE)
  
  # reimport
  file_in <- paste0(ch_dir_extdata, "out/gimms_mk", gsub("0\\.", "", p), "_8213.tif")
  rst_gimms_trend <- raster(file_in)
  
  png(paste0(ch_dir_extdata, "vis/gimms_mk", gsub("0\\.", "", p), ".png"), 
      width = 9, height = 9, units = "cm", res = 500)
  p <- spplot(mask(rst_gimms_trend, spy_iran), col.regions = cols(1000),
              scales = list(draw = TRUE), maxpixels = 100000, 
              at = seq(-.6, .6, .01)) +
    latticeExtra::layer(sp.polygons(spy_iran, col = "black"))
  print(p)
  dev.off()
  
  return(rst_gimms_trend)
}


################################################################################
## mann-kendall trend test (p < 0.001) per month
################################################################################

## monthly composites
# indices <- rep(1:(length(fls_gimms_dsn)/2), each = 2)
# rst_gimms_dsn_mnth <- stackApply(rst_gimms_dsn, indices, fun = mean)
# saveRDS(rst_gimms_dsn_mnth, file = "inst/extdata/dsn_mnth.rds")
rst_gimms_dsn_mnth <- readRDS("inst/extdata/dsn_mnth.rds")

for (p in c(0.05, 0.01, 0.001)) {
  
  #   lst_mnth_trends <-
  #     foreach(i = 1:12, j = month.abb,
  #             .packages = c("raster", "rgdal", "Kendall"),
  #             .export = "significantTau") %dopar% {
  #               
  #               rst_mnth <- rst_gimms_dsn_mnth[[seq(i, nlayers(rst_gimms_dsn_mnth), 12)]]
  #               
  #               calc(rst_mnth, fun = function(x) significantTau(x, p = p),
  #                    filename = paste0(ch_dir_extdata, "out/monthly/gimms_mk", gsub("0\\.", "", p), "_", j, "_8213"),
  #                    format = "GTiff", overwrite = TRUE)
  #             }
  #   
  #   rst_mnth_trends <- stack(lst_mnth_trends)
  
  # reimport
  lst_mnth_trends <- lapply(month.abb, function(i) {
    file_in <- paste0(ch_dir_extdata, paste0("out/monthly/gimms_mk", gsub("0\\.", "", p), "_", i, "_8213.tif"))
    raster(file_in)
  })
  rst_mnth_trends <- stack(lst_mnth_trends)
  
  cols <- colorRampPalette(brewer.pal(5, "BrBG"))
  png(paste0(ch_dir_extdata, "vis/gimms_mk", gsub("0\\.", "", p), "_monthly.png"), 
      width = 20, height = 24, units = "cm", res = 500)
  p <- spplot(mask(rst_mnth_trends, spy_iran), col.regions = cols(1000), 
              at = seq(-.85, .85, .01), layout = c(3, 4)) +
    latticeExtra::layer(sp.polygons(spy_iran, col = "black"))
  print(p)
  dev.off()
}

################################################################################
## mann-kendall trend test (p < 0.001) per season
################################################################################

fls_gimms_dsn <- rearrangeFiles(dsn = paste0(ch_dir_extdata, "/dsn"), 
                                pattern = "^DSN_.*VI3g.tif", full.names = TRUE, 
                                pos = c(4, 6, 11) + 8)

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
rst_gimms_dsn_ssn <- stackApply(rst_gimms_dsn, indices = indices, fun = mean)
saveRDS(rst_gimms_dsn_ssn, file = "inst/extdata/dsn_ssn.rds")
rst_gimms_dsn_ssn <- readRDS("inst/extdata/dsn_ssn.rds")

lst_ssn_trends <- foreach(p = c(0.05, 0.01, 0.001)) %do% {
  lst_ssn_trend <- 
    foreach(i = list(1, 2, 3, 4), j = list("DJF", "MAM", "JJA", "SON")) %dopar% {
      rst_ssn <- rst_gimms_dsn_ssn[[seq(i, nlayers(rst_gimms_dsn_ssn), 4)]]
      
      calc(rst_ssn, fun = function(x) significantTau(x, p = p),
           filename = paste0(ch_dir_extdata, "out/seasonal/gimms_mk", gsub("0\\.", "", p), "_", j, "_8213"),
           format = "GTiff", overwrite = TRUE)
    }
  
  rst_ssn_trends <- stack(lst_ssn_trend)
  
  png(paste0(ch_dir_extdata, "vis/gimms_mk", gsub("0\\.", "", p), "_ssn.png"), 
      width = 14, height = 14, units = "cm", res = 500)
  p <- spplot(mask(rst_ssn_trends, spy_iran), col.regions = cols(1000), at = seq(-.85, .85, .01)) +
    latticeExtra::layer(sp.polygons(spy_iran, col = "black"))
  print(p)
  dev.off()
  
  return(rst_ssn_trends)
}
  
    
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
