################################################################################
### global settings
################################################################################

## working directory
setwd("/media/fdetsch/modis_data/gimms_iran/")

## required packages
lib <- c("Rsenal", "gimms", "doParallel", "remote", "Kendall", "RColorBrewer", 
         "plyr", "grid", "latticeExtra")
jnk <- sapply(lib, function(x) library(x, character.only = TRUE))

# library(devtools)
# install_github('dutri001/bfastSpatial')
# library(bfastSpatial)

## initialize parallel backend
cores <- detectCores()
supcl <- makeCluster(cores)
registerDoParallel(supcl)


################################################################################
### data processing
################################################################################

## gimms raw data
gimms_files <- rearrangeFiles(dsn = "data", full.names = TRUE)
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

# reimport
fls_gimms_crop <- rearrangeFiles(dsn = "data/crp", full.names = TRUE, 
                                 pattern = "CRP_.*.tif$", pos = c(4, 6, 8) + 4)
gimms_rasters_crop <- stack(fls_gimms_crop)

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
         filename = paste0("data/longterm_means/mean_", j, k),
         format = "GTiff", overwrite = TRUE)
  }

rst_gimms_means <- stack(lst_gimms_means)

# reimport
fls_gimms_means <- list.files("data/longterm_means", full.names = TRUE)
int_id <- do.call("c", lapply(month.abb[c(7:12, 1:6)], function(i) grep(i, fls_gimms_means)))
fls_gimms_means <- fls_gimms_means[int_id]
rst_gimms_means <- stack(fls_gimms_means)

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
                    filename = paste0("data/dsn/DSN_",
                                      names(gimms_rasters_crop[[i]]), ".tif"),
                    format = "GTiff", overwrite = TRUE)
            
          }

rst_gimms_dsn <- stack(lst_gimms_dsn)

# reimport()

################################################################################
### apply mann-kendall trend test (Mann, 1945) on a pixel basis
################################################################################

## 1982-2013
fls_gimms_dsn <- rearrangeFiles(dsn = "data/dsn", 
                                pattern = "^DSN_.*VI3g.tif", full.names = TRUE, 
                                pos = c(4, 6, 11) + 8)

id_82 <- grep("82jan15a", fls_gimms_dsn)
id_13 <- grep("13dec15b", fls_gimms_dsn)
fls_gimms_dsn <- fls_gimms_dsn[id_82:id_13]
rst_gimms_dsn <- stack(fls_gimms_dsn)

cols <- colorRampPalette(brewer.pal(5, "BrBG"))

lst_gimms_trend <- foreach(p = c(0.05, 0.01, 0.001)) %do% {
  
#   rst_gimms_trend <- calc(rst_gimms_dsn,
#                           fun = function(x, n = 1) {
#                             Rsenal::significantTau(x, p = p, prewhitening = TRUE, 
#                                                    conf.intervals = FALSE)
#                           }, filename = paste0("data/out/gimms_mk", gsub("0\\.", "", p), "_8213"),
#                           format = "GTiff", overwrite = TRUE)
  
  # reimport
  file_in <- paste0("data/out/gimms_mk", gsub("0\\.", "", p), "_8213.tif")
  rst_gimms_trend <- raster(file_in)
  
  png(paste0("data/vis/overall/gimms_mk", gsub("0\\.", "", p), ".png"), 
      width = 14, height = 14, units = "cm", res = 500)
  p <- spplot(mask(rst_gimms_trend, spy_iran), col.regions = cols(1000),
              scales = list(draw = TRUE), maxpixels = 100000, 
              at = seq(-.4, .4, .01)) +
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
# saveRDS(rst_gimms_dsn_mnth, file = "data/dsn_mnth.rds")
rst_gimms_dsn_mnth <- readRDS("data/dsn_mnth.rds")

lst_trends <- foreach(p = c(0.05, 0.01, 0.001)) %do% {
  
#     lst_mnth_trends <-
#       foreach(i = 1:12, j = month.abb, .packages = lib) %dopar% {
#                 
#                 rst_mnth <- rst_gimms_dsn_mnth[[seq(i, nlayers(rst_gimms_dsn_mnth), 12)]]
#                 
#                 calc(rst_mnth, fun = function(x, n = 1) {
#                   Rsenal::significantTau(x, p = p, prewhitening = TRUE, 
#                                          conf.intervals = FALSE)
#                 }, filename = paste0("data/out/monthly/gimms_mk", gsub("0\\.", "", p), "_", j, "_8213"),
#                      format = "GTiff", overwrite = TRUE)
#               }
#     
#     rst_mnth_trends <- stack(lst_mnth_trends)
  
  # reimport
  lst_mnth_trends <- lapply(month.abb, function(i) {
    file_in <- paste0("data/out/monthly/gimms_mk", gsub("0\\.", "", p), "_", i, "_8213.tif")
    raster(file_in)
  })
  rst_mnth_trends <- stack(lst_mnth_trends[c(12, 1:11)])
  
  cols <- colorRampPalette(brewer.pal(5, "BrBG"))
  png(paste0("data/vis/monthly/gimms_mk", gsub("0\\.", "", p), "_monthly.png"), 
      width = 20, height = 24, units = "cm", res = 500)
  p <- spplot(mask(rst_mnth_trends, spy_iran), col.regions = cols(1000), 
              at = seq(-.85, .85, .01), layout = c(3, 4)) +
    latticeExtra::layer(sp.polygons(spy_iran, col = "black"))
  print(p)
  dev.off()
  
  return(rst_mnth_trends)
}

################################################################################
## mann-kendall trend test (p < 0.001) per season
################################################################################

## 1982-2013
fls_gimms_dsn <- rearrangeFiles(dsn = "data/dsn", 
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

# yrmn <- zoo::as.yearmon(substr(basename(fls_gimms_dsn), 12, 16), "%y%b")
# indices <- rep(1:(length(yrmn) / 6), each = 6)
# rst_gimms_dsn_ssn <- stackApply(rst_gimms_dsn, indices = indices, fun = mean)
# saveRDS(rst_gimms_dsn_ssn, file = "data/dsn_ssn.rds")
rst_gimms_dsn_ssn <- readRDS("data/dsn_ssn.rds")

lst_ssn_trends <- foreach(p = c(0.05, 0.01, 0.001)) %do% {
  lst_ssn_trend <- 
    foreach(i = list(1, 2, 3, 4), j = list("DJF", "MAM", "JJA", "SON")) %dopar% {
      rst_ssn <- rst_gimms_dsn_ssn[[seq(i, nlayers(rst_gimms_dsn_ssn), 4)]]
      
      calc(rst_ssn, fun = function(x) {
        Rsenal::significantTau(x, p = p, prewhitening = TRUE, 
                               conf.intervals = FALSE)
        }, filename = paste0("data/out/seasonal/gimms_mk", gsub("0\\.", "", p), "_", j, "_8213"),
           format = "GTiff", overwrite = TRUE)
    }
  
  rst_ssn_trends <- stack(lst_ssn_trend)
  
  png(paste0("data/vis/seasonal/gimms_mk", gsub("0\\.", "", p), "_ssn.png"), 
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

## import bfast data
fls_bfast <- list.files("out", pattern = "bfast", full.names = TRUE)
lst_bfast <- lapply(fls_bfast, function(i) {
  do.call("rbind.fill", readRDS(i))
})
dat_bfast <- do.call("rbind.fill", lst_bfast)

## raster template
rst_template <- bfastTemplate()
rst_breaks_counts <- rst_breaks_maxtime <- rst_breaks_maxmagn <- rst_template

## number of breakpoints
rst_breaks_counts[] <- dat_bfast[, "bp_vt_ct"]
rst_breaks_counts <- mask(rst_breaks_counts, spy_iran)

col_counts <- brewer.pal(5, "YlGnBu")
p_counts <- spplot(rst_breaks_counts, col.regions = col_counts, 
                   at = seq(.5, 5.5, 1), scales = list(draw = TRUE), 
                   colorkey = FALSE) + 
  latticeExtra::layer(sp.polygons(spy_iran))

## most intense breakpoint
rst_breaks_maxmagn[] <- dat_bfast[, "bp_vt_maxmagn"]
rst_breaks_maxmagn <- mask(rst_breaks_maxmagn, spy_iran)

col_maxmagn <- colorRampPalette(brewer.pal(9, "BrBG"))

int_id <- which(rst_breaks_maxmagn[] > .25)
spt_maxval <- xyFromCell(rst_breaks_maxmagn, int_id, spatial = TRUE)
rst_breaks_maxmagn[int_id] <- NA

brks <- quantile(rst_breaks_maxmagn, seq(0, 1, length.out = 256))
p_maxmagn <- spplot(rst_breaks_maxmagn, col.regions = col_maxmagn(1000),
                    at = c(-.3, brks, .3), scales = list(draw = TRUE)) + 
  latticeExtra::layer(sp.polygons(spy_iran))

rst_breaks_maxtime[] <- dat_bfast[, "bp_vt_maxtime"]
rst_breaks_maxtime <- mask(rst_breaks_maxtime, spy_iran)

col_maxtime <- envinmrPalette(1000)[100:1000]
p_maxtime <- spplot(rst_breaks_maxtime, col.regions = col_maxtime, 
                    at = seq(1985, 2010, 1), scales = list(draw = TRUE), 
                    colorkey = FALSE) + 
  latticeExtra::layer(sp.polygons(spy_iran))

p_comb <- latticeCombineGrid(list(p_counts, p_maxtime, p_maxmagn), 
                             layout = c(3, 1))

# store
png("out/bfast_trend.png", width = 21, height = 13, units = "cm", res = 500)
grid.newpage()

vp0 <- viewport(x = 0, y = 0, width = 1, height = .9, 
                just = c("left", "bottom"), name = "vp_count")
pushViewport(vp0)
print(p_comb, newpage = FALSE)

# key left
upViewport()
vp1 <- viewport(x = .06, y = .72, width = 1/3, height = .1, 
                just = c("left", "bottom"), name = "vp_key1")
pushViewport(vp1)
draw.colorkey(key = list(col = col_counts, width = 1, height = .5,
                         at = seq(.5, 5.5, 1), 
                         space = "bottom"), draw = TRUE)
grid.text("Trend breakpoints", x = 0.5, y = 1.12, just = c("centre", "top"), 
          gp = gpar(font = 2, cex = .85))

# key middle
upViewport()
vp2 <- viewport(x = 1/3 + .02, y = .72, width = 1/3, height = .1, 
                just = c("left", "bottom"), name = "vp_key2")
pushViewport(vp2)
draw.colorkey(key = list(col = col_maxtime, 
                         width = 1, height = .5,
                         at = 1985:2010, space = "bottom"), draw = TRUE)
grid.text("Timing of biggest change", x = 0.5, y = 1.12, 
          just = c("centre", "top"), gp = gpar(font = 2, cex = .85))

# key right
upViewport()
vp3 <- viewport(x = 2/3 - .02, y = .72, width = 1/3, height = .1, 
                just = c("left", "bottom"), name = "vp_key3")
pushViewport(vp3)
draw.colorkey(key = list(col = col_maxmagn(1000), width = 1, height = .5,
                         at = c(-.3, brks, .3), space = "bottom"), draw = TRUE)
grid.text("Magnitude of biggest change", x = 0.5, y = 1.12, 
          just = c("centre", "top"), gp = gpar(font = 2, cex = .85))

dev.off()

## deregister parallel backend
stopCluster(cl)
