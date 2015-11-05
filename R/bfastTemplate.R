bfastTemplate <- function() {
  ## gimms raw data
  gimms_files <- gimms::rearrangeFiles(dsn = "data", full.names = TRUE)
  gimms_raster <- raster::stack(gimms_files[1])
  
  ## reference extent
  library(rworldmap)
  data("countriesCoarse")
  spy_iran <- subset(countriesCoarse, ADMIN == "Iran")
  
  ## crop global gimms images
  rst_template <- crop(gimms_rasters[[1]], spy_iran)
  rst_template[] <- NA
  return(rst_template)
}