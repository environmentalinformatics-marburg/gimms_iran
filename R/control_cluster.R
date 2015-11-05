fls_ndvi <- rearrangeFiles(dsn = "data/dsn", pattern = "^DSN.*VI3g.tif$", 
                           full.names = TRUE, pos = c(4, 6, 11) + 8)

fls_ndvi <- fls_ndvi[grep("82jan15a", fls_ndvi):
                       grep("13dec15b", fls_ndvi)]
rst_ndvi <- stack(fls_ndvi)

## denoising
rst_ndvi_dns <- remote::denoise(rst_ndvi, expl.var = .9)

## reference extent
dat_ndvi_dns <- extract(rst_ndvi_dns, spy_iran, df = TRUE, cellnumbers = TRUE)
dat_ndvi_dns <- dat_ndvi_dns[complete.cases(dat_ndvi_dns), ]

## clustering
dat_ndvi_dns <- as.data.frame(rst_ndvi_dns)
kmn_ndvi_dns <- kmeans(dat_ndvi_dns[, -c(1, 2)], 7, iter.max = 100, 
                       nstart = 5, algorithm = "Lloyd")

rst_ndvi_kmn <- rst_ndvi[[1]]
rst_ndvi_kmn[] <- NA
rst_ndvi_kmn[dat_ndvi_dns$cell] <- kmn_ndvi_dns$cluster

## eot
rst_ndvi_dns <- mask(rst_ndvi_dns, spy_iran)
rst_ndvi_eot <- eot(rst_ndvi_dns, n = 7)
  