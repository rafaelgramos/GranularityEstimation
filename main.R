setwd("/Users/rafaelramos/crimeanalysis/projects/robust_mapping")

require(rgdal)
#library(dplyr)
#library(spatstat)
library(MASS)
#library(raster)
source("robust_mapping.R")

#
# ======================================================================================================
#

# ===========================================
# Load data
# ===========================================

shape = readOGR(dsn = ".", layer = "burglary")
burglary = data.frame(lon=shape$lon_m,lat=shape$lat_m,furtoOuRoubo=shape$V12,data=shape$burg_date,time=shape$hora)

#
# ======================================================================================================
#

# ===========================================
# Running granularity analysis
# ===========================================

my_map = robust.quadcount(burglary,
                          tradeoff_crit = "product",
                          uniformity="Nearest-neighbor",
                          robustness_method="Poisson",
                          signif=0.99,
                          verbose=T)

mean(my_map$is.csr[!is.na(my_map$is.csr)])
mean(my_map$covar[!is.na(my_map$covar)])
my_map$opt_granularity

plot(my_map$count)
plot(my_map$is.csr)
plot(my_map$covar)
plot(my_map$covar*my_map$count)

plot(my_map$granularities.tested, my_map$uniformity.curve)
points(my_map$granularities.tested, my_map$robustness.curve,col="red")


