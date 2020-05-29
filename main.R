setwd("/Users/rafaelramos/crimeanalysis/projects/robust_mapping")

require(rgdal)
#library(dplyr)
#library(spatstat)
#library(MASS)
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

my_map = robust.quadcount(burglary,tradeoff_crit = "product",uniformity="Nearest-neighbor",robustness="Poisson")

plot(my_map$count)
