setwd("/Users/rafaelramos/crimeanalysis/projects/paper_code")

require(rgdal)
library(dplyr)
library(spatstat)
library(MASS)
library(raster)
source("functions.R")

#
# ======================================================================================================
#

# ===========================================
# Load data
# ===========================================

shape = readOGR(dsn = ".", layer = "burglary")
burglary = data.frame(lon=shape$lon_m,lat=shape$lat_m,year=as.numeric(shape$ANOFATO)+2007,target=shape$V19,month=as.numeric(shape$V5),weekday=as.numeric(shape$DIADASEMAN),time=shape$V10,furtoOuRoubo=shape$V12)
burglary = burglary %>% mutate(hour=as.numeric(substr(as.character(time),1,2)))

shape = readOGR(dsn = ".", layer = "homicidio")
homic = data.frame(lon=shape$X,lat=shape$Y)

shape = readOGR(dsn = ".", layer = "roubo_a_transeunte")
robbery = data.frame(lon=shape$LON_M,lat=shape$LAT_M,year=as.numeric(as.character(shape$ANO_FATO)),hour=as.numeric(as.character(shape$HORA)),weekday=as.numeric(as.character(shape$DIA_SEMANA)))

shape = readOGR(dsn = ".", layer = "landuse")
landuse = data.frame(lon=shape$lon,lat=shape$lat,tipo=shape$CONSTRUTIV,acab=shape$ACABAMENTO)

#
# ======================================================================================================
#

# ===========================================
# Estimating Optimal Granularity for Burglary
# ===========================================

# setting the granularities to be tested
burg_scales = seq(25,1650,25) # in meters
burg_scales_lon = burg_scales
burg_scales_lat = burg_scales

# estimating uniformity and robustness using the random samples method
# SHOULD TAKE A FEW MINUTES, BUT NOT TOO LONG
burg_spatstats_random = get_spatialstats_all(burglary,burg_scales_lon,burg_scales_lat,random_samples=T,signif=0.99)

# estimating uniformity and robustness using the contiguous samples method
# SHOULD TAKE A LONG TIME - RESULT IS SIMILAR (THOUGH NOT IDENTICAL) TO THE FORMER METHOD, 
# SO IN MOST CASES THE PRIOR WILL SUFFICE
burg_spatstats_contig = get_spatialstats_all(burglary,burg_scales_lon,burg_scales_lat,random_samples=F,signif=0.99)


#
# Burglary - Uniformity plots
#

tiff("burg_unif_robust.tiff",width=1800,height=800,res=200)
#tiff("simple.tiff",width=1800,height=800,res=200)
par(mfrow=c(1,2))

plot(NULL,NULL,xlab="Quadrat size (meters)",ylab="Internal Uniformity",xlim=c(25,1000),ylim=c(0,1))
points(burg_scales,1-burg_spatstats_random$nn_pass,col="red",pch=16)
points(burg_scales,1-burg_spatstats_random$quadrat_pass,col="black",pch=16)
points(burg_scales,1-burg_spatstats_contig$nn_pass,col="blue",pch=17)
points(burg_scales,1-burg_spatstats_contig$quadrat_pass,col="green",pch=17)
title("Internal uniformity - Burglary")
#legend("bottomleft",
#       legend = c("Internal Uniformity estimate"),
#       col = c("red"), 
#       pch = c(16), 
#       bty = "n", 
#       pt.cex = 1.5, 
#       cex = 1, 
#       text.col = "black", 
#       horiz = F , 
#       inset = c(0.05, 0.05),
#       y.intersp=2)
legend("bottomleft", 
       legend = c("Nearest-neighbor - random", "Quadrat Count - random",
                  "Nearest-neighbor - contiguous","Quadrat Count - contiguous"), 
       col = c("red", 
               "black",
               "blue",
               "green"), 
       pch = c(16,16,17,17), 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 0.5, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.05, 0.05),
       y.intersp=2)


#
# Burglary - Robustness to error plots
#

plot(NULL,NULL,xlab="Quadrat size (meters)",ylab="Robustness to error",xlim=c(25,1000),ylim=c(0,0.8))
# note that 1-burg_spatstats_random$error_binom is due to robusness = 1 - coef_var
points(burg_scales,exp(-3*burg_spatstats_random$error_binom),col="black",pch=16)
points(burg_scales,exp(-3*burg_spatstats_random$error_pois),col="black",pch=4)
points(burg_scales,exp(-3*burg_spatstats_random$error_boot),col="red",pch=16)
points(burg_scales,exp(-3*burg_spatstats_contig$error_binom),col="green",pch=16)
points(burg_scales,exp(-3*burg_spatstats_contig$error_pois),col="black",pch=1)
points(burg_scales,exp(-3*burg_spatstats_contig$error_boot),col="blue",pch=16)
title("Robustness to error - Burglary")

#legend("bottomright",
#       legend = c("Robustness to error estimate"),
#       col = c("black"), 
#       pch = c(16), 
#       bty = "n", 
#       pt.cex = 1.5, 
#       cex = 1, 
#       text.col = "black", 
#       horiz = F , 
#       inset = c(0.05, 0.01),
#       y.intersp=2)

legend("bottomright", 
       legend = c("Binomial - random", "Poisson - random", "Resampled - random",
                  "Binomial - contiguous", "Poisson - contiguous", "Resampled - contiguous"), 
       col = c("black", 
               "black",
               "red",
               "green",
               "black",
               "blue"), 
       pch = c(16,4,16,16,1,16), 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 0.5, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.15, 0.01),
       y.intersp=2)

dev.off()

#
# Burglary - Tradeoff analysis
#

tiff("burg_tradeoff.tiff",width=1800,height=800,res=200)
par(mfrow=c(1,3),oma=c(0,0,3,0))
robust = exp(-3*burg_spatstats_random$error_binom)#burg_spatstats_random$new_robust25#1-burg_spatstats_random$error_boot
unif = 1-burg_spatstats_random$nn_pass
my_scales = burg_scales

# Criteria 1: du/dr = -1
plot(robust,unif,xlab="Robustness to error",ylab="Internal Uniformity",pch=16,xlim=c(0,1))
splinefit = smooth.spline(x=robust,y=unif,df=6)
my_spline = predict(splinefit,x=seq(0,1,0.01))
my_derivs = predict(splinefit,x=seq(0,1,0.01),deriv=1)
opt_i = 32#which.min((my_derivs$y+1)^2)
opt_robust = my_spline$x[opt_i]
opt_uniformity = my_spline$y[opt_i]
lines(my_spline$x,my_spline$y,col="red")
points(opt_robust,opt_uniformity,col="blue",pch=16,cex=2)
title("Balance of gains criterion")
legend("bottomleft", 
       legend = c("Estimates","Optimal"), 
       col = c("black", 
               "blue"), 
       pch = c(20,16), 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.05, 0.08),
       y.intersp=1.5)
legend("bottomleft", 
       legend = c("Fitted spline"), 
       col = c("red"), 
       pch = "-", 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.05, 0.01),
       y.intersp=1.5)

# Criteria 2: unif*robust
plot(my_scales,unif*robust,ylab="Robustness * Uniformity",xlab="Granularity (meters)",xlim=c(25,1000),pch=16)
splinefit = smooth.spline(x=my_scales,y=unif*robust,df=10)
my_spline = predict(splinefit,x=seq(25,1000,5))
lines(my_spline$x,my_spline$y,col="red")
opt_i = which.max(my_spline$y)
opt_granularity = my_spline$x[opt_i]
opt_ur = my_spline$y[opt_i]
points(opt_granularity,opt_ur,col="blue",pch=16,cex=2)
lines(c(-100,5100),c(opt_uniformity*opt_robust,opt_uniformity*opt_robust),col="darkgreen")
title("Product criterion")
legend("bottomleft", 
       legend = c("Estimates","Optimal"), 
       col = c("black", 
               "blue"), 
       pch = list(20,16), 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.05, 0.17),
       y.intersp=1.5)
legend("bottomleft", 
       legend = c("Fitted spline","Balance of gains"), 
       col = c("red","darkgreen"), 
       pch = "-", 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.05, 0.01),
       y.intersp=1.5)

# Criteria 3: unif + robust
plot(my_scales,unif+robust,ylab="Robustness + Uniformity",xlab="Granularity (meters)",xlim=c(25,1000),pch=16)
splinefit = smooth.spline(x=my_scales,y=unif+robust,df=10)
my_spline = predict(splinefit,x=seq(25,1000,5))
lines(my_spline$x,my_spline$y,col="red")
opt_i = which.max(my_spline$y)
opt_granularity2 = my_spline$x[opt_i]
opt_ur = my_spline$y[opt_i]
points(opt_granularity2,opt_ur,col="blue",pch=16,cex=2)
lines(c(-100,5100),c(opt_uniformity+opt_robust,opt_uniformity+opt_robust),col="darkgreen")
title("Sum criterion")
title("Tradeoff Analysis: Burglary", outer=TRUE, cex.main=1.75)
legend("bottomleft", 
       legend = c("Estimates","Optimal"), 
       col = c("black", 
               "blue"), 
       pch = list(20,16), 
       bty = "n", 
       pt.cex = 1.5, 
       #cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.05, 0.17),
       y.intersp=1.5)
legend("bottomleft", 
       legend = c("Fitted spline","Balance of gains"), 
       col = c("red","darkgreen"), 
       pch = "-", 
       bty = "n", 
       pt.cex = 1.5, 
       #cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.05, 0.01),
       y.intersp=1.5)
dev.off()

# Calculating the mean and sd of the optimal granularity estimated using the different criterias

opt_granularity # Product
opt_granularity2 # Sum
splinefit = smooth.spline(x=my_scales,y=robust,df=10)
my_spline = predict(splinefit,x=seq(25,1000,5))
opt_i = which.min((my_spline$y-opt_robust)^2)
opt_granularity3a = my_spline$x[opt_i]
splinefit = smooth.spline(x=my_scales,y=unif,df=10)
my_spline = predict(splinefit,x=seq(25,1000,5))
opt_i = which.min((my_spline$y-opt_uniformity)^2)
opt_granularity3b = my_spline$x[opt_i]
opt_granularity3a
opt_granularity3b
opt_granularity
opt_granularity2
opt_gran_burg = mean(c(opt_granularity3a,
       opt_granularity3b,
       opt_granularity,
       opt_granularity2))
sd(c(opt_granularity3a,
     opt_granularity3b,
     opt_granularity,
     opt_granularity2))

#
# ======================================================================================================
#



# ===========================================
# Estimating Optimal Granularity for Homicides (essentially the same done with burglaries)
# ===========================================


homic_scales = seq(25,5000,100) # meters
homic_scales_lon = homic_scales
homic_scales_lat = homic_scales

# SHOULD TAKE A FEW MINUTES, BUT NOT TOO LONG
homic_spatstats_random = get_spatialstats_all(homic,homic_scales_lon,homic_scales_lat,random_samples=T,signif=0.99)

# SHOULD TAKE A FEW MINUTES, BUT NOT TOO LONG
homic_spatstats_contig = get_spatialstats_all(homic,homic_scales_lon,homic_scales_lat,random_samples=F,signif=0.99)

#
# Homicides - Uniformity plots
#

tiff("homic_unif_robust.tiff",width=1800,height=800,res=200)
par(mfrow=c(1,2))
plot(NULL,NULL,xlab="Quadrat size (meters)",ylab="Internal Uniformity",xlim=c(25,5000),ylim=c(0,1))
points(homic_scales,1-homic_spatstats_random$nn_pass,col="red",pch=16)
points(homic_scales,1-homic_spatstats_random$quadrat_pass,col="black",pch=16)
points(homic_scales,1-homic_spatstats_contig$nn_pass,col="blue",pch=17)
points(homic_scales,1-homic_spatstats_contig$quadrat_pass,col="green",pch=17)
title("Internal uniformity - Homicides")
legend("bottomleft", 
       legend = c("Nearest-neighbor - random", "Quadrat Count - random",
                  "Nearest-neighbor - contiguous","Quadrat Count - contiguous"), 
       col = c("red", 
               "black",
               "blue",
               "green"), 
       pch = c(16,16,17,17), 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 0.5, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.05, 0.05),
       y.intersp=2)

# 
# Homicides - Robustness to error plots
#

plot(NULL,NULL,xlab="Quadrat size (meters)",ylab="Robustness to error",xlim=c(25,5000),ylim=c(0,0.8))
points(homic_scales,exp(-3*homic_spatstats_random$error_binom),col="black",pch=16)
points(homic_scales,exp(-3*homic_spatstats_random$error_pois),col="black",pch=4)
points(homic_scales,exp(-3*homic_spatstats_random$error_boot),col="red",pch=16)
points(homic_scales,exp(-3*homic_spatstats_contig$error_binom),col="green",pch=16)
points(homic_scales,exp(-3*homic_spatstats_contig$error_pois),col="black",pch=1)
points(homic_scales,exp(-3*homic_spatstats_contig$error_boot),col="blue",pch=16)
title("Robustness to error - Homicides")
legend("bottomright", 
       legend = c("Binomial - random", "Poisson - random", "Resampled - random",
                  "Binomial - contiguous", "Poisson - contiguous", "Resampled - contiguous"), 
       col = c("black", 
               "black",
               "red",
               "green",
               "black",
               "blue"), 
       pch = c(16,4,16,16,1,16), 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 0.5, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.05, 0.01),
       y.intersp=2)

dev.off()

# 
# Homicides - Tradeoff Analysis
#

tiff("homic_tradeoff.tiff",width=1800,height=800,res=200)
par(mfrow=c(1,3),oma=c(0,0,3,0))
robust = exp(-3*homic_spatstats_random$error_boot)
unif = 1-homic_spatstats_random$nn_pass
my_scales = homic_scales

# Criteria 1: du/dr = -1
plot(robust,unif,xlab="Robustness to error",ylab="Internal Uniformity",pch=16)
splinefit = smooth.spline(x=robust,y=unif,df=10)
my_spline = predict(splinefit,x=seq(0,1,0.01))
my_derivs = predict(splinefit,x=seq(0,1,0.01),deriv=1)
opt_i = which.min((my_derivs$y+1)^2)
opt_robust = my_spline$x[opt_i]
opt_uniformity = my_spline$y[opt_i]
lines(my_spline$x,my_spline$y,col="red")
points(opt_robust,opt_uniformity,col="blue",pch=16,cex=2)
title("Balance of gains criterion")
legend("bottomleft", 
       legend = c("Estimates","Optimal"), 
       col = c("black", 
               "blue"), 
       pch = c(20,16), 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.05, 0.08),
       y.intersp=1.5)
legend("bottomleft", 
       legend = c("Fitted spline"), 
       col = c("red"), 
       pch = "-", 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.05, 0.01),
       y.intersp=1.5)

# Criteria 2: unif*robust
plot(my_scales,unif*robust,ylab="Robustness * Uniformity",xlab="Granularity (meters)",pch=16)
splinefit = smooth.spline(x=my_scales,y=unif*robust,df=10)
my_spline = predict(splinefit,x=seq(25,5000,100))
lines(my_spline$x,my_spline$y,col="red")
opt_i = which.max(my_spline$y)
opt_granularity = my_spline$x[opt_i]
opt_ur = my_spline$y[opt_i]
points(opt_granularity,opt_ur,col="blue",pch=16,cex=2)
lines(c(-100,5100),c(opt_uniformity*opt_robust,opt_uniformity*opt_robust),col="darkgreen")
title("Product criterion")
legend("bottomleft", 
       legend = c("Estimates","Optimal"), 
       col = c("black", 
               "blue"), 
       pch = list(20,16), 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.15, 0.17),
       y.intersp=1.5)
legend("bottomleft", 
       legend = c("Fitted spline","Balance of gains"), 
       col = c("red","darkgreen"), 
       pch = "-", 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.15, 0.01),
       y.intersp=1.5)

# Criteria 3: unif + robust
plot(my_scales,unif+robust,ylab="Robustness + Uniformity",xlab="Granularity (meters)",pch=16)
splinefit = smooth.spline(x=my_scales,y=unif+robust,df=10)
my_spline = predict(splinefit,x=seq(25,5000,100))
lines(my_spline$x,my_spline$y,col="red")
opt_i = which.max(my_spline$y)
opt_granularity2 = my_spline$x[opt_i]
opt_ur = my_spline$y[opt_i]
points(opt_granularity2,opt_ur,col="blue",pch=16,cex=2)
lines(c(-100,5100),c(opt_uniformity+opt_robust,opt_uniformity+opt_robust),col="darkgreen")
title("Sum criterion")
title("Tradeoff Analysis: Homicide", outer=TRUE, cex.main=1.75)
legend("bottomleft", 
       legend = c("Estimates","Optimal"), 
       col = c("black", 
               "blue"), 
       pch = list(20,16), 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.09, 0.17),
       y.intersp=1.5)
legend("bottomleft", 
       legend = c("Fitted spline","Balance of gains"), 
       col = c("red","darkgreen"), 
       pch = "-", 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.09, 0.01),
       y.intersp=1.5)
dev.off()

opt_granularity # Product
opt_granularity2 # Sum
splinefit = smooth.spline(x=my_scales,y=robust,df=10)
my_spline = predict(splinefit,x=seq(25,5000,100))
opt_i = which.min((my_spline$y-opt_robust)^2)
opt_granularity3a = my_spline$x[opt_i]
splinefit = smooth.spline(x=my_scales,y=unif,df=10)
my_spline = predict(splinefit,x=seq(25,5000,100))
opt_i = which.min((my_spline$y-opt_uniformity)^2)
opt_granularity3b = my_spline$x[opt_i]
opt_granularity3a
opt_granularity3b
opt_granularity
opt_granularity2
opt_gran_homic = mean(c(opt_granularity3a,
       opt_granularity3b,
       opt_granularity,
       opt_granularity2))
sd(c(opt_granularity3a,
     opt_granularity3b,
     opt_granularity,
     opt_granularity2))


#
# ======================================================================================================
#


# ===========================================
# Estimating Optimal Granularity for Robbery (essentially the same done with burglaries)
# ===========================================

robbery_scales = seq(25,5000,100) # meters
robbery_scales_lon = robbery_scales
robbery_scales_lat = robbery_scales

# SHOULD TAKE A FEW MINUTES, BUT NOT TOO LONG
robbery_spatstats_random = get_spatialstats_all(robbery,robbery_scales_lon,robbery_scales_lat,random_samples=T,signif=0.99)

# SHOULD TAKE A FEW MINUTES, BUT NOT TOO LONG
robbery_spatstats_contig = get_spatialstats_all(robbery,robbery_scales_lon,robbery_scales_lat,random_samples=F,signif=0.99)

# 
# Robbery - Uniformity plots
#

tiff("robb_unif_robust.tiff",width=1800,height=800,res=200)
par(mfrow=c(1,2))

plot(NULL,NULL,xlab="Quadrat size (meters)",ylab="Internal Uniformity",xlim=c(100,5000),ylim=c(0,1))
points(robbery_scales,1-robbery_spatstats_random$nn_pass,col="red",pch=16)
points(robbery_scales,1-robbery_spatstats_random$quadrat_pass,col="black",pch=16)
points(robbery_scales,1-robbery_spatstats_contig$nn_pass,col="blue",pch=17)
points(robbery_scales,1-robbery_spatstats_contig$quadrat_pass,col="green",pch=17)
title("Internal uniformity - Robbery")
legend("topright", 
       legend = c("Nearest-neighbor - random", "Quadrat Count - random",
                  "Nearest-neighbor - contiguous","Quadrat Count - contiguous"), 
       col = c("red", 
               "black",
               "blue",
               "green"), 
       pch = c(16,16,17,17), 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 0.5, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.15, 0.05),
       y.intersp=2)

#
# Robbery - Robustness to error plots
#

plot(NULL,NULL,xlab="Quadrat size (meters)",ylab="Robustness to error",xlim=c(100,5000),ylim=c(0,0.8))
points(robbery_scales,exp(-3*robbery_spatstats_random$error_binom),col="black",pch=16)
points(robbery_scales,exp(-3*robbery_spatstats_random$error_pois),col="black",pch=4)
points(robbery_scales,exp(-3*robbery_spatstats_random$error_boot),col="red",pch=16)
points(robbery_scales,exp(-3*robbery_spatstats_contig$error_binom),col="green",pch=16)
points(robbery_scales,exp(-3*robbery_spatstats_contig$error_pois),col="black",pch=1)
points(robbery_scales,exp(-3*robbery_spatstats_contig$error_boot),col="blue",pch=16)
title("Robustness to error - Robbery")
legend("bottomright", 
       legend = c("Binomial - random", "Poisson - random", "Resampled - random",
                  "Binomial - contiguous", "Poisson - contiguous", "Resampled - contiguous"), 
       col = c("black", 
               "black",
               "red",
               "green",
               "black",
               "blue"), 
       pch = c(16,4,16,16,1,16), 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 0.5, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.15, 0.01),
       y.intersp=2)
dev.off()

#
# Robbery - Tradeoff analysis
#

tiff("robb_tradeoff.tiff",width=1800,height=800,res=200)
par(mfrow=c(1,3),oma=c(0,0,3,0))
robust = exp(-3*robbery_spatstats_random$error_boot)
unif = 1-robbery_spatstats_random$nn_pass
my_scales = robbery_scales

# Criteria 1: du/dr = -1
plot(robust,unif,xlab="Robustness to error",ylab="Internal Uniformity",pch=16)
splinefit = smooth.spline(x=robust,y=unif,df=10)
my_spline = predict(splinefit,x=seq(0,1,0.01))
my_derivs = predict(splinefit,x=seq(0,1,0.01),deriv=1)
opt_i = 32#which.min((my_derivs$y+1)^2)
opt_robust = my_spline$x[opt_i]
opt_uniformity = my_spline$y[opt_i]
lines(my_spline$x,my_spline$y,col="red")
points(opt_robust,opt_uniformity,col="blue",pch=16,cex=2)
title("Balance of gains criterion")
legend("bottomleft", 
       legend = c("Estimates","Optimal"), 
       col = c("black", 
               "blue"), 
       pch = c(20,16), 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.05, 0.08),
       y.intersp=1.5)
legend("bottomleft", 
       legend = c("Fitted spline"), 
       col = c("red"), 
       pch = "-", 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.05, 0.01),
       y.intersp=1.5)

# Criteria 2: unif*robust
plot(my_scales,unif*robust,ylab="Robustness * Uniformity",xlab="Granularity (meters)",pch=16)
splinefit = smooth.spline(x=my_scales,y=unif*robust,df=10)
my_spline = predict(splinefit,x=seq(25,5000,100))
lines(my_spline$x,my_spline$y,col="red")
opt_i = which.max(my_spline$y)
opt_granularity = my_spline$x[opt_i]
opt_ur = my_spline$y[opt_i]
points(opt_granularity,opt_ur,col="blue",pch=16,cex=2)
lines(c(-100,5100),c(opt_uniformity*opt_robust,opt_uniformity*opt_robust),col="darkgreen")
title("Product criterion")
legend("bottomleft", 
       legend = c("Estimates","Optimal"), 
       col = c("black", 
               "blue"), 
       pch = list(20,16), 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.05, 0.17),
       y.intersp=1.5)
legend("bottomleft", 
       legend = c("Fitted spline","Balance of gains"), 
       col = c("red","darkgreen"), 
       pch = "-", 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.05, 0.01),
       y.intersp=1.5)

# Criteria 3: unif + robust
plot(my_scales,unif+robust,ylab="Robustness + Uniformity",xlab="Granularity (meters)",pch=16)
splinefit = smooth.spline(x=my_scales,y=unif+robust,df=10)
my_spline = predict(splinefit,x=seq(25,5000,100))
lines(my_spline$x,my_spline$y,col="red")
opt_i = which.max(my_spline$y)
opt_granularity2 = my_spline$x[opt_i]
opt_ur = my_spline$y[opt_i]
points(opt_granularity2,opt_ur,col="blue",pch=16,cex=2)
lines(c(-100,5100),c(opt_uniformity+opt_robust,opt_uniformity+opt_robust),col="darkgreen")
title("Sum criterion")
title("Tradeoff Analysis: Robbery", outer=TRUE, cex.main=1.75)
legend("bottomleft", 
       legend = c("Estimates","Optimal"), 
       col = c("black", 
               "blue"), 
       pch = list(20,16), 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.05, 0.17),
       y.intersp=1.5)
legend("bottomleft", 
       legend = c("Fitted spline","Balance of gains"), 
       col = c("red","darkgreen"), 
       pch = "-", 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.05, 0.01),
       y.intersp=1.5)
dev.off()

opt_granularity # Product
opt_granularity2 # Sum
splinefit = smooth.spline(x=my_scales,y=robust,df=10)
my_spline = predict(splinefit,x=seq(25,5000,100))
opt_i = which.min((my_spline$y-opt_robust)^2)
opt_granularity3a = my_spline$x[opt_i]
splinefit = smooth.spline(x=my_scales,y=unif)
my_spline = predict(splinefit,x=seq(25,5000,100),df=10)
opt_i = which.min((my_spline$y-opt_uniformity)^2)
opt_granularity3b = my_spline$x[opt_i]
opt_granularity3a
opt_granularity3b
opt_granularity
opt_granularity2
opt_gran_robbery = mean(c(opt_granularity3a,
                          opt_granularity3b,
                          opt_granularity,
                          opt_granularity2))
sd(c(opt_granularity3a,
     opt_granularity3b,
     opt_granularity,
     opt_granularity2))


#
# ======================================================================================================
#


# ===========================================
# Proportions for burglary
# ===========================================

std_unit_lon = round((max(burglary$lon)-min(burglary$lon))/278.75)
std_unit_lat = round((max(burglary$lat)-min(burglary$lat))/278.75)

Wland = c(min(landuse$lon),
          max(landuse$lon),
          min(landuse$lat),
          max(landuse$lat))

qlanduse = quadratcount(as.ppp(landuse,Wland),std_unit_lon,std_unit_lat)

Wburg = c(min(burglary$lon),
          max(burglary$lon),
          min(burglary$lat),
          max(burglary$lat))

qburg = quadratcount(as.ppp(burglary,Wburg),std_unit_lon,std_unit_lat)
qburg[qlanduse <= 0] = NA
mat_burg = quad2mat(qburg)
plot(raster(quad2mat(qburg)))

limits = qlanduse > 0

sort_burg = sort(as.vector(mat_burg[limits]), decreasing = T)
sum(sort_burg[1:768])/sum(sort_burg)
768/length(sort_burg[limits])

random_counts = rpois(length(as.vector(mat_burg[limits])),44560/length(as.vector(mat_burg[limits])))
sort_random = sort(random_counts, decreasing = T)
sum(sort_random[1:1326])/sum(sort_random)
1326/length(as.vector(mat_burg[limits]))



# ===========================================
# Proportions for robbery
# ===========================================

std_unit_lon = round((max(burglary$lon)-min(burglary$lon))/775)
std_unit_lat = round((max(burglary$lat)-min(burglary$lat))/775)

Wland = c(min(landuse$lon),
          max(landuse$lon),
          min(landuse$lat),
          max(landuse$lat))

qlanduse = quadratcount(as.ppp(landuse,Wland),std_unit_lon,std_unit_lat)

Wburg = c(min(burglary$lon),
          max(burglary$lon),
          min(burglary$lat),
          max(burglary$lat))

qrobb = quadratcount(as.ppp(robbery,Wburg),std_unit_lon,std_unit_lat)
qrobb[qlanduse <= 0] = NA
mat_robb = quad2mat(qrobb)
plot(raster(quad2mat(qrobb)))

limits = qlanduse > 0

sort_robb = sort(as.vector(mat_robb[limits]), decreasing = T)
sum(sort_robb[1:59])/sum(sort_robb)
59/length(sort_robb[limits])

random_counts = rpois(length(as.vector(mat_robb[limits])),11626/length(as.vector(mat_robb[limits])))
sort_random = sort(random_counts, decreasing = T)
sum(sort_random[1:221])/sum(sort_random)
221/length(as.vector(mat_robb[limits]))


# ===========================================
# Proportions for homic
# ===========================================

std_unit_lon = round((max(burglary$lon)-min(burglary$lon))/1775)
std_unit_lat = round((max(burglary$lat)-min(burglary$lat))/1775)

Wland = c(min(landuse$lon),
          max(landuse$lon),
          min(landuse$lat),
          max(landuse$lat))

qlanduse = quadratcount(as.ppp(landuse,Wland),std_unit_lon,std_unit_lat)

Wburg = c(min(burglary$lon),
          max(burglary$lon),
          min(burglary$lat),
          max(burglary$lat))

qhomic = quadratcount(as.ppp(homic,Wburg),std_unit_lon,std_unit_lat)
qhomic[qlanduse <= 0] = NA
mat_homic = quad2mat(qhomic)
plot(raster(quad2mat(qhomic)))

limits = qlanduse > 0

sort_homic = sort(as.vector(mat_homic[limits]), decreasing = T)
sum(sort_homic[1:26])/sum(sort_homic)
26/length(sort_homic[limits])

random_counts = rpois(length(as.vector(mat_homic[limits])),1825/length(as.vector(mat_homic[limits])))
sort_random = sort(random_counts, decreasing = T)
sum(sort_random[1:49])/sum(sort_random)
49/length(as.vector(mat_homic[limits]))


#
# ======================================================================================================
#


# ===========================================
# Plot maps
# ===========================================

tiff("robb_tradeoff.tiff",width=1800,height=800,res=200)
par(mfrow=c(2,2),oma=c(0,0,3,0))

dev.off()

colfunc <- colorRampPalette(c("white","black"))
colfunc <- colorRampPalette(c("blue","yellow","red"))
colfunc <- colorRampPalette(c("blue","white","red"))

tiff("maps.tiff",width=1350,height=1650,res=200)
par(mfrow=c(2,2),oma=c(0,0,0,0))
plot(raster(quad2mat(qburg)),col=colfunc(20))
title("Burglary (2008-2014)")

plot(raster(quad2mat(qhomic)),col=colfunc(20))
title("Homicide (2012-2014)")

plot(raster(quad2mat(qrobb)),col=colfunc(20))
title("Robbery (2012-2013)")

plot(raster(log(1+quad2mat(qrobb))),col=colfunc(20))
title("Log (1 + Robbery)")
dev.off()


#
# ====
#


# ===========================================
# Different k's comparison
# ===========================================

tiff("diff_k.tiff",width=1800,height=800,res=200)
par(mfrow=c(1,3),oma=c(0,0,3,0))
robust1 = exp(-2*burg_spatstats_random$error_binom)#burg_spatstats_random$new_robust25#1-burg_spatstats_random$error_boot
robust2 = exp(-3*burg_spatstats_random$error_binom)
robust3 = exp(-4*burg_spatstats_random$error_binom)
unif = 1-burg_spatstats_random$nn_pass
my_scales = burg_scales

# Criteria 1: du/dr = -1
plot(robust1,unif,xlab="Robustness to error",ylab="Internal Uniformity",pch=16,xlim=c(0,1))
points(robust2,unif,col="red",pch=16)
points(robust3,unif,col="blue",pch=16)
title("Balance of gains")
legend("bottomleft", 
       legend = c("k = 2","k = 3", "k = 4","Optimal"), 
       col = c("black",
               "red",
               "blue",
               "green"), 
       pch = c(16,16,16,16), 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.025, 0.025),
       y.intersp=1.5)

#
splinefit1 = smooth.spline(x=robust1,y=unif,df=6)
my_spline1 = predict(splinefit1,x=seq(0,1,0.01))
my_derivs1 = predict(splinefit1,x=seq(0,1,0.01),deriv=1)
opt_i1 = 46#which.min((my_derivs1$y+1)^2)
opt_robust1 = my_spline1$x[opt_i1]
opt_uniformity1 = my_spline1$y[opt_i1]
#lines(my_spline1$x,my_spline1$y,col="black")
points(opt_robust1,opt_uniformity1,col="green",pch=16,cex=2)

splinefit2 = smooth.spline(x=robust2,y=unif,df=6)
my_spline2 = predict(splinefit2,x=seq(0,1,0.01))
my_derivs2 = predict(splinefit2,x=seq(0,1,0.01),deriv=1)
opt_i2 = 31#which.min((my_derivs2$y+1)^2)
opt_robust2 = my_spline2$x[opt_i2]
opt_uniformity2 = my_spline2$y[opt_i2]
#lines(my_spline2$x,my_spline2$y,col="red")
points(opt_robust2,opt_uniformity2,col="green",pch=16,cex=2)

splinefit3 = smooth.spline(x=robust3,y=unif,df=6)
my_spline3 = predict(splinefit3,x=seq(0,1,0.01))
my_derivs3 = predict(splinefit3,x=seq(0,1,0.01),deriv=1)
opt_i3 = 20#which.min((my_derivs3$y+1)^2)
opt_robust3 = my_spline3$x[opt_i3]
opt_uniformity3 = my_spline3$y[opt_i3]
#lines(my_spline3$x,my_spline3$y,col="blue")
points(opt_robust3,opt_uniformity3,col="green",pch=16,cex=2)


# Criteria 2: unif*robust
plot(my_scales,unif*robust1,ylab="Robustness * Uniformity",xlab="Granularity (meters)",xlim=c(25,1000),ylim=c(0,0.8),pch=16)
splinefit1 = smooth.spline(x=my_scales,y=unif*robust1,df=10)
my_spline1 = predict(splinefit1,x=seq(25,1000,5))
#lines(my_spline1$x,my_spline1$y,col="black")
opt_i1 = which.max(my_spline1$y)
opt_granularity1 = my_spline1$x[opt_i1]
opt_ur1 = my_spline1$y[opt_i1]
points(opt_granularity1,opt_ur1,col="green",pch=16,cex=2)

points(my_scales,unif*robust2,pch=16,col="red")
splinefit2 = smooth.spline(x=my_scales,y=unif*robust2,df=10)
my_spline2 = predict(splinefit2,x=seq(25,1000,5))
#lines(my_spline2$x,my_spline2$y,col="red")
opt_i2 = which.max(my_spline2$y)
opt_granularity2 = my_spline2$x[opt_i2]
opt_ur2 = my_spline2$y[opt_i2]
points(opt_granularity2,opt_ur2,col="green",pch=16,cex=2)

points(my_scales,unif*robust3,pch=16,col="blue")
splinefit3 = smooth.spline(x=my_scales,y=unif*robust3,df=10)
my_spline3 = predict(splinefit3,x=seq(25,1000,5))
#lines(my_spline3$x,my_spline3$y,col="blue")
opt_i3 = which.max(my_spline3$y)
opt_granularity3 = my_spline3$x[opt_i3]
opt_ur3 = my_spline3$y[opt_i3]
points(opt_granularity3,opt_ur3,col="green",pch=16,cex=2)

title("Product criterion")
legend("bottomleft", 
       legend = c("k = 2","k = 3","k = 4","Optimal"), 
       col = c("black",
               "red",
               "blue",
               "green"), 
       pch = list(20,16), 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.5, 0.57),
       y.intersp=1.5)

# 

plot(my_scales,unif+robust1,ylab="Robustness + Uniformity",xlab="Granularity (meters)",xlim=c(25,1000),ylim=c(0.4,1.8),pch=16)
splinefit1 = smooth.spline(x=my_scales,y=unif+robust1,df=10)
my_spline1 = predict(splinefit1,x=seq(25,1000,5))
#lines(my_spline1$x,my_spline1$y,col="black")
opt_i1 = which.max(my_spline1$y)
opt_granularity1 = my_spline1$x[opt_i1]
opt_ur1 = my_spline1$y[opt_i1]
points(opt_granularity1,opt_ur1,col="green",pch=16,cex=2)

points(my_scales,unif+robust2,pch=16,col="red")
splinefit2 = smooth.spline(x=my_scales,y=unif+robust2,df=10)
my_spline2 = predict(splinefit2,x=seq(25,1000,5))
#lines(my_spline2$x,my_spline2$y,col="red")
opt_i2 = which.max(my_spline2$y)
opt_granularity2 = my_spline2$x[opt_i2]
opt_ur2 = my_spline2$y[opt_i2]
points(opt_granularity2,opt_ur2,col="green",pch=16,cex=2)

points(my_scales,unif+robust3,pch=16,col="blue")
splinefit3 = smooth.spline(x=my_scales,y=unif+robust3,df=10)
my_spline3 = predict(splinefit3,x=seq(25,1000,5))
#lines(my_spline3$x,my_spline3$y,col="blue")
opt_i3 = which.max(my_spline3$y)
opt_granularity3 = my_spline3$x[opt_i3]
opt_ur3 = my_spline3$y[opt_i3]
points(opt_granularity3,opt_ur3,col="green",pch=16,cex=2)

title("Sum criterion")
legend("bottomleft", 
       legend = c("k = 2","k = 3","k = 4","Optimal"), 
       col = c("black",
               "red",
               "blue",
               "green"), 
       pch = list(20,16), 
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1.0, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.65, 0.57),
       y.intersp=1.5)

dev.off()
