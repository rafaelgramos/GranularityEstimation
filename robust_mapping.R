require(spatstat)
require(raster)

robust.quadcount<-function(point_set,
                           random_samples=T,
                           nsamples=500,
                           signif=0.99,
                           tradeoff_crit=c("sum","product"),#,"derivative"), #NOT doing derivative right now
                           uniformity_method=c("Quadratcount","Nearest-neighbor"),
                           robustness_method=c("Poisson","Binomial","Resampling"),
                           robustness_k = -3,
                           verbose=F){
  
  W = c(min(point_set[,1]),max(point_set[,1]),min(point_set[,2]),max(point_set[,2]))
  
  #HOW is this initialized?
  my_scales = seq(25,1650,25)
  
  # calculating robustness and uniformity for each granularity in 'my_scales'
  my_spatialstats = get_spatialstats_all(point_set,
                                        my_scales,
                                        my_scales,
                                        random_samples=random_samples,
                                        signif=signif,
                                        uniformity_method=uniformity_method,
                                        robustness_method=robustness_method,
                                        robustness_k = robustness_k,
                                        verbose=verbose)
  
  robust = my_spatialstats$robustness
  unif = my_spatialstats$uniformity
  
  # tradeoff analysis
  
  if(tradeoff_crit == "sum") {
    
    splinefit = smooth.spline(x=my_scales,y=unif+robust,df=10)
    my_spline = predict(splinefit,x=seq(25,1000,5))
    opt_i = which.max(my_spline$y)
    opt_granularity = my_spline$x[opt_i]
    opt_ur = my_spline$y[opt_i]
    
  } else if(tradeoff_crit == "product") {
    
    splinefit = smooth.spline(x=my_scales,y=unif*robust,df=10)
    my_spline = predict(splinefit,x=seq(25,1000,5))
    opt_i = which.max(my_spline$y)
    opt_granularity = my_spline$x[opt_i]
    opt_ur = my_spline$y[opt_i]
    
  #} else if(tradeoff_crit == "derivative") {
    
  #  splinefit = smooth.spline(x=robust,y=unif,df=6)
  #  my_spline = predict(splinefit,x=seq(0,1,0.01))
  #  my_derivs = predict(splinefit,x=seq(0,1,0.01),deriv=1)
  #  opt_i = which.min((my_derivs$y+1)^2)
  #  opt_robust = my_spline$x[opt_i]
  #  opt_uniformity = my_spline$y[opt_i]
    
  #  splinefit = smooth.spline(x=my_scales,y=robust,df=10)
  #  my_spline = predict(splinefit,x=seq(25,1000,5))
  #  opt_i = which.min((my_spline$y-opt_robust)^2)
  #  opt_granularity_1 = my_spline$x[opt_i]
  #  splinefit = smooth.spline(x=my_scales,y=unif,df=10)
  #  my_spline = predict(splinefit,x=seq(25,1000,5))
  #  opt_i = which.min((my_spline$y-opt_uniformity)^2)
  #  opt_granularity_2 = my_spline$x[opt_i]
  #  opt_granularity = 0.5*(opt_granularity_1 + opt_granularity_2)
    
  } else {
    print("Tradeoff criteria not recognized. Options allowed: sum, product.")#, derivative.")
    return(NULL)
  }
  
  map = list()
  nlon = floor((W[2] - W[1])/opt_granularity)
  nlat = floor((W[4] - W[3])/opt_granularity)
  q = quadratcount(as.ppp(point_set,W),nlon,nlat)
  map$counts = raster(matrix(data=q[],nrow=nrow(q),ncol=ncol(q)))
  map$counts = raster::setExtent(map$counts,ext=raster::extent(W))
  map$opt_granularity = opt_granularity

  final_sample = create_samples(point_set=point_set,
                              random_samples = F,
                              window_w = opt_granularity,
                              window_h = opt_granularity)
  
  stats_final_sample = get_spatialstats_sample(sample_set = final_sample,
                                            uniformity_method = uniformity_method,
                                            signif=signif,
                                            robustness_method = robustness_method)
  
  tmp = matrix(data=stats_final_sample$samples_covar,nrow=nrow(map$counts),ncol=ncol(map$counts))
  map$covar = raster::flip(raster(tmp),direction=2)
  map$covar = setExtent(map$covar,extent(map$counts))
  
  tmp = matrix(data=stats_final_sample$samples_csr_pass,nrow=nrow(map$counts),ncol=ncol(map$counts))
  map$is.csr = !(raster::flip(raster(tmp),direction=2))
  map$is.csr = setExtent(map$is.csr,extent(map$counts))
  
  robust = my_spatialstats$robustness
  unif = my_spatialstats$uniformity
  
  map$uniformity.curve = unif
  map$granularities.tested = my_scales
  #map$unif.score = 
  map$robustness.curve =robust
  #map$unit.size = 
  return(map)
}


#
#  Function that estimates the internal uniformity and robustness to error for 'point_set',
#  according to the granularity dimensions listed in 'scales_lon' and 'scales_lat' (these should be the
#  same size; if one is larger than the other, the extra elements in the largest are discarded)
#
#  - random_sample' determinines whether the metric should be estimated from random samples
#    taken from point_set, or from a contiguous regular partion of point_set.
# 
#  - 'n_samples' determines how many samples are taken, if 'random_sample' == T.
# 
#  - 'signif' determines the threshold for considering a sample to be spatially random enough
#    for the sake of calculating internal uniformity.
#

get_spatialstats_all<-function(point_set,
                               scales_lon,
                               scales_lat,
                               random_samples=T,
                               nsamples=500,
                               signif=0.90,
                               robustness_method=c("Poisson","Binomial","Resampling"),
                               robustness_k = -3,
                               uniformity_method=c("Quadratcount","Nearest-neighbor"),
                               verbose=verbose){
  
  if(robustness_k >= 0) {
    print("robust_k parameter must be less than zero. See documentation for details.")
    return(NULL)
  }
  
  len_scales = min(c(length(scales_lon),length(scales_lat)))
  total_points = nrow(point_set)
  
  # allocating vectors for our outputs
  
  # avg number of points per sample (and var)
  scales_count = vector(length=len_scales)
  scales_var = vector(length=len_scales)
  # mean and median for the nearest neighbor p-value test
  scales_csr_p = vector(length=len_scales)
  scales_csr_median_p = vector(length=len_scales)
  # proportion of samples that passed the threshold for CSR using the nearest neighbor approach
  scales_csr_pass = vector(length=len_scales)
  # mean values for the coef_of var (later used to calculate robustness)
  scales_covar = vector(length=len_scales)
  # number of zeros
  scales_zeros = vector(length=len_scales)
  
  for(i in 1:len_scales) {
    if(verbose) print(paste(100*(i/len_scales),"%"))
    
    # generating set of quadrats of the current granularity for taking the samples
    my_samples = create_samples(point_set=point_set,
                                random_samples = random_samples,
                                nsamples = nsamples,
                                window_w = scales_lon[i],
                                window_h = scales_lat[i])
    
    my_stats_sample = get_spatialstats_sample(sample_set = my_samples,
                                              signif = signif,
                                              uniformity_method = uniformity_method,
                                              robustness_method = robustness_method)
    
    # avg number of points per sample (and var)
    scales_count[i] = mean(my_stats_sample$samples_count[!is.na(my_stats_sample$samples_count)])
    scales_var[i] = var(my_stats_sample$samples_count[!is.na(my_stats_sample$samples_count)])
    
    # mean and median of the p-value for the CSR tests
    scales_csr_p[i] = mean(my_stats_sample$samples_csr_p[!is.na(my_stats_sample$samples_csr_p)])
    scales_csr_median_p[i] = median(my_stats_sample$samples_csr_p[!is.na(my_stats_sample$samples_csr_p)])
    # proportion of samples that passed the threshold for CSR
    scales_csr_pass[i] = mean(my_stats_sample$samples_csr_pass[!is.na(my_stats_sample$samples_csr_pass)])
    
    # mean values for the coef_of var
    scales_covar[i] = mean(my_stats_sample$samples_covar[!is.na(my_stats_sample$samples_covar)])
    # number of zeros
    scales_zeros[i] = my_stats_sample$samples_zeros/my_samples$nsample
  }
  
  scales_robustness = exp(-3*scales_covar)
  scales_uniformity = 1-scales_csr_pass
  
  out = list(count = scales_count,
                  var = scales_var,
                  #csr_p = scales_csr_p,
                  #csr_median_p = scales_csr_median_p,
                  uniformity = scales_uniformity,
                  robustness = scales_robustness,
                  zerors = scales_zeros
  )
  return(out)
}

create_samples<-function(point_set,random_samples,nsamples,window_w,window_h) {
  max_b_lon = max(point_set$lon)
  max_b_lat = max(point_set$lat)
  min_b_lon = min(point_set$lon)
  min_b_lat = min(point_set$lat)
  w_b = max_b_lon-min_b_lon
  h_b = max_b_lat-min_b_lat
  
  if(random_samples == T) {
    # quadrats taken randomly
    offset_lon = ((w_b-window_w)*runif(nsamples))+min_b_lon
    offset_lat = ((h_b-window_h)*runif(nsamples))+min_b_lat
  }
  else {
    # contiguous quadrats
    nlon = floor(w_b/window_w)
    nlat = floor(h_b/window_h)
    nsamples = nlon*nlat
    offset_lon = matrix(nrow=nlat,ncol=nlon)
    offset_lat = matrix(nrow=nlat,ncol=nlon)
    for(k in 1:nlat) offset_lon[k,] = (0:(nlon-1))*window_w+min_b_lon
    for(k in 1:nlon) offset_lat[,k] = (0:(nlat-1))*window_h+min_b_lat
    offset_lon = as.vector(offset_lon)
    offset_lat = as.vector(offset_lat)
  }
  out = list()
  out$point_set = point_set
  out$offset_lon = offset_lon
  out$offset_lat = offset_lat
  out$window_w = window_w
  out$window_h = window_h
  out$random_samples = random_samples
  out$npopulation = nrow(point_set)
  out$nsamples = nsamples
  
  return(out)
}

get_spatialstats_sample<-function(sample_set,
                                  signif,
                                  robustness_method=c("Poisson","Binomial","Resampling"),
                                  uniformity_method=c("Quadratcount","Nearest-neighbor")){
  
  nsamples = length(sample_set$offset_lon)
  
  # allocating temporary structures
  samples_count = vector(length=nsamples)
  samples_csr_p = vector(length=nsamples)
  samples_csr_pass = vector(length=nsamples)
  samples_covar = vector(length=nsamples)
  samples_zeros = 0
  
  # for each quadrat, take a the points inside and test for robustness and internal uniformity
  for(j in 1:nsamples) {
    
    # extracting the points inside the quadrat
    W_ext = c(sample_set$offset_lon[j],
              sample_set$offset_lon[j]+sample_set$window_w,
              sample_set$offset_lat[j],
              sample_set$offset_lat[j]+sample_set$window_h)
    
    sub_points = sample_set$point_set[inside.owin(x=sample_set$point_set$lon,
                                                  y=sample_set$point_set$lat,
                                                  w=W_ext),]
    
    # estimate robustness and uniformity
    if(nrow(sub_points)==0) {
      # if there are no points inside, the metrics are NA
      samples_zeros = samples_zeros + 1
      samples_count[j] = NA
      samples_csr_p[j] = NA
      samples_csr_pass[j] = NA
      samples_covar[j] = NA
    }
    else {
      # in case we have points, calculate the metric
      
      # allocating a slightly larger quadrat just to avoid that our sampled points to fall
      # exactly in the boundary of the quadrat.
      W_disp = c(sample_set$offset_lon[j]-0.0001,
                 sample_set$offset_lon[j]+sample_set$window_w+0.0001,
                 sample_set$offset_lat[j]-0.0001,
                 sample_set$offset_lat[j]+sample_set$window_h+0.0001)
      
      # testing Complete Spatial Randomness (CSR) with Clark-Evans nearest neighbor test (part of estimating uniformity)
      samples_count[j] = nrow(sub_points)
      
      if(uniformity_method == "Nearest-neighbor") {
        if(samples_count[j] > 20) {
          my_clarkevans = clarkevans.test(as.ppp(sub_points,W_disp),alternative="two.sided")#,correction = "Donnelly")
        } else {
          my_clarkevans = clarkevans.test(as.ppp(sub_points,W_disp),alternative="two.sided",nsim=100)#,correction = "Donnelly")
        }
        samples_csr_p[j] =  my_clarkevans$p.value
        # assign T to quadrats that pass the threshold 'signif' for CSR
        samples_csr_pass[j] =  my_clarkevans$p.value < (1-signif)
      }
      else if(uniformity_method == "Quadratcount") {
        # test CST with the quadrat test
        my_quadrattest = quadrat.test(as.ppp(sub_points,W_disp),nx=5)
        samples_csr_p[j] = my_quadrattest$p.value
        samples_csr_pass[j] = my_quadrattest$p.value < (1-signif)
      }
      else {
        samples_csr_p[j] = NA
        samples_csr_pass[j] = NA
      }
      
      samples_covar[j] = calc_covar(nrow(sub_points),sample_set$npopulation,robustness_method)
      # Calculating robustness.
      # First we calculate the expected coefficient of variation for the samples, using different methods
      #prob_event = nrow(sub_points)/sample_set$npopulation
      #if(robustness == "Binomial") {
      #  # calculating coef of var using the Binomial estimation method (see paper)
      #  samples_covar[j] = sqrt(prob_event*(1-prob_event)*sample_set$npopulation)/nrow(sub_points)#sd(sim_rates)/mean(sim_rates)
      #}
      #else if(robustness == "Poisson") {
      #  # calculating coef of var using the Poisson estimation method (see paper)
      #  tmp = fitdistr(nrow(sub_points),"Poisson")
      #  samples_covar[j] = tmp$sd/tmp$estimate
      #}
      #else if(robustness == "Resampling") {
      #  # calculating coef of var using the resampling estimation method (see paper)
      #  sim_rates = rbinom(1000,sample_set$npopulation,prob_event)
      #  samples_covar[j] = mean(sqrt(1/sim_rates[sim_rates!=0]))
      #}
      
    }
  }
  out = list()
  out$samples_count = samples_count
  out$samples_csr_p = samples_csr_p
  out$samples_csr_pass = samples_csr_pass
  out$samples_covar = samples_covar
  out$samples_zeros = samples_zeros
  return(out)
}


calc_covar<-function(nsub,ntotal,robustness_method){
  # Calculating robustness.
  # First we calculate the expected coefficient of variation for the samples, using different methods
  prob_event = nsub/ntotal
  if(robustness_method == "Binomial") {
    # calculating coef of var using the Binomial estimation method (see paper)
    samples_covar = sqrt(prob_event*(1-prob_event)*ntotal)/nsub#sd(sim_rates)/mean(sim_rates)
  }
  else if(robustness_method == "Poisson") {
    # calculating coef of var using the Poisson estimation method (see paper)
    tmp = fitdistr(nsub,"Poisson")
    samples_covar = tmp$sd/tmp$estimate
  }
  else if(robustness_method == "Resampling") {
    # calculating coef of var using the resampling estimation method (see paper)
    sim_rates = rbinom(1000,ntotal,prob_event)
    samples_covar = mean(sqrt(1/sim_rates[sim_rates!=0]))
  }
  return(samples_covar)
}





########################3




#
#   Utility functions, for ploting quadrat counts etc.
#

quad2mat<-function(quadset) {
  tmp = matrix(data=quadset[],nrow=nrow(quadset),ncol=ncol(quadset))
  #tmp2 = apply(tmp, 2, rev)
  #tmp3 = t(tmp2)
  return(tmp)
}

points2quad<-function(point_set,my_scale,mask=NULL){
  max_b_lon = max(point_set$lon)
  max_b_lat = max(point_set$lat)
  min_b_lon = min(point_set$lon)
  min_b_lat = min(point_set$lat)
  W = c(min(point_set$lon),
        max(point_set$lon),
        min(point_set$lat),
        max(point_set$lat))
  nlon = floor((max_b_lon - min_b_lon)/my_scale)
  nlat = floor((max_b_lat - min_b_lat)/my_scale)
  my_ret = quadratcount(as.ppp(point_set,W),nlon,nlat)
  
  if(!is.null(mask)){
    Wmask = c(min(mask$lon),
          max(mask$lon),
          min(mask$lat),
          max(mask$lat))
    my_mask = quadratcount(as.ppp(mask,Wmask),nlon,nlat) 
    my_ret[my_mask <= 0] = NA
  }
  
  return(my_ret)
}