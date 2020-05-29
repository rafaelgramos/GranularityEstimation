
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

get_spatialstats_all<-function(point_set,scales_lon,scales_lat,random_samples=T,nsamples=500,signif=0.90){
  len_scales = min(c(length(scales_lon),length(scales_lat)))
  total_points = nrow(point_set)
  
  # set window of interest that covers all points, but not more
  max_b_lon = max(point_set$lon)
  max_b_lat = max(point_set$lat)
  min_b_lon = min(point_set$lon)
  min_b_lat = min(point_set$lat)
  w_b = max_b_lon-min_b_lon
  h_b = max_b_lat-min_b_lat
  
  # allocating vectors for our outputs
  
  # avg number of points per sample (and var)
  scales_count = vector(length=len_scales)
  scales_var = vector(length=len_scales)
  # mean and median for the nearest neighbor p-value test
  scales_nn_p = vector(length=len_scales)
  scales_nn_median_p = vector(length=len_scales)
  # proportion of samples that passed the threshold for CSR using the nearest neighbor approach
  scales_nn_pass = vector(length=len_scales)
  # mean and median for the p-value of the quadrat test
  scales_quadrat_p = vector(length=len_scales)
  # proportion of samples that passed the threshold for CSR using the quadrat test
  scales_quadrat_median_p = vector(length=len_scales)
  scales_quadrat_pass = vector(length=len_scales)
  # mean values for the coef_of var (later used to calculate robustness)
  scales_error_binom = vector(length=len_scales)
  scales_error_pois = vector(length=len_scales)
  scales_error_boot = vector(length=len_scales)
  scales_error_comp = vector(length=len_scales)
  scales_new_robust10 = vector(length=len_scales)
  scales_new_robust25 = vector(length=len_scales)
  scales_new_robust5 = vector(length=len_scales)
  scales_robust_pass = vector(length=len_scales)
  # number of zeros
  scales_zeros = vector(length=len_scales)
  
  for(i in 1:len_scales) {
    print(paste(100*(i/len_scales),"%"))
    window_w = scales_lon[i]
    window_h = scales_lat[i]
    
    # generating set of quadrats of the current granularity for taking the samples
    if(random_samples == T) {
      # quadrats taken randomly
      iter2=nsamples
      offset_lon = ((w_b-window_w)*runif(nsamples))+min_b_lon
      offset_lat = ((h_b-window_h)*runif(nsamples))+min_b_lat
    }
    else {
      # contiguous quadratas
      nlon = floor(w_b/scales_lon[i])
      nlat = floor(h_b/scales_lat[i])
      iter2=nlon*nlat
      offset_lon = matrix(nrow=nlat,ncol=nlon)
      offset_lat = matrix(nrow=nlat,ncol=nlon)
      for(k in 1:nlat) offset_lon[k,] = (0:(nlon-1))*scales_lon[i]+min_b_lon
      for(k in 1:nlon) offset_lat[,k] = (0:(nlat-1))*scales_lat[i]+min_b_lat
      offset_lon = as.vector(offset_lon)
      offset_lat = as.vector(offset_lat)
    }
    
    # allocating temporary structures
    samples_count = vector(length=iter2)
    samples_nn_p = vector(length=iter2)
    samples_nn_pass = vector(length=iter2)
    samples_quadrat_p = vector(length=iter2)
    samples_quadrat_pass = vector(length=iter2)
    samples_error_binom = vector(length=iter2)
    samples_error_pois = vector(length=iter2)
    samples_error_boot = vector(length=iter2)
    samples_error_comp = vector(length=iter2)
    samples_new_robust10 = vector(length=iter2)
    samples_new_robust25 = vector(length=iter2)
    samples_new_robust5 = vector(length=iter2)
    samples_robust_pass = vector(length=iter2)
    samples_var = vector(length=iter2)
    samples_zeros = 0
    
    # for each quadrat, take a the points inside and test for robustness and internal uniformity
    for(j in 1:iter2) {
      # extracting the points inside the quadrat
      W_ext = c(offset_lon[j],
                offset_lon[j]+window_w,
                offset_lat[j],
                offset_lat[j]+window_h)
      sub_points = point_set[inside.owin(x=point_set$lon,y=point_set$lat,w=W_ext),]
      
      # estimate robustness and uniformity
      nrow(sub_points)
      if(nrow(sub_points)==0) {
        # if there are no points inside, the metrics are NA
        samples_zeros = samples_zeros + 1
        samples_count[j] = NA
        samples_nn_p[j] = NA
        samples_nn_pass[j] = NA
        samples_quadrat_p[j] = NA
        samples_quadrat_pass[j] = NA
        samples_error_pois[j] = NA
        samples_error_binom[j] = NA
        samples_error_boot[j] = NA
        samples_error_comp[j] = NA
        samples_new_robust10[j] = NA
        samples_new_robust25[j] = NA
        samples_new_robust5[j] = NA
        samples_robust_pass[j] = NA
        samples_var[j] = NA
      }
      else {
        # in case we have points, calculate the metric
  
        # allocating a slightly larger quadrat just to avoid that our sampled points to fall
        # exactly in the boundary of the quadrat.
        W_disp = c(offset_lon[j]-0.0001,
                   offset_lon[j]+window_w+0.0001,
                   offset_lat[j]-0.0001,
                   offset_lat[j]+window_h+0.0001)
        
        # testing Complete Spatial Randomness (CSR) with Clark-Evans nearest neighbor test (part of estimating uniformity)
        samples_count[j] = nrow(sub_points)
        if(samples_count[j] > 20) {
          my_clarkevans = clarkevans.test(as.ppp(sub_points,W_disp),alternative="two.sided")#,correction = "Donnelly")
        } else {
          my_clarkevans = clarkevans.test(as.ppp(sub_points,W_disp),alternative="two.sided",nsim=100)#,correction = "Donnelly")
        }
        samples_nn_p[j] =  my_clarkevans$p.value
        # assign T to quadrats that pass the threshold 'signif' for CSR
        samples_nn_pass[j] =  my_clarkevans$p.value < (1-signif)
        
        # test CST with the quadrat test
        my_quadrattest = quadrat.test(as.ppp(sub_points,W_disp),nx=5)
        samples_quadrat_p[j] = my_quadrattest$p.value
        samples_quadrat_pass[j] = my_quadrattest$p.value < (1-signif)
        
        # Calculating robustness.
        # First we calculate the expected coefficient of variation for the samples, using different methods
        prob_event = nrow(sub_points)/total_points
        # calculating coef of var using the Binomial estimation method (see paper)
        samples_error_binom[j] = sqrt(prob_event*(1-prob_event)*total_points)/nrow(sub_points)#sd(sim_rates)/mean(sim_rates)
        
        # calculating coef of var using the Poisson estimation method (see paper)
        tmp = fitdistr(nrow(sub_points),"Poisson")
        samples_error_pois[j] = tmp$sd/tmp$estimate
        
        # calculating coef of var using the resampling estimation method (see paper)
        sim_rates = rbinom(1000,total_points,prob_event)
        samples_error_boot[j] = mean(sqrt(1/sim_rates[sim_rates!=0]))
        
        my_test = poisson.test(round((1+0.10)*tmp$estimate),round(tmp$estimate))
        samples_new_robust10[j] = 1-my_test$p.value
        
        my_test = poisson.test(round((1+0.25)*tmp$estimate),round(tmp$estimate))
        samples_new_robust25[j] = 1-my_test$p.value
        
        my_test = poisson.test(round((1+0.05)*tmp$estimate),round(tmp$estimate))
        samples_new_robust5[j] = 1-my_test$p.value
        
        samples_robust_pass[j] = samples_new_robust10[j] > signif
        
        # aditional metrics for debugging: difference between the binomial and poisson estimate, and estimated variance using a binomial model
        samples_error_comp[j] = ((samples_error_binom[j]-samples_error_pois[j])/(samples_error_binom[j]+samples_error_pois[j]))^2
        samples_var[j] = prob_event*(1-prob_event)*total_points
      }
    }
    
    # avg number of points per sample (and var)
    scales_count[i] = mean(samples_count[!is.na(samples_count)])
    scales_var[i] = var(samples_count[!is.na(samples_count)])
    
    # mean and median for the nearest neighbor p-value test
    scales_nn_p[i] = mean(samples_nn_p[!is.na(samples_nn_p)])
    scales_nn_median_p[i] = median(samples_nn_p[!is.na(samples_nn_p)])
    
    # proportion of samples that passed the threshold for CSR using the nearest neighbor approach
    scales_nn_pass[i] = mean(samples_nn_pass[!is.na(samples_nn_pass)])
    
    # mean and median for the p-value of the quadrat test
    scales_quadrat_p[i] = mean(samples_quadrat_p[!is.na(samples_quadrat_p)])
    scales_quadrat_median_p[i] = median(samples_quadrat_p[!is.na(samples_quadrat_p)])
    
    # proportion of samples that passed the threshold for CSR using the quadrat test
    scales_quadrat_pass[i] = mean(samples_quadrat_pass[!is.na(samples_quadrat_pass)])
    
    # mean values for the coef_of var
    scales_error_binom[i] = mean(samples_error_binom[!is.na(samples_error_binom)])
    scales_error_pois[i] = mean(samples_error_pois[!is.na(samples_error_pois)])
    scales_error_boot[i] = mean(samples_error_boot[!is.na(samples_error_boot)])
    scales_error_comp[i] = mean(samples_error_comp[!is.na(samples_error_comp)])
    scales_new_robust10[i] = mean(samples_new_robust10[!is.na(samples_new_robust10)])
    scales_new_robust25[i] = mean(samples_new_robust25[!is.na(samples_new_robust25)])
    scales_new_robust5[i] = mean(samples_new_robust5[!is.na(samples_new_robust5)])
    scales_robust_pass[i] = mean(samples_robust_pass[!is.na(samples_robust_pass)])
    
    # number of zeros
    scales_zeros[i] = samples_zeros/iter2
  }
  ret_list = list(count = scales_count,
                  var = scales_var,
                  nn_p = scales_nn_p,
                  nn_median_p = scales_nn_median_p,
                  nn_pass = scales_nn_pass,
                  quadrat_p = scales_quadrat_p,
                  quadrat_median_p = scales_quadrat_median_p,
                  quadrat_pass = scales_quadrat_pass,
                  error_binom = scales_error_binom,
                  error_pois = scales_error_pois,
                  error_boot = scales_error_boot,
                  error_comp = scales_error_comp,
                  new_robust10 = scales_new_robust10,
                  new_robust25 = scales_new_robust25,
                  new_robust5 = scales_new_robust5,
                  robust_pass = scales_robust_pass,
                  zerors = scales_zeros
  )
  return(ret_list)
}

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