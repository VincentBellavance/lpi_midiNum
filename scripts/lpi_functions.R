# Supporting functions to demonstrate the calculation of the Living Planet Index 


# function to generate population abundances through time
sim_pop <- function(r, N0, n_years){
  
  # r = annual growth rate
  # N0 = initial populations size
  # n_years = number of years to simulate
  
  # set up vector to hold time series
  N <- matrix(NA, nrow = nsteps, ncol = length(r))
  N[1,] <- N0 # assign initial population value to first time step
  
  # calculate
  for(t in 2:nsteps){  
    Nt = N[c(t-1),p]*10^r[t]  # this is N(t) = N(t-1)*10^dt
    # store each value in matrix
    N[t, p] = Nt
  }
}



# function to calculate annual growth rate per population from GAM predictions
get_dt <- function(gam.pred_ls, time){
  # gam.pred_ls: list of predictions from GAMs (one per population)
  # time: the time of each prediction (i.e. the year)
  
  # extract predicted population size values and reverse log10 transformation
  N = 10^gam.pred_ls$fit
  
  # initialize df to store dts
  dt_df = data.frame(time = time, dt = NA)
  for(i in 2:length(N)){
    # calculate dt
    dt = log10(N[i]/N[i-1])
    dt_df[i, "dt"] = dt # save in the table
  }
  return("dt" = dt_df)
}

# geometric mean function
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# function to calculate geometric mean that will be used in bootstrapping
geoMean_boot <- function(dt, i){
  d <- dt[i] # indexing happens here
  avg <- gm_mean(d, na.rm = TRUE)
  return(avg)
}

# function to get geometric mean via bootstraping 
dt_boot = function(dt){
  # bootstrap resampling
  dt_r = boot::boot(dt, statistic = geoMean_boot, R = 1000)
  # calculate 95% confidence intervals
  dt_ci = boot::boot.ci(dt_r, type = "basic", R = 1000)
  # wrangle into table for output
  dtboot_df = data.frame(gm = dt_ci$t0, 
                         cilo = dt_ci$basic[4], 
                         cihi = dt_ci$basic[5])
  return(dtboot_df)
}
