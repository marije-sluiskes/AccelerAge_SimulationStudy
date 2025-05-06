
################################################################
#### Set-up
################################################################

library(parallel)

# Set the number of cores
no_cores <- 7
#no_cores <- detectCores() - 1

################################################################
#### Prediction function
################################################################

# Function to lapply over (easy to parallelize). Functions can be found in bioage_estimate_median.R.
GetEstimatesList <- function(df_train, df_test, M, lt, df_splines){
  
  train <- df_train
  test <- df_test
  
  # get biological age estimates + RMSE
  aft_gomp_res <- GetB_AFTmrl_Gompertz(train, test, M, lt)
  grimage_res <- GetB_GrimAge(train, test)
  coxph_res <- GetB_CoxPHmrl(train, test, M, lt)
  aft_weib_res <- GetB_AFTmrl_Weibull(train, test, M, lt)
  aft_semipar_res <- GetB_AFTmrl_semipar(train, test, M, lt)
  aft_flexpar_res <- GetB_AFTmrl_flexpar(train, test, M, lt, df_splines)
  
  # make list with df_train, df_test and results from all methods
  list_iter <- list(aft_gomp_res = aft_gomp_res, grimage_res = grimage_res, 
                    coxph_res = coxph_res, aft_weib_res = aft_weib_res, 
                    aft_semipar_res = aft_semipar_res, aft_flexpar_res = aft_flexpar_res)
  return(list_iter)
  
}

# Catch errors (only relevant for Gompertz AFT, due to instability of model-fitting procedure)
tryGetEstimatesList <- function (df_train, df_test, M, lt, df_splines) {
  return(tryCatch(GetEstimatesList(df_train, df_test, M, lt, df_splines), error=function(e) NULL))
}

################################################################
#### Scenario 1: Gompertz AFT 
################################################################

# There is instability in the eha::aftreg(dist = "gompertz") function. 
# Sometimes it returns that the Hessian is singular, in which case I ignore that iteration with tryGetEstimatesList.
# Sometimes it just doesn't properly estimate the parameters (sigma and tau) accurately at all (they blow up completely).
# The author of the eha::aftreg function has confirmed that the AFT Gompertz function is 'numerically delicate to estimate', resulting in these issues.
# I save the indices of the iterations where these issues occur, and ignore them when plotting the results in the next script. 

# load data 
load("output/scenarios/Gompertz_AFT/simulated/params")
load("output/scenarios/Gompertz_AFT/simulated/lifetable")
load("output/scenarios/Gompertz_AFT/simulated/df_test")

# Initiate cluster
cl <- makeCluster(no_cores) 

# Add survival library and own prediction functions to all clusters
clusterCall(cl, function() library(survival))
clusterCall(cl, function() library(eha))
clusterEvalQ(cl, source("src/bioage_estimate_median.R"))
clusterExport(cl, "GetEstimatesList")

for (k in 1:length(params$n_obs)) { # looping over different sample sizes 
  
  load(paste0("output/scenarios/Gompertz_AFT/simulated/list_df_iter-n_obs", params$n_obs[k]))
  list_all_results <- parLapply(cl, list_df_iter, tryGetEstimatesList, M = params$M, lt = lt, df_test = df_test, df_splines = 3)
  
  no_err <- which(sapply(list_all_results, length) != 0)
  
  save(no_err, file = paste0("output/scenarios/Gompertz_AFT/analyzed/no_err", params$n_obs[k]))
  save(list_all_results, file = paste0("output/scenarios/Gompertz_AFT/analyzed/list_all_results-n_obs", params$n_obs[k]))
  
}

stopCluster(cl)

################################################################
#### Scenario 2: Gompertz PH
################################################################

# load data 
load("output/scenarios/Gompertz_PH/simulated/params")
load("output/scenarios/Gompertz_PH/simulated/lifetable")
load("output/scenarios/Gompertz_PH/simulated/df_test")

# Initiate cluster
cl <- makeCluster(no_cores)

# Add survival library and own prediction functions to all clusters
clusterCall(cl, function() library(survival))
clusterCall(cl, function() library(eha))
clusterEvalQ(cl, source("src/bioage_estimate_median.R"))
clusterExport(cl, "GetEstimatesList")

for (k in 1:length(params$n_obs)) { # looping over different sample sizes 
  
  load(paste0("output/scenarios/Gompertz_PH/simulated/list_df_iter-n_obs", params$n_obs[k]))
  list_all_results <- parLapply(cl, list_df_iter, tryGetEstimatesList, M = params$M, lt = lt, df_test = df_test, df_splines = 3)
  
  no_err <- which(sapply(list_all_results, length) != 0)
  
  save(no_err, file = paste0("output/scenarios/Gompertz_PH/analyzed/no_err", params$n_obs[k]))
  save(list_all_results, file = paste0("output/scenarios/Gompertz_PH/analyzed/list_all_results-n_obs", params$n_obs[k]))
  
}

stopCluster(cl)

################################################################
#### Scenario 3: Weibull
################################################################

# load data 
load("output/scenarios/Weibull/simulated/params")
load("output/scenarios/Weibull/simulated/lifetable")
load("output/scenarios/Weibull/simulated/df_test")

# Initiate cluster
cl <- makeCluster(no_cores)

# Add survival library and own prediction functions to all clusters
clusterCall(cl, function() library(survival))
clusterCall(cl, function() library(eha))
clusterEvalQ(cl, source("src/bioage_estimate_median.R"))
clusterExport(cl, "GetEstimatesList")

start_time <- Sys.time()

for (k in 1:length(params$n_obs)) { # looping over different sample sizes 
  
  load(paste0("output/scenarios/Weibull/simulated/list_df_iter-n_obs", params$n_obs[k]))
  list_all_results <- parLapply(cl, list_df_iter, tryGetEstimatesList, M = params$M, lt = lt, df_test = df_test, df_splines = 3)
  
  no_err <- which(sapply(list_all_results, length) != 0)
  
  save(no_err, file = paste0("output/scenarios/Weibull/analyzed/no_err", params$n_obs[k]))
  save(list_all_results, file = paste0("output/scenarios/Weibull/analyzed/list_all_results-n_obs", params$n_obs[k]))
  
}

end_time <- Sys.time()
end_time - start_time

stopCluster(cl)
