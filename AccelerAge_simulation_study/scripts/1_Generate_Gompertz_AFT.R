
################################################################
#### Set-up
################################################################

library(survival)
library(scam)
source("src/gompertz_draw.R")

set.seed(336) # run entire script at once for reproducible results

################################################################
#### Parameters
################################################################

a <- exp(-9)
b <- 0.085
sigma <- 1/b
tau <- a/b

followup = 20
M <- 2 # number of predictors (all standard normally distributed)
betas <- c(0.05, 0.05) # vector of beta's (should be of length M)

N_pop <- 1e5 # population size (from which lifetable is determined)
n_test <- 5e3 # number of samples in test data set 
n_obs <- c(500, 2500, 5000, 7500, 10000) # number of observations 
n_sim <- 200 # number of iterations for each n_obs scenario

params <- list(a = a, b = b, sigma = sigma, tau = tau,
               followup = followup, M = M, betas = betas,
               N_pop = N_pop, n_obs = n_obs, n_sim = n_sim)
save(params, file = "output/scenarios/Gompertz_AFT/simulated/params")

################################################################
#### Create lifetable for population
################################################################

X <- matrix( rnorm(N_pop*M,mean=0,sd=1), N_pop, M) # matrix of predictors 

t <- vector(length = N_pop)                                                
for (i in 1:N_pop){
  t[i] <- rgompertz_aft(1, sigma = sigma, tau = tau, linpred = sum(betas * X[i,])) 
}

lifetable_pop <- as.data.frame(cbind(X, t))
lifetable_pop <- lifetable_pop[order(lifetable_pop$t),]

medrl <- vector(length = N_pop)
for (j in 1:N_pop){
  medrl[j] <- median(lifetable_pop$t[j:N_pop]) - lifetable_pop$t[j]
}

lifetable_pop$medrl <- medrl

# smoothen 
fit4 <- scam(medrl ~ s(t, bs = "mpd"), data = lifetable_pop)
xx <- seq(0,max(lifetable_pop$t), by = 0.1)
lt <- as.data.frame(cbind(t = xx, medrl =  predict(fit4, data.frame(t=xx))))

save(lt, file = "output/scenarios/Gompertz_AFT/simulated/lifetable")

################################################################
#### Create test data set 
################################################################

n_gen <- 2 * n_test # twice as many to ensure I generate enough, because for some T < C => not observed

X <- matrix( rnorm(n_gen*M,mean=0,sd=1), n_gen, M) 
c <- runif(n_gen, 20, 80)
linpred <- rowSums(sweep(X, 2, betas, "*"))

# Get age of death
age_death <- vector(length = n_gen)
for (i in 1:n_gen){
  age_death[i] <- rgompertz_aft(1, sigma = sigma, linpred = linpred[i], tau = tau)
}

# Remove observations that are left-truncated
indx_obs <- which(c < age_death)[1:n_test]  # [1:n_test] to obtain the intended sample size 
df_sim <- as.data.frame(cbind(X, age_death, c, linpred))[indx_obs,]

# Get median residual life
medrl <- vector(length = nrow(df_sim))

for (i in 1:nrow(df_sim)){
  surv_prob <- gomp_baseline_surv(df_sim$c[i] * exp(df_sim$linpred[i]), a = a, b = b)
  med_prob <- surv_prob / 2
  t_unadj <- inverse_gomp_baseline_surv(med_prob, a=a, b=b)
  t_adj <- t_unadj / exp(df_sim$linpred[i])
  medrl[i] <-  t_adj - df_sim$c[i]
}
df_sim$medrl <- medrl

# Get biological age (via population lifetable)
bio_age <- vector(length = nrow(df_sim))

for (i in 1:nrow(df_sim)){
  bio_age[i] <- lt$t[ which.min(abs(lt$medrl - df_sim$medrl[i])) ]
}
df_sim$b <- bio_age

# Add censoring 
df_sim$yrs_rem <- df_sim$age_death - df_sim$c
wh <- which(df_sim$yrs_rem > followup) 
df_sim$status <- 1
df_sim$status[wh] <- 0
df_sim$follow_up_time <- df_sim$yrs_rem
df_sim$follow_up_time[wh] <- followup
df_sim$age_end <- df_sim$c + df_sim$follow_up_time

df_test <- df_sim
save(df_test, file = "output/scenarios/Gompertz_AFT/simulated/df_test")

################################################################
#### Generate sample 
################################################################

for (k in 1:length(n_obs)){ # looping over n_obs (number of observations within each dataset)

  list_df_iter <- list()
  n_gen <- 2 * n_obs[k] # twice as many to ensure I generate enough, because for some T < C => not observed
  
  for (j in 1:n_sim){ # looping over n_sim (number of iterations)
  
    X <- matrix( rnorm(n_gen*M,mean=0,sd=1), n_gen, M) 
    c <- runif(n_gen, 20, 80)
    linpred <- rowSums(sweep(X, 2, betas, "*"))
    
    # Get age of death
    age_death <- vector(length = n_gen)
    for (i in 1:n_gen){
      age_death[i] <- rgompertz_aft(1, sigma = sigma, linpred = linpred[i], tau = tau)
    }
    
    # Remove observations that are left-truncated
    indx_obs <- which(c < age_death)[1:n_obs[k]]  # [1:n_obs[k]] to obtain the intended sample size 
    df_sim <- as.data.frame(cbind(X, age_death, c, linpred))[indx_obs,]
    
    # Get median residual life
    medrl <- vector(length = nrow(df_sim))
    
    for (i in 1:nrow(df_sim)){
      surv_prob <- gomp_baseline_surv(df_sim$c[i] * exp(df_sim$linpred[i]), a = a, b = b)
      med_prob <- surv_prob / 2
      t_unadj <- inverse_gomp_baseline_surv(med_prob, a=a, b=b)
      t_adj <- t_unadj / exp(df_sim$linpred[i])
      medrl[i] <-  t_adj - df_sim$c[i]
    }
    df_sim$medrl <- medrl
    
    # Get biological age (via population lifetable)
    bio_age <- vector(length = nrow(df_sim))
    
    for (i in 1:nrow(df_sim)){
      bio_age[i] <- lt$t[ which.min(abs(lt$medrl - df_sim$medrl[i])) ]
    }
    df_sim$b <- bio_age
    
    # Add censoring 
    df_sim$yrs_rem <- df_sim$age_death - df_sim$c
    wh <- which(df_sim$yrs_rem > followup) 
    df_sim$status <- 1
    df_sim$status[wh] <- 0
    df_sim$follow_up_time <- df_sim$yrs_rem
    df_sim$follow_up_time[wh] <- followup
    df_sim$age_end <- df_sim$c + df_sim$follow_up_time
    
    list_df_iter[[j]] <- df_sim
  
  } 
  
  save(list_df_iter, file = paste0("output/scenarios/Gompertz_AFT/simulated/list_df_iter-n_obs", n_obs[k]))

}
