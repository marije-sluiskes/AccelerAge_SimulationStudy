
################################################################
#### Set-up
################################################################

library(ggplot2)
library(dplyr)

################################################################
#### Scenario 1: Gompertz AFT 
################################################################

# load data 
load("output/scenarios/Gompertz_AFT/simulated/params")
load("output/scenarios/Gompertz_AFT/simulated/lifetable")

list_rmse_aft_gomp <- list_rmse_grimage <- list_rmse_coxph <- list_rmse_aft_weib <- list_rmse_aft_semipar <- list_rmse_aft_flexpar <- list()

for (k in 1:length(params$n_obs)) {
  
  load(paste0("output/scenarios/Gompertz_AFT/analyzed/list_all_results-n_obs", params$n_obs[k]))
  load(paste0("output/scenarios/Gompertz_AFT/analyzed/no_err", params$n_obs[k]))
  
  # remove iterations  that threw an error, see explanation script 2 at Gompertz AFT
  list_all_results <- list_all_results[no_err]
  
  remove <- vector(length = length(list_all_results))
  # remove iterations with very high error
  for (i in 1:length(list_all_results)){
    if(list_all_results[[i]]$aft_gomp_res$rmse > 20){ # removing unstable cases, see explanation script 2 at Gompertz AFT 
      remove[i] <- TRUE
    }
  }
  
  print(paste(params$n_obs[k], ":", sum(remove)))
  
  list_all_results <- list_all_results[!remove]
  
  rmse_aft_gomp <- rmse_grimage <- rmse_coxph <- rmse_aft_weib <- rmse_aft_semipar <- rmse_aft_flexpar <- vector(length = length(list_all_results))
  
  for (i in 1:length(list_all_results)){
    
    rmse_aft_gomp[i] <- list_all_results[[i]]$aft_gomp_res$rmse
    rmse_grimage[i] <- list_all_results[[i]]$grimage_res$rmse
    rmse_coxph[i] <- list_all_results[[i]]$coxph_res$rmse
    rmse_aft_weib[i] <- list_all_results[[i]]$aft_weib_res$rmse
    rmse_aft_semipar[i] <- list_all_results[[i]]$aft_semipar_res$rmse
    rmse_aft_flexpar[i] <- list_all_results[[i]]$aft_flexpar_res$rmse
  }
  
  list_rmse_aft_gomp[[k]] <- rmse_aft_gomp
  list_rmse_grimage[[k]] <- rmse_grimage
  list_rmse_coxph[[k]] <- rmse_coxph
  list_rmse_aft_weib[[k]] <- rmse_aft_weib
  list_rmse_aft_semipar[[k]] <- rmse_aft_semipar
  list_rmse_aft_flexpar[[k]] <- rmse_aft_flexpar
}

names(list_rmse_aft_gomp) = names(list_rmse_grimage) = names(list_rmse_coxph) = names(list_rmse_aft_weib) = names(list_rmse_aft_semipar) = names(list_rmse_aft_flexpar) <- params$n_obs

df_rmse_aft_gomp <- as.data.frame(cbind(rmse_mean = sapply(list_rmse_aft_gomp, mean), rmse_sd = sapply(list_rmse_aft_gomp, sd), n_obs = params$n_obs, method = "AFT_Gompertz-mrl"))
df_rmse_grimage <- as.data.frame(cbind(rmse_mean = sapply(list_rmse_grimage, mean), rmse_sd = sapply(list_rmse_grimage, sd), n_obs = params$n_obs, method = "GrimAge"))
df_rmse_coxph <- as.data.frame(cbind(rmse_mean = sapply(list_rmse_coxph, mean), rmse_sd = sapply(list_rmse_coxph, sd), n_obs = params$n_obs, method = "Cox_PH-mrl"))
df_rmse_aft_weib <- as.data.frame(cbind(rmse_mean = sapply(list_rmse_aft_weib, mean), rmse_sd = sapply(list_rmse_aft_weib, sd), n_obs = params$n_obs, method = "AFT_Weibull-mrl"))
df_rmse_aft_semipar <- as.data.frame(cbind(rmse_mean = sapply(list_rmse_aft_semipar, mean), rmse_sd = sapply(list_rmse_aft_semipar, sd), n_obs = params$n_obs, method = "AFT_semiparametric-mrl"))
df_rmse_aft_flexpar <- as.data.frame(cbind(rmse_mean = sapply(list_rmse_aft_flexpar, mean), rmse_sd = sapply(list_rmse_aft_flexpar, sd), n_obs = params$n_obs, method = "AFT_flexibleparametric-mrl"))

df_rmse_all <- rbind(df_rmse_aft_gomp, df_rmse_grimage, df_rmse_coxph)
df_rmse_all$n_obs <- as.numeric(df_rmse_all$n_obs)
df_rmse_all$rmse_mean <- as.numeric(df_rmse_all$rmse_mean)
df_rmse_all$rmse_sd <- as.numeric(df_rmse_all$rmse_sd)

# plot 
p.gomp_aft <- ggplot(df_rmse_all, aes(x = n_obs, y = rmse_mean, color = method)) +
  geom_line(linewidth = 1) + 
  geom_point() +
  ggtitle("Gompertz-AFT") +
  labs(x = "number of observations", y = "root-mean-square error", color = "Predictor") +
  coord_cartesian(ylim = c(0, 3)) +
  scale_color_manual("Predictor", values = c( "#56B4E9", "#CC79A7", "#000000"),
                     labels = c("AccelerAge-Gompertz", "PH-semipar", "GrimAge-type procedure")) + # color-blind friendly 
  theme_bw() 
p.gomp_aft

# save plot
ggsave("output/scenarios/Gompertz_AFT/plots/Gompertz_AFT.pdf", p.gomp_aft, width = 7, height = 5)
ggsave("output/scenarios/Gompertz_AFT/plots/Gompertz_AFT.png", p.gomp_aft, width = 7, height = 5, units = "in", dpi = 600)



################################################################
#### Scenario 2: Gompertz PH
################################################################

# load data 
load("output/scenarios/Gompertz_PH/simulated/params")
load("output/scenarios/Gompertz_PH/simulated/lifetable")

list_rmse_aft_gomp <- list_rmse_grimage <- list_rmse_coxph <- list_rmse_aft_weib <- list_rmse_aft_semipar <- list_rmse_aft_flexpar <- list()

for (k in 1:length(params$n_obs)) {
  
  load(paste0("output/scenarios/Gompertz_PH/analyzed/list_all_results-n_obs", params$n_obs[k]))
  load(paste0("output/scenarios/Gompertz_PH/analyzed/no_err", params$n_obs[k]))
  
  # remove iterations  that threw an error, see explanation script 2 at Gompertz AFT
  list_all_results <- list_all_results[no_err]
  
  remove <- vector(length = length(list_all_results))
  # remove iterations with very high error
  for (i in 1:length(list_all_results)){
    if(list_all_results[[i]]$aft_gomp_res$rmse > 20){ # removing unstable cases, see explanation script 2 at Gompertz AFT 
      remove[i] <- TRUE
    }
  }
  
  print(paste(params$n_obs[k], ":", sum(remove)))
  
  list_all_results <- list_all_results[!remove]
  
  rmse_aft_gomp <- rmse_grimage <- rmse_coxph <- rmse_aft_weib <- rmse_aft_semipar <- rmse_aft_flexpar <- vector(length = length(list_all_results))
  
  for (i in 1:length(list_all_results)){
    
    rmse_aft_gomp[i] <- list_all_results[[i]]$aft_gomp_res$rmse
    rmse_grimage[i] <- list_all_results[[i]]$grimage_res$rmse
    rmse_coxph[i] <- list_all_results[[i]]$coxph_res$rmse
    rmse_aft_weib[i] <- list_all_results[[i]]$aft_weib_res$rmse
    rmse_aft_semipar[i] <- list_all_results[[i]]$aft_semipar_res$rmse
    rmse_aft_flexpar[i] <- list_all_results[[i]]$aft_flexpar_res$rmse
  }
  
  list_rmse_aft_gomp[[k]] <- rmse_aft_gomp
  list_rmse_grimage[[k]] <- rmse_grimage
  list_rmse_coxph[[k]] <- rmse_coxph
  list_rmse_aft_weib[[k]] <- rmse_aft_weib
  list_rmse_aft_semipar[[k]] <- rmse_aft_semipar
  list_rmse_aft_flexpar[[k]] <- rmse_aft_flexpar
}

names(list_rmse_aft_gomp) = names(list_rmse_grimage) = names(list_rmse_coxph) = names(list_rmse_aft_weib) = names(list_rmse_aft_semipar) = names(list_rmse_aft_flexpar) <- params$n_obs

df_rmse_aft_gomp <- as.data.frame(cbind(rmse_mean = sapply(list_rmse_aft_gomp, mean), rmse_sd = sapply(list_rmse_aft_gomp, sd), n_obs = params$n_obs, method = "AFT_Gompertz-mrl"))
df_rmse_grimage <- as.data.frame(cbind(rmse_mean = sapply(list_rmse_grimage, mean), rmse_sd = sapply(list_rmse_grimage, sd), n_obs = params$n_obs, method = "GrimAge"))
df_rmse_coxph <- as.data.frame(cbind(rmse_mean = sapply(list_rmse_coxph, mean), rmse_sd = sapply(list_rmse_coxph, sd), n_obs = params$n_obs, method = "Cox_PH-mrl"))
df_rmse_aft_weib <- as.data.frame(cbind(rmse_mean = sapply(list_rmse_aft_weib, mean), rmse_sd = sapply(list_rmse_aft_weib, sd), n_obs = params$n_obs, method = "AFT_Weibull-mrl"))
df_rmse_aft_semipar <- as.data.frame(cbind(rmse_mean = sapply(list_rmse_aft_semipar, mean), rmse_sd = sapply(list_rmse_aft_semipar, sd), n_obs = params$n_obs, method = "AFT_semiparametric-mrl"))
df_rmse_aft_flexpar <- as.data.frame(cbind(rmse_mean = sapply(list_rmse_aft_flexpar, mean), rmse_sd = sapply(list_rmse_aft_flexpar, sd), n_obs = params$n_obs, method = "AFT_flexibleparametric-mrl"))

df_rmse_all <- rbind(df_rmse_aft_gomp, df_rmse_grimage, df_rmse_coxph)
df_rmse_all$n_obs <- as.numeric(df_rmse_all$n_obs)
df_rmse_all$rmse_mean <- as.numeric(df_rmse_all$rmse_mean)
df_rmse_all$rmse_sd <- as.numeric(df_rmse_all$rmse_sd)

# plot 
p.gomp_ph <- ggplot(df_rmse_all, aes(x = n_obs, y = rmse_mean, color = method)) +
  geom_line(linewidth = 1) + 
  geom_point() +
  ggtitle("Gompertz-PH") +
  labs(x = "number of observations", y = "root-mean-square error", color = "Predictor") +
  coord_cartesian(ylim = c(0, 3)) +
  scale_color_manual("Predictor", values = c("#56B4E9",  "#CC79A7", "#000000"),
                     labels = c( "AccelerAge-Gompertz", "PH-semipar", "GrimAge-type procedure")) + # color-blind friendly 
  theme_bw() 
p.gomp_ph

# save plot
ggsave("output/scenarios/Gompertz_PH/plots/Gompertz_PH.png", p.gomp_ph,  width = 7, height = 5, units = "in", dpi = 600)
ggsave("output/scenarios/Gompertz_PH/plots/Gompertz_PH.pdf", p.gomp_ph,  width = 7, height = 5)


################################################################
#### Scenario 3: Weibull
################################################################

# load data 
load("output/scenarios/Weibull/simulated/params")
load("output/scenarios/Weibull/simulated/lifetable")


list_rmse_aft_gomp <- list_rmse_grimage <- list_rmse_coxph <- list_rmse_aft_weib <- list_rmse_aft_semipar <- list_rmse_aft_flexpar <- list()

for (k in 1:length(params$n_obs)) {
  
  load(paste0("output/scenarios/Weibull/analyzed/list_all_results-n_obs", params$n_obs[k]))
  load(paste0("output/scenarios/Weibull/analyzed/no_err", params$n_obs[k]))
  
  # remove iterations  that threw an error, see explanation script 2 at Gompertz AFT
  list_all_results <- list_all_results[no_err]
  
  remove <- vector(length = length(list_all_results))
  # remove iterations with very high error
  for (i in 1:length(list_all_results)){
    if(list_all_results[[i]]$aft_gomp_res$rmse > 20){ # removing unstable cases, see explanation script 2 at Gompertz AFT 
      remove[i] <- TRUE
    }
  }
  
  list_all_results <- list_all_results[!remove]
  
  print(paste(params$n_obs[k], ":", sum(remove)))
  
  rmse_aft_gomp <- rmse_grimage <- rmse_coxph <- rmse_aft_weib <- rmse_aft_semipar <- rmse_aft_flexpar <- vector(length = length(list_all_results))
  
  for (i in 1:length(list_all_results)){
    
    rmse_aft_gomp[i] <- list_all_results[[i]]$aft_gomp_res$rmse
    rmse_grimage[i] <- list_all_results[[i]]$grimage_res$rmse
    rmse_coxph[i] <- list_all_results[[i]]$coxph_res$rmse
    rmse_aft_weib[i] <- list_all_results[[i]]$aft_weib_res$rmse
    rmse_aft_semipar[i] <- list_all_results[[i]]$aft_semipar_res$rmse
    rmse_aft_flexpar[i] <- list_all_results[[i]]$aft_flexpar_res$rmse
  }
  
  list_rmse_aft_gomp[[k]] <- rmse_aft_gomp
  list_rmse_grimage[[k]] <- rmse_grimage
  list_rmse_coxph[[k]] <- rmse_coxph
  list_rmse_aft_weib[[k]] <- rmse_aft_weib
  list_rmse_aft_semipar[[k]] <- rmse_aft_semipar
  list_rmse_aft_flexpar[[k]] <- rmse_aft_flexpar
}

names(list_rmse_aft_gomp) = names(list_rmse_grimage) = names(list_rmse_coxph) = names(list_rmse_aft_weib) = names(list_rmse_aft_semipar) = names(list_rmse_aft_flexpar) <- params$n_obs

df_rmse_aft_gomp <- as.data.frame(cbind(rmse_mean = sapply(list_rmse_aft_gomp, mean), rmse_sd = sapply(list_rmse_aft_gomp, sd), n_obs = params$n_obs, method = "AFT_Gompertz-mrl"))
df_rmse_grimage <- as.data.frame(cbind(rmse_mean = sapply(list_rmse_grimage, mean), rmse_sd = sapply(list_rmse_grimage, sd), n_obs = params$n_obs, method = "GrimAge"))
df_rmse_coxph <- as.data.frame(cbind(rmse_mean = sapply(list_rmse_coxph, mean), rmse_sd = sapply(list_rmse_coxph, sd), n_obs = params$n_obs, method = "Cox_PH-mrl"))
df_rmse_aft_weib <- as.data.frame(cbind(rmse_mean = sapply(list_rmse_aft_weib, mean), rmse_sd = sapply(list_rmse_aft_weib, sd), n_obs = params$n_obs, method = "AFT_Weibull-mrl"))
df_rmse_aft_semipar <- as.data.frame(cbind(rmse_mean = sapply(list_rmse_aft_semipar, mean), rmse_sd = sapply(list_rmse_aft_semipar, sd), n_obs = params$n_obs, method = "AFT_semiparametric-mrl"))
df_rmse_aft_flexpar <- as.data.frame(cbind(rmse_mean = sapply(list_rmse_aft_flexpar, mean), rmse_sd = sapply(list_rmse_aft_flexpar, sd), n_obs = params$n_obs, method = "AFT_flexibleparametric-mrl"))

df_rmse_all <- rbind(df_rmse_aft_gomp, df_rmse_grimage, df_rmse_coxph)
df_rmse_all$n_obs <- as.numeric(df_rmse_all$n_obs)
df_rmse_all$rmse_mean <- as.numeric(df_rmse_all$rmse_mean)
df_rmse_all$rmse_sd <- as.numeric(df_rmse_all$rmse_sd)

# plot 
p.weib <- ggplot(df_rmse_all, aes(x = n_obs, y = rmse_mean, color = method)) +
  geom_line(linewidth = 1) + 
  geom_point() +
  ggtitle("Weibull") +
  labs(x = "number of observations", y = "root-mean-square error", color = "Predictor") +
  coord_cartesian(ylim = c(0, 3)) +
  scale_color_manual("Predictor", values = c("#56B4E9",  "#CC79A7", "#000000"),
                     labels = c( "AccelerAge-Gompertz", "PH-semipar", "GrimAge-type procedure")) + # color-blind friendly 
  theme_bw() 
p.weib

# save plot
ggsave("output/scenarios/Weibull/plots/Weibull.png", p.weib, width = 7, height = 5, units = "in", dpi = 600)
ggsave("output/scenarios/Weibull/plots/Weibull.pdf", p.weib, width = 7, height = 5)


