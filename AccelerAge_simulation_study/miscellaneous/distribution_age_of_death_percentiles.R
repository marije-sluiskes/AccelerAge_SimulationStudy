################################################################
#### Set-up
################################################################

library(survival)
library(scam)
library(dplyr)
library(ggplot2)
library(gridExtra)
source("src/gompertz_draw.R")
source("src/weibull_draw.R")

set.seed(336) # run entire script at once for reproducible results

################################################################
#### Parameters
################################################################

a <- exp(-9)
b <- 0.085
sigma <- 1/b
tau <- a/b

lambda = 34^-10
nu = 8

M <- 2 # number of predictors (all standard normally distributed)
n_gen<- 10e4 # number of samples in test data set 

################################################################
#### Gompertz PH 
################################################################

betas <- c(0.3, 0.3) # vector of beta's (should be of length M)

X <- matrix( rnorm(n_gen*M,mean=0,sd=1), n_gen, M) 
c <- rep(0, n_gen)
linpred <- rowSums(sweep(X, 2, betas, "*"))

# Get age of death
age_death <- vector(length = n_gen)
for (i in 1:n_gen){
  age_death[i] <- rgompertz_ph(1, a= a, b = b, linpred = linpred[i])         
}

df_sim <- as.data.frame(cbind(X, age_death, c, linpred))

p.ph <- df_sim %>% 
  mutate(linpred_fac = cut_number(linpred, 5)) %>% 
  ggplot(aes(x = age_death)) +
  geom_density(aes(col = linpred_fac)) +
  scale_color_manual(values = c("springgreen4", "springgreen1", "gold2", "orange", "red"), labels = c("0-20", "20-40", "40-60", "60-80", "80-100")) +
  labs(title = "Gompertz PH", 
       x = "chronological age", 
       col = "Percentile") +
  coord_cartesian(xlim = c(0, 110), ylim = c(0,0.04)) +
  theme_bw()

################################################################
#### Gompertz AFT 
################################################################

betas <- c(0.05, 0.05) # vector of beta's (should be of length M)

X <- matrix( rnorm(n_gen*M,mean=0,sd=1), n_gen, M) 
c <- runif(n_gen, 20, 80)
linpred <- rowSums(sweep(X, 2, betas, "*"))

# Get age of death
age_death <- vector(length = n_gen)
for (i in 1:n_gen){
  age_death[i] <- rgompertz_aft(1, sigma = sigma, linpred = linpred[i], tau = tau)
}

df_sim <- as.data.frame(cbind(X, age_death, c, linpred))

p.aft <- df_sim %>% 
  mutate(linpred_fac = cut_number(linpred, 5)) %>% 
  ggplot(aes(x = age_death)) +
  geom_density(aes(col = linpred_fac)) +
  scale_color_manual(values = c("springgreen4", "springgreen1", "gold2", "orange", "red"), labels = c("0-20", "20-40", "40-60", "60-80", "80-100")) +
  labs(title = "Gompertz AFT", 
       x = "chronological age", 
       col = "Percentile") +
  coord_cartesian(xlim = c(0, 110), ylim = c(0,0.04)) +
  theme_bw()

################################################################
#### Weibull
################################################################

betas <- c(0.35,0.35) # vector of beta's (should be of length M)

X <- matrix( rnorm(n_gen*M,mean=0,sd=1), n_gen, M) 
c <- runif(n_gen, 20, 80)
linpred <- rowSums(sweep(X, 2, betas, "*"))

# Get age of death
age_death <- vector(length = n_gen)                                       
for (i in 1:n_gen){
  age_death[i] <- rweibull(1, lambda = lambda, nu = nu, linpred = linpred[i])
}

df_sim <- as.data.frame(cbind(X, age_death, c, linpred))

p.weib <- df_sim %>% 
  mutate(linpred_fac = cut_number(linpred, 5)) %>% 
  ggplot(aes(x = age_death)) +
  geom_density(aes(col = linpred_fac)) +
  scale_color_manual(values = c("springgreen4", "springgreen1", "gold2", "orange", "red"), labels = c("0-20", "20-40", "40-60", "60-80", "80-100")) +
  labs(title = "Weibull", 
       x = "chronological age", 
       col = "Percentile") +
  coord_cartesian(xlim = c(0, 110), ylim = c(0,0.04)) +
  theme_bw()

p.all <- grid.arrange(p.ph, p.aft, p.weib)
ggsave("output/quantile.png", p.all, width = 6, height =6, units = "in", dpi = 600)
ggsave("output/quantile.pdf", p.all, width = 6, height =6)
