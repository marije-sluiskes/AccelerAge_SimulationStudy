
## One plot, one iteration

library(dplyr)
library(ggplot2)
library(gridExtra)

### Gompertz AFT

load("output/scenarios/Gompertz_AFT/simulated/params")
load("output/scenarios/Gompertz_AFT/simulated/lifetable")
load("output/scenarios/Gompertz_AFT/simulated/df_test")

p.pop.gomp_aft <- df_test %>% 
  ggplot(aes(x = c, y = b)) +
  geom_point(alpha = 0.25) +
  geom_abline(slope = 1, col = "orange", linewidth = 1) + 
  labs(title = "Gompertz-AFT",
       x = "Chronological age",
       y = "Biological age") +
  coord_cartesian(ylim = c(0, 110)) +
  theme_bw()
p.pop.gomp_aft

ggsave("output/scenarios/Gompertz_AFT/plots/Gompertz_AFT_population.png", p.pop.gomp_aft, width = 7, height =5, units = "in", dpi = 600)
ggsave("output/scenarios/Gompertz_AFT/plots/Gompertz_AFT_population.pdf", p.pop.gomp_aft, width = 7, height =5)

p.bioage.example <- lt %>% 
  ggplot(aes(x = t, y = medrl)) +
  geom_line(size = 1) + 
  labs(title = "Population lifetable",
       x = "Chronological age",
       y = "Mean residual life") +
  geom_segment(aes(x = 0, y = 20, xend = 55, yend = 20),
                 arrow = arrow(length = unit(0.5, "cm")), col = "orange") +
  geom_segment(aes(x = 55, y = 20, xend = 55, yend = 0),
               arrow = arrow(length = unit(0.5, "cm")), col = "orange") +
  theme_bw() +
  xlim(c(0, 100))
p.bioage.example

ggsave("output/bioage_exmple.png", p.bioage.example, width = 7, height =5, units = "in", dpi = 600)
ggsave("output/bioage_exmple.pdf", p.bioage.example, width = 7, height =5)

### Gompertz PH

load("output/scenarios/Gompertz_PH/simulated/params")
load("output/scenarios/Gompertz_PH/simulated/lifetable")
load("output/scenarios/Gompertz_PH/simulated/df_test")

p.pop.gomp_ph <- df_test %>% 
  ggplot(aes(x = c, y = b)) +
  geom_point(alpha = 0.25) +
  geom_abline(slope = 1, col = "orange", linewidth = 1) + 
  labs(title = "Gompertz-PH",
       x = "Chronological age",
       y = "Biological age")  +
  coord_cartesian(ylim = c(0, 110)) +
  theme_bw()
p.pop.gomp_ph

ggsave("output/scenarios/Gompertz_AFT/plots/Gompertz_PH_population.jpg", p.pop)


### Weibull AFT/PH

load("output/scenarios/Weibull/simulated/params")
load("output/scenarios/Weibull/simulated/lifetable")
load("output/scenarios/Weibull/simulated/df_test")

p.pop.weib <- df_test %>% 
  ggplot(aes(x = c, y = b)) +
  geom_point(alpha = 0.25) +
  geom_abline(slope = 1, col = "orange", linewidth = 1) + 
  labs(title = "Weibull",
       x = "Chronological age",
       y = "Biological age") +
  coord_cartesian(ylim = c(0, 110)) +
  theme_bw()
p.pop.weib

ggsave("output/scenarios/Gompertz_AFT/plots/Weibull_population.jpg", p.pop)

p.pop.all <- grid.arrange(p.pop.gomp_ph, p.pop.gomp_aft, p.pop.weib, nrow = 3)
ggsave("output/population_all.png", p.pop.all, width = 6, height = 10, units = "in", dpi = 600)
ggsave("output/population_all.pdf", p.pop.all, width = 6, height = 10)

