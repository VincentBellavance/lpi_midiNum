## LPI demonstration

# load required packages
library(tidyverse)
library(mgcv)
library(boot)

# load supporting functions
source("scripts/lpi_functions.R")

# set ggplot theme
theme_set(theme_linedraw() +
            theme(panel.grid = element_blank()
                  #axis.text = element_text(size = 9),
                  #axis.title = element_text(size = 10),
                  #strip.text = element_text(face = "bold", size = 10),
                  #legend.text = element_text(size = 10))
                  )
)

# load population time series
df_l <- readRDS("outputs/timeseries_long.RDS")

# save a vector of the years
years <- unique(df_l$year)

# plot time series
ggplot(df_l) +
  geom_line(aes(x = year, y = N, col = population)) +
  facet_grid(taxa ~ system) +
  coord_cartesian(ylim = c(0,max(df_l$N+10))) +
  labs(y = "Abundance", x = "",
       col = "Populations", 
       title = "Step 1: Population time series") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "none")
ggsave("images/fig1_timeseries.svg", width = 5.54, height = 3.62, units = "in", dpi = 150)

# Fit model on each population's abundance through time ------------------------

# log-transform population abundance
df_l$N_log10 <- log10(df_l$N)
pops <- split(df_l, df_l$population)

# create list to store models
m <- vector("list", length(pops))
# fit models
for(i in 1:length(m)){
  m[[i]] <- gam(N_log10 ~ s(year, k = 10), # k = half the length of the time series
                   data = pops[[i]],
                   family = gaussian(), 
                   fx = TRUE, 
                   method = "REML")
}
names(m) <- names(pops)



# Predict abundance over time period covered by the index ----------------------

pred_ls = lapply(m, predict.gam, type = "response", se.fit = TRUE)

# wrangle into long format
pred = pred_ls %>%
  lapply(bind_cols) %>%
  lapply(mutate, year = unique(df_l$year)) %>%
  bind_rows(.id = "population")
# join to population time series dataframe
df_l = full_join(df_l, pred, by = c("population", "year"))

# plot the predicted trend over the original
ggplot(df_l, aes(x = year, group = population)) +
  # original time series
  geom_line(aes(y = N_log10, col = population), 
            lty = 2, lwd = .3) +
  # GAM prediction
  geom_ribbon(aes(ymin = fit - se.fit,
                  ymax = fit + se.fit, fill = population), alpha = .3) +
  geom_line(aes(y = fit)) +
  labs(y = "Abundances (log10)", x = "",
       col = "Abundances", 
       fill = "GAM predictions",
       title = "Step 2: Predict abundances ~ time from GAM",
       caption = "Ribbon shows standard error from GAM predictions.") + 
  facet_grid(taxa ~ system) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.position = "none")
ggsave("images/fig2_gam.svg", width = 5.54, height = 3.62, units = "in", dpi = 150)


# Calculate each population's growth rate (dt) through time --------------------

# calculate dt from GAM predictions
dt <- lapply(pred_ls, get_dt, time = unique(df_l$year)) %>% 
  bind_rows(.id = "population") %>%
  rename("year" = time)
# join to population time series dataframe
df_l = full_join(df_l, dt, by = c("population", "year"))

# plot the growth rates
ggplot(df_l, aes(x = year, col = population)) +
  geom_line(aes(y = 10^dt)) +
  labs(y = "Annual growth rate", x = "",
       col = "Populations", 
       title = "Step 3: Calculate each population's annual growth rate") +
  facet_grid(taxa ~ system) +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "none")
ggsave("images/fig3_growthrates.svg", width = 5.54, height = 3.62, units = "in", dpi = 150)

# Take average annual growth rate per taxonomic group per system ---------------

# create df to store results
mammals_marine <- data.frame(dt_mean = 1, cilo = 1, cihi = 1)
birds_marine <- data.frame(dt_mean = 1, cilo = 1, cihi = 1)
mammals_terrestrial <- data.frame(dt_mean = 1, cilo = 1, cihi = 1)
birds_terrestrial <- data.frame(dt_mean = 1, cilo = 1, cihi = 1)

# calculate geometric mean growth rate 
# with 95% CI from bootstrapping
for(t in 2:length(years)){
  
  # birds, marine
  dt_temp <- na.omit(df_l) %>% 
    filter(taxa == "birds", system == "marine", year == years[t]) 
  birds_marine[t,] <- dt_boot(10^dt_temp$dt) 
  
  # mammals, marine
  dt_temp <- na.omit(df_l) %>% 
    filter(taxa == "mammals", system == "marine", year == years[t]) 
  mammals_marine[t,] <- dt_boot(10^dt_temp$dt) 
  
  # birds, terrestrial
  dt_temp <- na.omit(df_l) %>% 
    filter(taxa == "birds", system == "terrestrial", year == years[t]) 
  birds_terrestrial[t,] <- dt_boot(10^dt_temp$dt) 
  
  # mammals, terrestrial
  dt_temp <- na.omit(df_l) %>% 
    filter(taxa == "mammals", system == "terrestrial", year == years[t]) 
  mammals_terrestrial[t,] <- dt_boot(10^dt_temp$dt) 
  
}

# join the data frames together
birds_marine <-  mutate(birds_marine, 
                        year = years, taxa = "birds", system = "marine")
mammals_marine <- mutate(mammals_marine, 
                         year = years, taxa = "mammals", system = "marine")
birds_terrestrial <- mutate(birds_terrestrial, 
                            year = years, taxa = "birds", system = "terrestrial")
mammals_terrestrial <- mutate(mammals_terrestrial, 
                              year = years, taxa = "mammals", system = "terrestrial")
taxa_system <- rbind(birds_marine, mammals_marine, birds_terrestrial, mammals_terrestrial)

# plot the taxa_system group means 
ggplot(taxa_system, aes(x = year)) +
  geom_ribbon(aes(ymin = cilo, ymax = cihi, fill = taxa), alpha = .3) +
  geom_line(aes(y = dt_mean)) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(y = "Mean annual growth rate", x = "",
       col = "Taxonomic group", fill = "Taxonomic group", 
       title = "Step 4: Calculate the mean annual growth rate (Part 1)",
       subtitle = "Mean growth rate within each taxonomic group in each system") +
  facet_grid(taxa ~ system) +
  coord_cartesian(ylim = c(0.9, 1.1)) +
  theme(legend.position = "none")
ggsave("images/fig4_mean_taxasystem.svg", width = 5.54, height = 3.62, units = "in", dpi = 150)


# Take mean annual growth rate per system --------------------------------------

# create df to store results
marine <- data.frame(dt_mean = 1, cilo = 1, cihi = 1)
terrestrial <- data.frame(dt_mean = 1, cilo = 1, cihi = 1)

# calculate geometric mean growth rate 
for(t in 1:length(years)){
  
  # marine
  dt_temp <- taxa_system %>% filter(system == "marine", year == years[t]) 
  marine[t,1:3] <- apply(select(dt_temp, c(dt_mean, cilo, cihi)), 2, function(x) gm_mean(10^x))
  
  # terrestrial
  dt_temp <- taxa_system %>% filter(system == "terrestrial", year == years[t]) 
  terrestrial[t,1:3] <- apply(select(dt_temp, c(dt_mean, cilo, cihi)), 2, function(x) gm_mean(10^x))
}

# join the data frames together
marine <-  mutate(marine, year = years, system = "marine")
terrestrial <- mutate(terrestrial,  year = years, system = "terrestrial")
system_df <- rbind(marine, terrestrial)

# plot the system group means 
ggplot(system_df, aes(x = year, group = system)) +
  geom_ribbon(aes(ymin = cilo, ymax = cihi, fill = system), alpha = .3) +
  geom_line(aes(y = dt_mean), lwd = .7) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  labs(y = "Mean annual growth rate", x = "",
       col = "System", 
       title = "Step 4: Calculate the mean annual growth rate (Part 2)",
       subtitle = "Mean growth rate within each system") +
  facet_wrap(~system) +
  theme(legend.position = "none")
ggsave("images/fig5_mean_system.svg", width = 5.54, height = 3.62, units = "in", dpi = 150)

  
# Take average annual growth rate globally -------------------------------------

# create df to store results
all <- data.frame(dt_mean = 1, cilo = 1, cihi = 1)

# calculate geometric mean growth rate 
for(t in 1:length(years)){
    dt_temp <- system_df %>% filter(year == years[t]) 
  all[t,1:3] <- apply(select(dt_temp, c(dt_mean, cilo, cihi)), 2, function(x) gm_mean(log10(x)))
}

# join the data frames together
all <-  mutate(all, year = years)

# plot the system group means 
ggplot(all, aes(x = year)) +
  geom_ribbon(aes(ymin = cilo, ymax = cihi), alpha = .5, fill = "#66a61e") +
  geom_line(aes(y = dt_mean)) +
  labs(y = "Mean annual growth rate", x = "",
       title = "Step 4: Calculate the mean annual growth rate (Part 3)",
       subtitle = "Mean growth rate across systems") +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.9,1.1))
ggsave("images/fig6_mean_all.svg", width = 5.54, height = 3.62, units = "in", dpi = 150)


# Calculate Living Planet Index (LPI baseline = 1) -----------------------------

# calculate LPI with confidence intervals
lpi <- data.frame(
  year = years, 
  lpi = calclpi(log10(all$dt_mean)),
  cilo = calclpi(log10(all$cilo)),
  cihi = calclpi(log10(all$cihi))
)

# plot
ggplot(lpi, aes(x = year)) +
  geom_ribbon(aes(ymin = cilo, ymax = cihi), alpha = .5, fill = "#66a61e") +
  geom_line(aes(y = lpi)) +
  geom_hline(aes(yintercept = 1), lwd = .2, lty = 2) +
  labs(y = "Living Planet Index", x = "",
       title = "Step 5: Calculate the Living Planet Index") +
  coord_cartesian(ylim = c(0,2))
ggsave("images/fig7_lpi.svg", width = 5.54, height = 3.62, units = "in", dpi = 150)
