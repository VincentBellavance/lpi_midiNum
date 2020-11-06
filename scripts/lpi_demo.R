## LPI demonstration

# load required packages
library(tidyverse)
library(mgcv)
library(boot)

# load supporting functions
source("scripts/lpi_functions.R")

# set ggplot theme
theme_set(theme_linedraw() +
            theme(panel.grid = element_blank()))

# load population time series
df_l <- readRDS("outputs/timeseries_long.RDS")

# save a vector of the years
years <- unique(df_l$year)

# plot time series
ggplot(df_l) +
  geom_line(aes(x = year, y = N, col = population)) +
  facet_grid(taxa ~ system) +
  coord_cartesian(ylim = c(0,max(df_l$N+10))) +
  labs(y = "N(t)", 
       col = "Populations", 
       title = "Step 1: Population time series") +
  theme(legend.position = "bottom")

# Fit model on each population's abundance through time ------------------------

# log-transform population abundance
df_l$N_log10 <- log10(df_l$N)
pops <- split(df_l, df_l$population)

# create list to store models
m <- vector("list", length(pops))
# fit models
for(i in 1:length(gams)){
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
  labs(y = "N(t) (log10)", 
       col = "Abundances", 
       fill = "GAM predictions",
       title = "Step 2: Predict abundances ~ time from GAM",
       caption = "Ribbon shows standard error from GAM predictions.") + 
  facet_grid(taxa ~ system) +
  theme(legend.position = "bottom")


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
  labs(y = "Annual growth rate", 
       col = "Populations", 
       title = "Step 3: Calculate each population's annual growth rate") +
  facet_grid(taxa ~ system) +
  theme(legend.position = "bottom")


# Take average annual growth rate per taxonomic group --------------------------

# subset to mammals and split by year
mammals <- df_l[which(df_l$taxa == "mammals"),]

# subset to birds and split by year
birds <- df_l[which(df_l$taxa == "birds"),]

# create df to store results
mean_mammals <- data.frame(dt_mean = 1, cilo = 1, cihi = 1)
mean_birds <- data.frame(dt_mean = 1, cilo = 1, cihi = 1)

# calculate geometric mean growth rate 
# with 95% CI from bootstrapping
for(t in 2:length(years)){
  
  # birds, marine
  dt_birds_marine <- filter(birds, year == years[t]) %>% 
    filter(system == "marine") %>% 
    na.omit()
  birds_marine[t,] <- dt_boot(10^dt_birds$dt) 
  
  
  # mammals
  dt_mammals <- filter(birds, year == years[t]) %>% na.omit()
  mean_mammals[t,] <- dt_boot(10^dt_mammals$dt) 
  
}

# dt_gm = bind_rows(dt_gm) %>% log10() %>% 
#   mutate(time = unique(dt_df$time))

# Calculate index based on baseline (LPI = 1 in 1970) --------------------------

