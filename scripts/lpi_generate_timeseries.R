# Script to generate population time series for LPI demonstration

# load packages
library(tidyverse)

# load supporting functions
source("scripts/lpi_functions.R")

# Generate population time series ----------------------------------------------

# create data frame of identifier variables
# (ignoring realms for simplicity)
# (assuming there is only 1 population per species)
df_id <- data.frame(
  "system" = c(rep("marine", 4), rep("terrestrial", 4)),
  "taxa" = c(rep(c("birds", "mammals"), 4)),
  "population" = c(LETTERS[1:8])
)

# create data frame to store values
df_pops <- matrix(nrow = 20, ncol = 8, 
               dimnames = list(2000:2019, df_id$population)) %>% 
  as.data.frame()
# generate random abundances 
for(i in 1:ncol(df_pops)){
  df_pops[[i]] <- round(runif(nrow(df_pops), min = 20, max = 40), digits = 0)
}
# transpose
df_pops <- t(df_pops)

# combine the two data frames into one
df <- cbind(df_id, df_pops)

# convert to long format
df_long <- pivot_longer(df, 
                        cols = -c(system, taxa, population),
                        names_to = "year", values_to = "N")
df_long$year <- as.integer(df_long$year)
  
# write file
saveRDS(df, "outputs/timeseries.RDS")
saveRDS(df_long, "outputs/timeseries_long.RDS")
