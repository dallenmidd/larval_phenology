## compare timing to oviposition and ecolosion
## between degree day method and ogden 2005 method


require(tidyverse)
require(bbmle)
require(cowplot)
require(lubridate)
library(mgcv)
set.seed(31417)

leaf_litter_temp <- read_csv('data/leaf_litter_temp_BRF2_separate.csv')

foote <- leaf_litter_temp %>% filter(site == 'Foote')

params_mean <- list(
  ovi_m = 187,
  ovi_sd = 50,
  ecl_m = 533,
  ecl_sd = 28,
  diapause = 0.5,
  hardening = 7, 
  start_quest = 10,
  max_quest = 25,
  host_find = 0.02, 
  mort = 0.01, 
  overwinter_surv = 0.435
)


dd_ovi_day <- function(tmean,params)
{
  tot_day <- length(tmean)
  dd6 <- cumsum(ifelse(tmean > 6, tmean - 6, 0))
  # fraction of cohort oviposited on each day
  ovi <- pnorm(dd6, params$ovi_m, params$ovi_sd) - pnorm(lag(dd6, n = 1, default = 0), params$ovi_m, params$ovi_sd)             
  
  sum(ovi * 1:366)
}


dd_ecl_day <- function(tmean,params)
{
  tot_day <- length(tmean)
  dd6 <- cumsum(ifelse(tmean > 6, tmean - 6, 0))
  # fraction of cohort oviposited on each day
  ovi <- pnorm(dd6, params$ovi_m, params$ovi_sd) - pnorm(lag(dd6, n = 1, default = 0), params$ovi_m, params$ovi_sd)             
  
  # fraction of cohort eclosed on each day
  ecl <- rep(0, tot_day)
  for (i in 1:(tot_day-1))
  {
    i1 <- i+1
    dd11 <- c(rep(0,i),cumsum(ifelse(tmean[i1:tot_day]>11,tmean[i1:tot_day]-11,0)))
    # ecl_temp is prob of ecolsion on each day if oviposition is on day i
    ecl_prob <- pnorm(dd11,params$ecl_m,params$ecl_sd) - pnorm(lag(dd11,n=1,default = 0),params$ecl_m,params$ecl_sd)
    ecl <- ecl + ovi[i] * ecl_prob
  }
  
  sum(ecl * 1:366)/sum(ecl)
}



#### Not all ticks ecl by end of season at high elev
#### Not sure what happens to them in my exiting model
#### Might need to add them to next season ticks or maybe alreay doing that



#### Now need to do it for odgen exponential method

## from Odgen 2004
# days to oviposition = 1300*temp^-1.42
# so daily development = 1/1300 * temp^1.42

# days to eclosion = 34234*temp^-2.27
# daily deve = 1/34234*temp^2.27

ogd_ovi_day <- function(tmean)
{
  ovidailyrate <- 1/1300 * ifelse(tmean<0,0,tmean)^1.42
  ovidev <- cumsum(ovidailyrate)
  days <- 1:366
  min(days[ovidev >1])
}

ogd_ecl_day <- function(tmean)
{
  ovidailyrate <- 1/1300 * ifelse(tmean<0,0,tmean)^1.42
  ovidev <- cumsum(ovidailyrate)
  days <- 1:366
  oviday <- min(days[ovidev >1])
  
  tmean_mod <- c(rep(0,oviday),tmean[(oviday+1):366])
  ecldailyrate <- 1/34234 * ifelse(tmean_mod<0,0,tmean_mod)^2.27
  ecldev <- cumsum(ecldailyrate)
  days <- 1:366
  eclday <- min(days[ecldev >1])
}


result_table <- tibble(
  site = leaf_litter_temp %>% pull(site) %>% unique(),
  dd_ovi = NA,
  dd_ecl = NA,
  ogd_ovi = NA,
  ogd_ecl = NA
)
for_ind <- 0

for(whichsite in leaf_litter_temp %>% pull(site) %>% unique())
{
  for_ind <- for_ind + 1
  tmean <- leaf_litter_temp %>% filter(site == whichsite) %>% pull(tmean)
  result_table$dd_ovi[for_ind] <- dd_ovi_day(tmean, params_mean)
  result_table$dd_ecl[for_ind] <- dd_ecl_day(tmean, params_mean)
  result_table$ogd_ovi[for_ind] <- ogd_ovi_day(tmean)
  result_table$ogd_ecl[for_ind] <- ogd_ecl_day(tmean)
}

