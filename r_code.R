# Analysis Code for 
# "A mechanistic model explains variation in larval  
# tick questing phenology along an elevation gradient"
# Submitted ____????

require(tidyverse)
require(bbmle)
require(cowplot)
require(lubridate)
library(mgcv)
set.seed(31417)

# Descriptive stats in results
{
  samples <- read_csv('data/drag_sampling.csv')  
  samples %>% summarise(n(), median(larva), mean(larva), sd(larva) )
}

# Define functions for the two-peak phenology curve
# doesn't do anything on its own but these functions 
# are necessary for the next two sections to run
{
  # this one takes in parameters and data and outputs the negative log likelihood of observing those data with provided parameters
  twoPeak <- function(peak_e, tau_e, mu_e, peak_l, tau_l, mu_l, sigma_l, k, day, tickNum)
  {
    if (peak_e > 0 & tau_e > 70 & tau_e < 200 & mu_e > 0 & mu_e < 75 & peak_l > 0 & tau_l > tau_e + mu_e & tau_l < 275 & mu_l > 0 & sigma_l > 0.1 & sigma_l < 1.5)
    {
      expectedNum<- peak_e * exp(-0.5* ((day-tau_e)/mu_e)^2 ) + ifelse(day<=tau_l,0, peak_l * exp(-0.5 * (log((day-tau_l)/mu_l)/sigma_l)^2 ))
      nll <- -sum(dnbinom(x = tickNum, mu = expectedNum, size = k, log = TRUE))
    } else
    {
      nll <- 99999999
    }
    return(nll )
  }
  
  # this curve takes in the day and paramters and gives the number of larvae 
  twoPeakCurve <- function(x,
                           peak_e = coef(fit1)['peak_e'], 
                           tau_e = coef(fit1)['tau_e'], 
                           mu_e = coef(fit1)['mu_e'], 
                           peak_l = coef(fit1)['peak_l'], 
                           tau_l = coef(fit1)['tau_l'], 
                           mu_l = coef(fit1)['mu_l'], 
                           sigma_l = coef(fit1)['sigma_l'])
  {
    peak_e * exp(-0.5* ((x-tau_e)/mu_e)^2 ) + ifelse(x<=tau_l,0, peak_l * exp(-0.5 * (log((x-tau_l)/mu_l)/sigma_l)^2 ))
  }
  
  fit_phenology_nll_fun <- function(peak_e, tau_e, mu_e, peak_l, tau_l, mu_l, sigma_l, k, day, tickNum)
  {
    expectedNum<- peak_e * exp(-0.5* ((day-tau_e)/mu_e)^2 ) + ifelse(day<=tau_l,0, peak_l * exp(-0.5 * (log((day-tau_l)/mu_l)/sigma_l)^2 ))
    nll <- -sum(dnbinom(x = tickNum, mu = expectedNum, size = k, log = TRUE))
    return(nll )
  }
  
  given_phenology_nll_fun <- function(peak_mult, k, tickNum, day, phenology_pred)
  {
    pred_larva <- peak_mult * phenology_pred + 0.001
    nll <- -sum(dnbinom(x = tickNum, mu = pred_larva, size = k, log =TRUE))
    return(nll )
  }
  
  # This function is used to find confidence intervals around fit parameters
  paramRangeTwoPeak<- function(fit, param, dataList = dataList)
  {
    paramNames <- names(coef(fit))
    paramBest <- coef(fit)[param]
    startParam <- as.list(coef(fit))
    
    paramValsHigh <- c(paramBest)
    newFit <- fit
    nllsHigh <- c(-logLik(newFit) )
    while(max(nllsHigh) < (min(nllsHigh) + qchisq(p=0.95,df=1)/2) )
    {
      paramValsHigh <- c(paramValsHigh,max(paramValsHigh) + paramBest/100 )
      fixedList <- list( )
      fixedList[[param]] <- max(paramValsHigh)
      newStart <- as.list(coef(newFit))
      newFit <- mle2(twoPeak,start=newStart, data=dataList, method = 'BFGS', fixed=fixedList)
      nllsHigh <- c(nllsHigh, -logLik(newFit))
      
      
    }
    paramValsLow <- c(paramBest)
    newFit <- fit
    nllsLow <- c(-logLik(newFit) )
    while(max(nllsLow) < (min(nllsLow) + qchisq(p=0.95,df=1)/2) )
    {
      paramValsLow <- c(paramValsLow,min(paramValsLow) - paramBest/100 )
      fixedList <- list( )
      fixedList[[param]] <- min(paramValsLow)
      newStart <- as.list(coef(newFit))
      newFit <- mle2(twoPeak,start=newStart, data=dataList, method = 'BFGS', fixed=fixedList)
      nllsLow <- c(nllsLow, -logLik(newFit))
      
    }
    
    return(c(unname(paramBest),range(paramValsHigh,paramValsLow))) 
    
  }
  
}

# Define the mechanistic model and specify its parameters
# Then get mechanistic model predictions 
# Both averaged across all sites and for each site
# Input: data/leaf_litter_temp.csv
# Output: results/mech_pheno.csv
{
  # Parameters for model, See table 1
  params_mean <- list(
    ovi_m = 187,
    ovi_sd = 50,
    ecl_m = 533,
    ecl_sd = 28,
    diapause = 0.5,
    hardening = 14, 
    start_quest = 10,
    max_quest = 25,
    host_find = 0.02, 
    mort = 0.01, 
    overwinter_surv = 0.435
  )
  
  # takes in daily average temp at leaf litter and parameters
  # outputs fraction of larvae oviposited, ecolosed, questing on each day
  # See figure A1 for flow diagram
  larval_quest <- function(tmean, params)
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
    
    sum_sol <- 172 # julian day of summer solstice
    # without diapause fraction of ticks becoming active each day 
    # this takes into account hardening period
    act_prediapause <- lag(ecl, params$hardening, 0) 
    # no ticks eclosing before summer solstice enter diapause
    # fraction of them ecolsing afterwards do (Ogden et al. 2018)
    diapause_frac_daily <- c(rep(0, sum_sol + params$hardening), rep(params$diapause, tot_day - (sum_sol + params$hardening)))
    # fraction of cohort starting activity on each day
    start_act <- (1 - diapause_frac_daily) * act_prediapause
    # fraction of cohort that enters diapause after ecolosion
    diapause_frac_total <- sum(diapause_frac_daily*act_prediapause)
    
    # this section calculates the fraction of cohort questing on each day
    questing <- rep(0, tot_day)
    #here by 'winter' means when temperature drop below the min temp for larvae to quest
    days_to_winter <- min(which(tmean<params$start_quest)[which(tmean<params$start_quest) > sum_sol]) 
    active <- 0 # tracks number of active ticks
    
    # fraction of ticks active that will quest based on temp 
    quest_slope <- 1/(params$max_quest-params$start_quest)
    f_quest <- ifelse(tmean < params$max_quest, 
                      quest_slope*(tmean-params$start_quest),
                      -quest_slope*(tmean-params$max_quest) + 1 )
    f_quest <- ifelse(f_quest>1,1,ifelse(f_quest<0,0,f_quest))
    
    # now actually calculate the questing ticks each day
    for (i in 1:days_to_winter)
    {
      # add number becoming active on day i to active pool
      active <- active + start_act[i] 
      # fraction questing on day i 
      questing[i] <- active * f_quest[i]
      # remove ticks from active cohort if they find a host or die
      active <- (1 - params$host_find * f_quest[i] - params$mort) * active
    }
    # fraction of the cohort which survives overwinters
    # overwinter survival times three groups: still active at end of questing season; those starting activity after questing season; those that entered diapause earlier
    larvae_overwinter <- params$overwinter_surv * (active + sum(start_act[(days_to_winter+1):tot_day]) + diapause_frac_total)
    
    # add overwintered larvae to questing vector
    # these are early-summer questing ticks
    for (i in 1:days_to_winter)
    {
      questing[i] <- questing[i] + larvae_overwinter * f_quest[i]
      larvae_overwinter <- (1 - params$host_find * f_quest[i] - params$mort) * larvae_overwinter
    }
    
    output_tibble <- tibble(
      ovi = ovi, # fraction of cohort oviposited each day
      ecl = ecl, # fraction of cohort eclosed each day
      qst = questing, # fraction of cohort questing each day (sums > 1 since most ticks quest on more than one day)
      qst_norm = questing/sum(questing) # fraction of questing larvae on a given day, for comparison to real data
    )
    
    return(output_tibble)
  }
  
  samples <- read_csv('data/drag_sampling.csv')  
  leaf_litter_temp <- read_csv('data/leaf_litter_temp.csv')
  
  mech_model_pred <- tibble(
    julian = numeric(),
    larva_frac = numeric(),
    site = character())
  
  #filter two sites that never had any larvae
  mysites <- samples %>% 
    filter(site != 'Snowbowl', site != 'Crystal') %>%
    pull(site) %>% 
    unique()  
  
  for (which_site in c(mysites, 'mean'))
  {
    if (which_site == 'mean') {
      elev_cat_temp <- leaf_litter_temp %>%
        group_by(jday) %>%
        summarise(tmean = mean(tmean))
    } else {
      elev_cat_temp <- leaf_litter_temp %>%
        filter(site == which_site) %>%
        group_by(jday) %>%
        summarise(tmean = mean(tmean)) 
    }
    temp_pred <- tibble(
      julian = 1:365,
      larva_frac = larval_quest(elev_cat_temp$tmean[1:365], params_mean)$qst_norm,
      site = which_site)
    mech_model_pred<- rbind(mech_model_pred,temp_pred)
  }
  
  write_csv(mech_model_pred, file = 'results/mech_pheno.csv')
}

# Compare the three models
# Inputs: data/drag_sampling.csv, results/mech_pheno.csv
# Output: results/model_fits.RData, results/full_model_pred.csv
{
  samples <- read_csv('data/drag_sampling.csv')  
  
  ## Parameter starting points for mle2 search
  {
    param_guess_ph <- list()
    param_guess_ph[['Foote']] <- list(peak_e=30, tau_e=160, mu_e=20, peak_l=70, tau_l=200, mu_l=20, sigma_l=0.1,k=0.2)
    param_guess_ph[['Chipman']] <- list(peak_e=30, tau_e=160, mu_e=20, peak_l=70, tau_l=200, mu_l=20, sigma_l=0.1,k=0.2)
    param_guess_ph[['Major']] <- list(peak_e=30, tau_e=160, mu_e=20, peak_l=70, tau_l=200, mu_l=20, sigma_l=0.1,k=0.2)
    param_guess_ph[['Lourie']] <- list(peak_e=30, tau_e=160, mu_e=20, peak_l=70, tau_l=200, mu_l=20, sigma_l=0.1,k=0.2)
    
    param_guess_ph[['Chipman2']] <- list(peak_e=25, tau_e=135, mu_e=35, peak_l=5, tau_l=200, mu_l=60, sigma_l=0.1,k=0.2)
    param_guess_ph[['Gorge']] <- list(peak_e=25, tau_e=135, mu_e=35, peak_l=5, tau_l=200, mu_l=60, sigma_l=0.1,k=0.2)
    param_guess_ph[['BRF']] <- list(peak_e=25, tau_e=135, mu_e=35, peak_l=5, tau_l=200, mu_l=60, sigma_l=0.1,k=0.2)
    
    param_guess_ph[['BRF2']] <- list(peak_e=7, tau_e=170, mu_e=11.5, peak_l=0.03, tau_l=203, mu_l=50, sigma_l=1.5,k=0.03) 
    param_guess_ph[['Frost']] <- list(peak_e=7, tau_e=170, mu_e=11.5, peak_l=0.03, tau_l=203, mu_l=50, sigma_l=1.5,k=0.03) 
    param_guess_ph[['SPIN']] <- list(peak_e=7, tau_e=170, mu_e=11.5, peak_l=0.03, tau_l=203, mu_l=50, sigma_l=1.5,k=0.03) 
    param_guess_ph[['Gilmore']] <- list(peak_e=7, tau_e=170, mu_e=11.5, peak_l=0.03, tau_l=203, mu_l=50, sigma_l=1.5,k=0.03) 
    
    param_guess_ph[['mean']] <- list(peak_e = 10, tau_e = 135, mu_e = 20, peak_l = 30, tau_l = 200, mu_l = 50, sigma_l =0.2, k = 0.1)
    
    param_guess_me <- list()
    param_guess_me[['Foote']] <- list(peak_mult=5000, k = 0.1)
    param_guess_me[['Chipman']] <- list(peak_mult=5000, k = 0.1)
    param_guess_me[['Major']] <- list(peak_mult=5000, k = 0.1)
    param_guess_me[['Lourie']] <- list(peak_mult=5000, k = 0.1)
    
    param_guess_me[['Chipman2']] <- list(peak_mult=750, k = 0.08)
    param_guess_me[['Gorge']] <- list(peak_mult=750, k = 0.08)
    param_guess_me[['BRF']] <- list(peak_mult=750, k = 0.08)
    
    param_guess_me[['BRF2']] <- list(peak_mult=50, k = 0.02)
    param_guess_me[['Frost']] <- list(peak_mult=50, k = 0.02) 
    param_guess_me[['SPIN']] <- list(peak_mult=50, k = 0.02)
    param_guess_me[['Gilmore']] <- list(peak_mult=50, k = 0.02) 
    
    make_parscale <- function(x)
    {
      my_parscale <- x
      for (i in 1:length(x)) my_parscale[[i]] <- x[[i]]/x[[1]]
      unlist(my_parscale)
    }
  }
  
  ## fit a single phenomological phenology across all sites 
  {
    data_list <- samples %>% 
      select(day = julian, tickNum = larva) %>%
      as.list()
    
    my_parscale <- make_parscale(param_guess_ph[['mean']])
    fit1 <- mle2(fit_phenology_nll_fun, start=param_guess_ph[['mean']], data=data_list, method='BFGS', control = list(parscale = my_parscale))
    
    phenomological_mean_pred <- tibble(
      julian = 1:365,
      pred_larva = twoPeakCurve(1:365)/sum(twoPeakCurve(1:365))
    )
  }
  
  # now compare the four models
  {
    all_mod_pred <- tibble(
      julian = numeric(),
      larva = numeric(),
      site = character(),
      mod_type = character()
    )
    
    phenomological_site_nll <- 0
    phenomological_mean_nll <- 0
    mechanistic_site_nll <- 0
    mechanistic_mean_nll <- 0
    
    fitList<-list()
    mech_model_pred <- read_csv('results/mech_pheno.csv')
  
    # loop through the sites and fit each model for each site
    # the mean-level models do not fit a phenology at each site
    # but just the height of the peaks
    for (which_site in mysites)
    {
      
      # site-level phenomological fits
      {
        data_list <- samples %>% 
          filter(site == which_site) %>%
          select(day = julian, tickNum = larva) %>%
          as.list()
        
        my_parscale <- make_parscale(param_guess_ph[[which_site]])
        fit1 <- mle2(twoPeak, start=param_guess_ph[[which_site]], data=data_list,method='BFGS')
        fitList[['phenom_site']][[which_site]] <- fit1
        
        phenomological_site_nll <- phenomological_site_nll - as.numeric(logLik(fit1))
        temp_pred <- tibble(julian = 1:365, larva = twoPeakCurve(1:365), site = which_site, mod_type = 'Site phenomological')
        all_mod_pred <- bind_rows(all_mod_pred, temp_pred)
      }
      
      # mean-level phenomological fits
      {
        data_list <- samples %>%
          filter(site == which_site) %>%
          left_join(phenomological_mean_pred, by = 'julian') %>%
          select(day = julian, tickNum = larva, phenology_pred = pred_larva) %>%
          as.list()
        
        my_parscale <- make_parscale(param_guess_me[[which_site]])
        fit1 <- mle2(given_phenology_nll_fun, start=param_guess_me[[which_site]], data=data_list, method='BFGS', control = list(parscale = my_parscale) )    
        fitList[['phenom_mean']][[which_site]] <- fit1
        
        phenomological_mean_nll <- phenomological_mean_nll - as.numeric(logLik(fit1))  
        temp_pred <- tibble(julian = 1:365, larva = unname(coef(fit1)['peak_mult']) * phenomological_mean_pred$pred_larva, site = which_site, mod_type = 'Mean phenomological')
        all_mod_pred <- bind_rows(all_mod_pred, temp_pred)
        
      }
      
      # site-level mechanistic fits
      {
        data_list <- samples %>%
          filter(site == which_site) %>%
          left_join(mech_model_pred, by = c('site', 'julian')) %>%
          select(day = julian, tickNum = larva, phenology_pred = larva_frac) %>%
          as.list()
        
        my_parscale <- make_parscale(param_guess_me[[which_site]])
        fit1 <- mle2(given_phenology_nll_fun, start=param_guess_me[[which_site]], data=data_list, method='BFGS', control = list(parscale = my_parscale) )    
        fitList[['mech_site']][[which_site]] <- fit1
        
        mechanistic_site_nll <- mechanistic_site_nll - as.numeric(logLik(fit1))    

        this_site_pheno <- mech_model_pred %>%
          filter(site == which_site) %>%
          pull(larva_frac)
        
        temp_pred <- tibble(julian = 1:365, larva = unname(coef(fit1)['peak_mult']) * this_site_pheno, site = which_site, mod_type = 'Site mechanistic')
        all_mod_pred <- bind_rows(all_mod_pred, temp_pred)
      }
      
      # mean-level mechanistic fits
      {
        mean_mech_model_pred <- mech_model_pred %>% filter(site == 'mean')
        
        data_list <- samples %>%
          filter(site == which_site) %>%
          left_join(mean_mech_model_pred, by =  'julian') %>%
          select(day = julian, tickNum = larva, phenology_pred = larva_frac) %>%
          as.list()
        
        fit1 <- mle2(given_phenology_nll_fun, 
                     start=param_guess_me[[which_site]], 
                     data=data_list, method='BFGS', 
                     control = list(parscale = c(peak_mult = param_guess_me[[which_site]]$peak_mult/param_guess_me[[which_site]]$k, k = 1) ) )
        fitList[['mech_mean']][[which_site]] <- fit1
        
        mechanistic_mean_nll <- mechanistic_mean_nll - as.numeric(logLik(fit1))    

        mean_pheno_just_nums <- mean_mech_model_pred %>%
          pull(larva_frac)
        
        temp_pred <- tibble(julian = 1:365, larva = unname(coef(fit1)['peak_mult']) * mean_pheno_just_nums, site = which_site, mod_type = 'Mean mechanistic')
        all_mod_pred <- bind_rows(all_mod_pred, temp_pred)
      }
      
      #smoothed 
      {
        smoothed_pheno <- samples %>%
          filter(site == which_site) %>%
          gam(larva~ s(julian, sp = 0.005), family = nb(), data = .)
        temp_pred <- tibble(
          julian = 1:365, 
          larva = predict(smoothed_pheno, newdata = data.frame(julian = 1:365), type = 'response'),
          site = which_site, 
          mod_type = 'Smoothed'
        )
        all_mod_pred <- bind_rows(all_mod_pred, temp_pred)
      }
      
      
    }  
  
  }
  
  
  phenomological_site_aic <- phenomological_site_nll + 2*8*length(mysites)
  phenomological_mean_aic <- phenomological_mean_nll + 2*6 + 2*2*length(mysites)
  mechanistic_site_aic <- mechanistic_site_nll + 2*2*length(mysites)
  mechanistic_mean_aic <- mechanistic_mean_nll + 2*2*length(mysites)
  
  save(fitList, file = "results/model_fits.RData")
  write_csv(x = all_mod_pred, file = 'results/full_model_pred.csv')
}


all_mod_pred %>% 
  mutate(mean_or_not = ifelse(str_detect(mod_type, "Mean"),'Y', 'N')) %>%
  ggplot(aes(julian, larva, color = mod_type, lty = mean_or_not)) + 
  geom_line() + 
  facet_wrap(~site, scales = 'free_y') +
  scale_color_manual(values = c("Site phenomological" = 'red', 
                                "Mean phenomological" = 'red',
                                'Site mechanistic' = 'blue',
                                'Mean mechanistic' = 'blue',
                                'Smoothed' = 'black')) +
  xlim(c(100,300))

#### # HERE ### # ## # # # # ###
#### # HERE ### # ## # # # # ###
#### # HERE ### # ## # # # # ###
#### # HERE ### # ## # # # # ###
#### # HERE ### # ## # # # # ###
#### # HERE ### # ## # # # # ###

# New analysis to do model comparison
{
  samples <- read_csv('data/drag_sampling.csv') %>% 
    mutate(elevCat = cut(elev,c(0,200,400,1000),c('low','mid','high'))) 

  
  # Run model on mean param
  {
    
    params_mean <- list(
      ovi_m = 187,
      ovi_sd = 50,
      ecl_m = 533,
      ecl_sd = 28,
      diapause = 0.5,
      hardening = 14, #### change here based on leal et al citaiton of 1-7d hardening
      start_quest = 10,
      max_quest = 25,
      host_find = 0.02, 
      mort = 0.0099, 
      overwinter_surv = 0.45
    )
    leaf_litter_temp <- read_csv(file = 'data/leaf_litter_temp.csv')
  
    mech_model_pred <- tibble(
      julian = numeric(),
      larva_frac = numeric(),
      elevCat = character())
  
    for (which_elev in c('low', 'mid', 'high', 'mean'))
      {
        if (which_elev == 'mean') {
          elev_cat_temp <- leaf_litter_temp %>%
            group_by(jday) %>%
            summarise(tmean = mean(tmean))
        } else {
          elev_cat_temp <- leaf_litter_temp %>%
            filter(elevCat == which_elev) %>%
            group_by(jday) %>%
            summarise(tmean = mean(tmean)) 
          }
        temp_pred <- tibble(
          julian = 1:365,
          larva_frac = larval_quest(elev_cat_temp$tmean[1:365], params_mean)$qst_norm,
          elevCat = which_elev)
        mech_model_pred<- rbind(mech_model_pred,temp_pred)
      }
    }
  

  mechanistic_nll_fun <- function(low_mult, mid_mult, high_mult, k_low, k_mid, k_high, tickNum, day, phenology_pred)
  {
    pred_larva <- ifelse(elevCat == 'low',low_mult, ifelse(elevCat == 'mid', mid_mult, high_mult)) * phenology_pred + 0.001
    k <- ifelse(elevCat == 'low', k_low, ifelse(elevCat == 'mid', k_mid, k_high))
    nll <- -sum(dnbinom(x = tickNum, mu = pred_larva, size = k, log =TRUE))
    return(nll )
  }
  
  phenomological_nll_fun <- function(mid_mult, high_mult, peak_e, tau_e, mu_e, peak_l, tau_l, mu_l, sigma_l, k_low, k_mid, k_high, day, tickNum, elevCat)
  {
    peak_e <- ifelse(elevCat == 'low', peak_e, ifelse(elevCat == 'mid', peak_e*mid_mult, peak_e*high_mult))
    peak_l <- ifelse(elevCat == 'low', peak_l, ifelse(elevCat == 'mid', peak_l*mid_mult, peak_l*high_mult))
    k <- ifelse(elevCat == 'low', k_low, ifelse(elevCat == 'mid', k_mid, k_high))
    if (peak_e > 0 & tau_e > 70 & tau_e < 200 & mu_e > 0 & mu_e < 75 & peak_l > 0 & tau_l > tau_e + mu_e & tau_l < 275 & mu_l > 0 & sigma_l > 0.1 & sigma_l < 1.5)
    {
      expectedNum<- peak_e * exp(-0.5* ((day-tau_e)/mu_e)^2 ) + ifelse(day<=tau_l,0, peak_l * exp(-0.5 * (log((day-tau_l)/mu_l)/sigma_l)^2 ))
      nll <- -sum(dnbinom(x = tickNum, mu = expectedNum, size = k, log = TRUE))
    } else
    {
      nll <- 99999999
    }
    return(nll )
  }
  
  data_list_mech_full <- samples %>%
    left_join(mech_model_pred, by = c('elevCat', 'julian')) %>%
    select(day = julian, tickNum = larva, elevCat, phenology_pred = larva_frac) %>%
    as.list()
  data_list_mech_mean<- samples %>%
    left_join(mech_model_pred %>% filter(elevCat == 'mean'), by = 'julian') %>%
    select(day = julian, tickNum = larva, elevCat = elevCat.x, phenology_pred = larva_frac) %>%
    as.list()
  data_list_phenomological <- samples %>% 
    select(day = julian, tickNum = larva, elevCat = elevCat) %>%
    as.list()
  
  guess_param_mech <- list(low_mult = 9500, mid_mult = 2000, high_mult = 50, k_low = 0.11, k_mid = 0.095, k_high = 0.02)
  guess_param_phenomological <- list(peak_e=35, mid_mult=0.54, high_mult=0.05, tau_e=160, mu_e=20, peak_l=100, tau_l=180, mu_l=50, sigma_l=0.2,k_low=0.2,k_mid=0.2,k_high=0.02)
  
  mechanistic_full_fit <- mle2(mechanistic_nll_fun, start=guess_param_mech, data=data_list_mech_full) 
  mechanistic_mean_fit <- mle2(mechanistic_nll_fun, start=guess_param_mech, data=data_list_mech_mean) 
  phenomological_fit <- mle2(phenomological_nll_fun, start = guess_param_phenomological, data = data_list_phenomological)
  
  
  AIC_full <- -as.numeric(logLik(mechanistic_full_fit)) + 2*6 + (2*6^2 + 2*6)/(length(data_list_mech_full$day) - 6 - 1)
  AIC_mean <- -as.numeric(logLik(mechanistic_mean_fit)) + 2*6 + (2*6^2 + 2*6)/(length(data_list_mech_mean$day) - 6 - 1)
  AIC_phenomological <- -as.numeric(logLik(phenomological_fit)) + 2*12 + (2*12^2 + 2*12)/(length(data_list_phenomological$day) - 12 - 1)
}


# get model predictions for graphing
{
  
    
  mean_pred <- mech_model_pred %>%
    filter(elevCat == 'mean') %>%
    pull(larva_frac)
  
  low_pred <- mech_model_pred %>%
    filter(elevCat == 'low') %>%
    pull(larva_frac)
  
  mid_pred <- mech_model_pred %>%
    filter(elevCat == 'mid') %>%
    pull(larva_frac)
  
  high_pred <- mech_model_pred %>%
    filter(elevCat == 'high') %>%
    pull(larva_frac)
  
 
 fit1 <- phenomological_fit   
    
  predicted_larvae <- 
    tibble(julian = rep(1:365,9),
           elevCat = rep(c(rep('low',365),rep('mid',365),rep('high',365)),3),
           larva = c(coef(mechanistic_mean_fit)['low_mult']*mean_pred,
                     coef(mechanistic_mean_fit)['mid_mult']*mean_pred,
                     coef(mechanistic_mean_fit)['high_mult']*mean_pred,
                     coef(mechanistic_full_fit)['low_mult']*low_pred,
                     coef(mechanistic_full_fit)['mid_mult']*mid_pred,
                     coef(mechanistic_full_fit)['high_mult']*high_pred,
                     twoPeakCurve(1:365),
                     coef(phenomological_fit)['mid_mult']*twoPeakCurve(1:365),
                     coef(phenomological_fit)['high_mult']*twoPeakCurve(1:365)),
           mod_type = c(rep('mechanistic_mean',365*3),rep('mechanistic_full',365*3),rep('phenomological',365*3))
    ) %>%
    mutate(elevCat = factor(elevCat, levels = c('low', 'mid', 'high')))
    
    samples_expanded <- tibble(
        julian = rep(samples$julian,3),
        larva = rep(samples$larva,3),
        elevCat = rep(samples$elevCat,3),
        mod_type = c(rep('mechanistic_mean',length(samples$julian)),rep('mechanistic_full',length(samples$julian)),rep('phenomological',length(samples$julian)))
      ) %>%
      mutate(elevCat = factor(elevCat, levels = c('low', 'mid', 'high')))
    
    
    # summaries for predictions
    mid_point <- predicted_larvae %>%
      group_by(elevCat, mod_type) %>%
      mutate(delta = larva - lag(larva, default = 0),
             mid_point = ifelse(delta > 0 & lag(delta) < 0, 1, 0) ) %>%
      filter(mid_point == 1) %>%
      select(-c(larva, mid_point, delta)) %>%
      rename(day_mid = julian)
    
    predicted_larvae %>%
      left_join(mid_point, by = c('mod_type', 'elevCat')) %>%
      mutate(early_larva = ifelse(julian< day_mid, larva, 0),
             late_larva = ifelse(julian< day_mid, 0, larva) ) %>%
      group_by(mod_type, elevCat) %>%
      summarise(frac_early = sum(early_larva)/sum(early_larva+late_larva))
  
    
    library(mgcv)
    
    smoothed_obs_low_gam <- gam(larva ~ s(julian, sp = 0.001), data = samples %>% filter(elevCat == 'low'), family = negbin(theta = 0.1))
    smoothed_obs_mid_gam <- gam(larva ~ s(julian, sp = 0.001), data = samples %>% filter(elevCat == 'mid'), family = negbin(theta = 0.1))
    smoothed_obs_high_gam <- gam(larva ~ s(julian, sp = 0.001), data = samples %>% filter(elevCat == 'high'), family = negbin(theta = 0.1))
  
    smoothed_obs_low <- predict(smoothed_obs_low_gam, newdata = data.frame(julian = 1:365), type = 'response')
    smoothed_obs_mid <- predict(smoothed_obs_mid_gam, newdata = data.frame(julian = 1:365), type = 'response')
    smoothed_obs_high <- predict(smoothed_obs_high_gam, newdata = data.frame(julian = 1:365), type = 'response')
    
    smoothed_larvae <- tibble(julian = rep(1:365,3), 
                              larva = c(smoothed_obs_low,smoothed_obs_mid, smoothed_obs_high),
                              elevCat = c(rep('low',365),rep('mid',365),rep('high',365)) )
    
    mid_point2<- smoothed_larvae %>%
      group_by(elevCat) %>%
      mutate(delta = larva - lag(larva, default = 0),
             mid_point = ifelse(delta > 0 & lag(delta) < 0, 1, 0) ) %>%
      filter(mid_point == 1) %>%
      select(-c(larva, mid_point, delta)) %>%
      rename(day_mid = julian)
    
    smoothed_larvae %>%
      left_join(mid_point2, by = 'elevCat') %>%
      mutate(early_larva = ifelse(julian< day_mid, larva, 0),
             late_larva = ifelse(julian< day_mid, 0, larva) ) %>%
      group_by(elevCat) %>%
      summarise(frac_early = sum(early_larva)/sum(early_larva+late_larva))
    
    
    smooth_expanded <-tibble(
      julian = rep(smoothed_larvae$julian,3),
      larva = rep(smoothed_larvae$larva,3),
      elevCat = rep(smoothed_larvae$elevCat,3),
      mod_type = c(rep('mechanistic_mean',length(smoothed_larvae$julian)),rep('mechanistic_full',length(smoothed_larvae$julian)),rep('phenomological',length(smoothed_larvae$julian)))
    ) %>%
      mutate(elevCat = factor(elevCat, levels = c('low', 'mid', 'high')))   
    
    smooth_expanded%>%
      ggplot(aes(julian,larva)) +
      geom_path() +
      #geom_smooth(data = samples_expanded,se = F) +
      geom_path(data= predicted_larvae, lty = 2) +
      facet_grid(elevCat ~ mod_type, scales = 'free_y') +
      coord_cartesian(xlim=c(100,300))
    
    
}

## Make a graph of the pheno model with parameters

peak_e_example <- 30
tau_e_example <- 140
mu_e_example <- 20
peak_l_example <- 60
tau_l_example <- 210
mu_l_example <- 20
sigma_l_example <- 0.5

pheno_example_data <- tibble(
  day = 50:300,
  larvae = peak_e_example * exp(-0.5* ((day-tau_e_example)/mu_e_example)^2 ) + ifelse(day<=tau_l_example,0, peak_l_example * exp(-0.5 * (log((day-tau_l_example)/mu_l_example)/sigma_l_example)^2 ))
)

pdf('figures/phenomenological_example.pdf',width=7.25,height=3)
  
  pheno_example_data %>%
    ggplot(aes(day,larvae)) +
    geom_line() +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank()) +
    xlab('Julian day') +
    ylab('Questing larvae') +
    annotate('segment', 
             x=50, xend = tau_e_example,
             y=peak_e_example+1, yend=peak_e_example+1, 
             arrow = arrow(ends='both', length = unit(0.1,'cm'))) +
    annotate('text',x=(50+tau_e_example)/2,y=peak_e_example+3,label=expression(tau['e'])) +
    annotate('segment',
             x = tau_e_example, xend = tau_e_example, 
             y = 0, yend = peak_e_example,
             arrow = arrow(ends='both', length = unit(0.1,'cm'))) +
    annotate('text',x=tau_e_example+4,y=peak_e_example/2,label=expression('H'['e'])) +
    annotate('text',x=tau_e_example,y=peak_e_example+6,label=expression(paste("Shape parameter ",sigma['e']))) +
    annotate('segment', 
             x=50, xend = tau_l_example,
             y=peak_l_example+1, yend=peak_l_example+1, 
             arrow = arrow(ends='both', length = unit(0.1,'cm'))) +
    annotate('text',x=(50+tau_l_example)/2,y=peak_l_example+3,label=expression(tau['l'])) +
    annotate('segment', 
             x=tau_l_example, xend = tau_l_example+mu_l_example,
             y=peak_l_example+1, yend=peak_l_example+1, 
             arrow = arrow(ends='both', length = unit(0.1,'cm'))) +
    annotate('text',x=(2*tau_l_example+mu_l_example)/2,y=peak_l_example+4,label=expression(mu['l'])) +
    annotate('segment',
             x = tau_l_example+mu_l_example, xend = tau_l_example+mu_l_example, 
             y = 0, yend = peak_l_example,
             arrow = arrow(ends='both', length = unit(0.1,'cm'))) +
    annotate('text',x=tau_l_example+mu_l_example+3,y=peak_l_example/2,label=expression('H'['l'])) +
    annotate('text',x=tau_l_example+mu_l_example+37,y=peak_l_example*0.75,label=expression(paste("Shape parameter ",sigma['l'])))
    
    
dev.off()

  

  



