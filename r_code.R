# Analysis Code for 
# "A mechanistic model explains variation in larval  
# tick questing phenology along an elevation gradient"
# Submitted Sept 2024, re-submitted Jan 2025

require(tidyverse)
require(bbmle)
require(cowplot)
library(mgcv)

# Descriptive stats in results
{
  samples <- read_csv('data/drag_sampling.csv')  
  samples %>% summarise(n(), median(larva), mean(larva), sd(larva) )
}

# Define the mechanistic model and specify its parameters
# Then get mechanistic model predictions 
# Both averaged across all sites and for each site
# Input: data/leaf_litter_temp.csv
# Output: results/mech_pheno.csv
{
  # Parameters for model, See table 1
  params_mean <- list(
    adult_start_quest = 3,
    adult_max_quest = 8,
    adult_host_find = 0.04,
    adult_mort = 0.006,
    ovi_m = 188,
    ovi_sd = 50.1,
    ecl_m = 532.1,
    ecl_sd = 38.5,
    diapause = 0.5,
    hardening = 7, 
    start_quest = 10,
    max_quest = 25,
    host_find = 0.02, 
    mort = 0.01, 
    overwinter_surv = 0.45
  )
  
  # takes in daily average temp at leaf litter and parameters
  # outputs fraction of larvae oviposited, ecolosed, questing on each day
  # See figure A1 for flow diagram
  larval_quest <- function(tmean, params)
  {
    tot_day <- length(tmean)
    
    # fraction of adults that are engroged on each day
    adult_questing_start_day <- 213 # Aug 1 first day possible (but will be too hot)
    adult_questing_end_day <- 182 # july 1
    adult_days <- c(adult_questing_start_day:tot_day,1:tot_day) # adults start questing on Oct 1 and then continue to the next year
    adult_engorged_long <- numeric(length(adult_days))
  
    adult_active <- 1 # start the total fractino of adults
    
    # fraction of ticks active that will quest based on temp 
    adult_quest_slope <- 1/(params$adult_max_quest-params$adult_start_quest)
    adult_f_quest <- ifelse(tmean < params$adult_max_quest, 
                          adult_quest_slope*(tmean-params$adult_start_quest),
                      -adult_quest_slope*(tmean-params$adult_max_quest) + 1 )
    adult_f_quest <- ifelse(adult_f_quest>1,1,ifelse(adult_f_quest<0,0,adult_f_quest))
    
    # now actually calculate the questing ticks each day
    dec31_ind <- tot_day - adult_questing_start_day + 1
    for (i in 1:length(adult_days))
    {
      adult_engorged_long[i] <- adult_active * adult_f_quest[adult_days[i]] * params$adult_host_find
      adult_active <- adult_active - (adult_active * adult_f_quest[adult_days[i]] * params$adult_host_find) - ifelse(adult_f_quest[adult_days[i]]>0, adult_active * adult_f_quest[adult_days[i]] * params$adult_mort,0)
      if (i == (dec31_ind + adult_questing_end_day)) adult_active<-0
    }

    adult_engorged <- c(sum(adult_engorged_long[1:dec31_ind]),adult_engorged_long[(dec31_ind+2):length(adult_engorged_long)])
    adult_engorged <- adult_engorged/sum(adult_engorged)
    
    adult_engorged <- lag(adult_engorged,7)
    adult_engorged <- ifelse(is.na(adult_engorged),0,adult_engorged)

    # fraction of cohort oviposited on each day
    ovi <- numeric(tot_day)
    for (i in 1:(tot_day-1))
    {
      i1 <- i+1
      dd6 <- c(rep(0,i),cumsum(ifelse(tmean[i1:tot_day]>6,tmean[i1:tot_day]-6,0)))
      # ovi_prob is prob of oviposition on each day if engorged is on day i
      ovi_prob <- pnorm(dd6,params$ovi_m,params$ovi_sd) - pnorm(lag(dd6,n=1,default = 0),params$ovi_m,params$ovi_sd)
      ovi <- ovi + adult_engorged[i] * ovi_prob
    }
   

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
    # plus all the ticks that didn't finish development and are kicked to next year
    diapause_frac_total <- sum(diapause_frac_daily*act_prediapause) + (1-sum(ecl))
    
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
    filter(site != 'Crystal') %>%
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

# Define functions for the phenomenological model
# and other functions for maximum likelihood fitting
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
  
  # this curve takes in the day and parameters and gives the number of larvae 
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
  
}

## Make a figure showing the example phenomenological model with parameters
## Output: figures/phenomenological_example.pdf
{
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
  
  
  pdf('/Users/dallen@middlebury.edu/My Drive/Research/lyme/2024/pheno_manu/figures/phenomenological_example.pdf',width=7.25,height=3)
  
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
    #annotate('segment',
    #         x = tau_e_example, xend = tau_e_example, 
    #         y = 0, yend = peak_e_example,
    #         arrow = arrow(ends='both', length = unit(0.1,'cm'))) +
    #annotate('text',x=tau_e_example+4,y=peak_e_example/2,label=expression('H'['e'])) +
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
}

# Compare the four models
# Inputs: data/drag_sampling.csv, results/mech_pheno.csv
# Outputs: results/model_fits.RData, results/full_model_pred.csv
{
  samples <- read_csv('data/drag_sampling.csv')  
  
  ## Parameter starting points for mle2 search
  {
    param_guess_ph <- list()
    param_guess_ph[['Foote']] <- list(peak_e=30, tau_e=160, mu_e=20, peak_l=200, tau_l=200, mu_l=20, sigma_l=0.1,k=0.2)
    param_guess_ph[['Chipman']] <- list(peak_e=30, tau_e=160, mu_e=20, peak_l=70, tau_l=200, mu_l=20, sigma_l=0.1,k=0.2)
    param_guess_ph[['Major']] <- list(peak_e=30, tau_e=160, mu_e=20, peak_l=70, tau_l=200, mu_l=20, sigma_l=0.1,k=0.2)
    param_guess_ph[['Lourie']] <- list(peak_e=25, tau_e=160, mu_e=20, peak_l=50, tau_l=220, mu_l=20, sigma_l=0.1,k=0.2)
    param_guess_ph[['Jackson']] <- list(peak_e=30, tau_e=160, mu_e=20, peak_l=200, tau_l=200, mu_l=20, sigma_l=0.1,k=0.2)
    
        
    param_guess_ph[['UpperChipman']] <- list(peak_e=25, tau_e=135, mu_e=35, peak_l=5, tau_l=200, mu_l=60, sigma_l=0.1,k=0.2)
    param_guess_ph[['Gorge']] <- list(peak_e=25, tau_e=135, mu_e=35, peak_l=5, tau_l=200, mu_l=60, sigma_l=0.1,k=0.2)
    param_guess_ph[['BRF']] <- list(peak_e=5, tau_e=150, mu_e=35, peak_l=5, tau_l=200, mu_l=50, sigma_l=0.1,k=0.2)
    
    param_guess_ph[['BRF2']] <- list(peak_e=2, tau_e=170, mu_e=11.5, peak_l=0.03, tau_l=203, mu_l=50, sigma_l=1.5,k=0.03) 
    param_guess_ph[['Frost']] <- list(peak_e=2, tau_e=170, mu_e=11.5, peak_l=0.03, tau_l=203, mu_l=50, sigma_l=1.5,k=0.03) 
    param_guess_ph[['SPIN']] <- list(peak_e=2, tau_e=170, mu_e=20, peak_l=0.03, tau_l=203, mu_l=50, sigma_l=1.5,k=0.03) 
    param_guess_ph[['Gilmore']] <- list(peak_e=2, tau_e=170, mu_e=11.5, peak_l=0.03, tau_l=203, mu_l=50, sigma_l=1.5,k=0.03) 
    param_guess_ph[['Snowbowl']] <- list(peak_e=0.01, tau_e=170, mu_e=11.5, peak_l=0.0001, tau_l=203, mu_l=50, sigma_l=1.5,k=0.03) 
    
        
    param_guess_ph[['mean']] <- list(peak_e = 10, tau_e = 150, mu_e = 15, peak_l = 27, tau_l = 200, mu_l = 50, sigma_l =0.2, k = 0.1)
    
    param_guess_me <- list()
    param_guess_me[['Foote']] <- list(peak_mult=5000, k = 0.1)
    param_guess_me[['Chipman']] <- list(peak_mult=5000, k = 0.1)
    param_guess_me[['Major']] <- list(peak_mult=5000, k = 0.1)
    param_guess_me[['Lourie']] <- list(peak_mult=5000, k = 0.1)
    param_guess_me[['Jackson']] <- list(peak_mult=5000, k = 0.1)
    
    param_guess_me[['UpperChipman']] <- list(peak_mult=750, k = 0.08)
    param_guess_me[['Gorge']] <- list(peak_mult=750, k = 0.08)
    param_guess_me[['BRF']] <- list(peak_mult=750, k = 0.08)
    
    param_guess_me[['BRF2']] <- list(peak_mult=50, k = 0.02)
    param_guess_me[['Frost']] <- list(peak_mult=50, k = 0.02) 
    param_guess_me[['SPIN']] <- list(peak_mult=50, k = 0.02)
    param_guess_me[['Gilmore']] <- list(peak_mult=50, k = 0.02) 
    param_guess_me[['Snowbowl']] <- list(peak_mult=10, k = 0.02) 
    
        
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
          rbind(data.frame(
                site = NA,
                elev = NA,
                date = NA,
                julian = c(rep(99,5),rep(75,5)),
                larva = rep(0,10)
          )) %>%
          gam(larva~ s(julian, sp = 0.001), family = nb(), data = .)
        temp_pred <- tibble(
          julian = 1:365, 
          larva = as.vector(predict(smoothed_pheno, newdata = data.frame(julian = 1:365), type = 'response')),
          site = which_site, 
          mod_type = 'Smoothed'
        )
        all_mod_pred <- bind_rows(all_mod_pred, temp_pred)
      }
      
      
    }  
  
  }
  
  
  phenomological_site_n_param <- 8*length(mysites) ## seven for phenom model plus dispersion (=8) for each site
  phenomological_mean_n_param <- 6 + 2*length(mysites) ## seven for single phenom model plus max and dispersion for each site
  mechanistic_site_n_param <- 2*length(mysites) ## max and dispersion for each site
  mechanistic_mean_n_param <- 2*length(mysites) ## max and dispersion for each site

  phenomological_site_aic <- phenomological_site_nll + 2*phenomological_site_n_param
  phenomological_mean_aic <- phenomological_mean_nll + 2*phenomological_mean_n_param
  mechanistic_site_aic <- mechanistic_site_nll + 2*mechanistic_site_n_param
  mechanistic_mean_aic <- mechanistic_mean_nll + 2*mechanistic_mean_n_param
  
  save(fitList, file = "results/model_fits.RData")
  write_csv(x = all_mod_pred, file = 'results/full_model_pred.csv')
}

## Make first results figure
## Inputs: data/drag_sampling.csv and results/full_model_pred.csv
## Output: figures/pheno_bysite.pdf
{
  ## relabel sites by elevation rather than site name
  samples <- read_csv('data/drag_sampling.csv')  
  siteelev <- samples %>%
    group_by(site) %>%
    summarise(elev = mean(elev)) %>% 
    ungroup()
  all_mod_pred <- read_csv('results/full_model_pred.csv') %>%
    left_join(siteelev,by='site') %>%
    mutate(elevlabel = paste(elev,' m'))

  yaxis <- expression(paste("Larvae (per 200 ", m^2, " sample)"))
  p1 <- all_mod_pred %>% 
    filter(julian >100, julian < 300) %>%
    filter(mod_type %in% c('Site mechanistic', 'Smoothed')) %>%
    ggplot(aes(julian, larva, linetype = mod_type)) + 
    geom_line() + 
    facet_wrap(~elevlabel, scales = 'free_y') +
    scale_linetype_manual(values = c( 'Site mechanistic' = 2,'Smoothed' = 1)) +
    scale_x_continuous(limits = c(100,300),
                       breaks =c(121,  182, 244),
                       labels=c('May 1', 'Jul 1','Sep 1')) +
    labs(x = '', y = yaxis) +
    theme_bw() +
    theme(legend.position="none")
  pdf('/Users/dallen@middlebury.edu/My Drive/Research/lyme/2024/pheno_manu/figures/pheno_bysite.pdf',width=7,height=8)  
   p1
  dev.off()
}

## Calculate fraction of larvae questing in early v late
## Both observed and predicted.
## Second results figure
## Inputs: results/full_model_pred.csv
## Output: figures/early_v_late.pdf
{
  all_mod_pred <- read_csv('results/full_model_pred.csv')
  
  switch_date <- all_mod_pred %>%
    group_by(site,mod_type) %>%
    mutate(delta = larva - lag(larva),
           inc = ifelse(delta > 0, 1,-1),
           crit = ifelse(inc*lag(inc) < 0, 'Yes', 'No' ),
           crit = ifelse(lag(larva)==0 | is.na(crit),'No', crit),
           crit_type = ifelse(inc <0, 'Peak', 'Valley')) %>%
    filter(crit == 'Yes', crit_type == 'Valley') %>%
    summarise(julian = max(julian)) %>% ## sometimes two valleys but the second is always the real one
    ungroup() %>%
    select(site, mod_type, switch_jul = julian)
  

  
  summary_data <- all_mod_pred %>%
    left_join(switch_date, by = c('site', 'mod_type')) %>%
    mutate(e_l = ifelse(julian <= switch_jul, 'early', 'late'),
           e_l = ifelse(is.na(e_l), 'early', e_l)) %>% #corrects for runs that don't have a late peak
    group_by(site, mod_type, e_l) %>%
    summarise(
      when = julian[which.max(larva)],
      larva = sum(larva) ) %>%
    pivot_wider(names_from = 'e_l', values_from = c('when', 'larva'), values_fill = 0) # fill in with 0 if no late peak rather than NA
  
  
  summary_data <- summary_data %>% mutate(late_frac = larva_late/(larva_late+larva_early) )
  
  siteelev <- samples %>%
    group_by(site) %>%
    summarise(elev = mean(elev))
  
  p2 <- summary_data %>%
    select(-when_early, -when_late, -larva_early, -larva_late) %>%
    pivot_wider(names_from = mod_type, values_from = late_frac) %>%
    ggplot(aes(Smoothed, `Site mechanistic`)) +
    geom_point() +
    geom_smooth(method = 'lm', se = F, lty = 2, color = 'black') +
    coord_cartesian(xlim=c(0,1), ylim = c(0,1)) +
    geom_line(data = tibble(Smoothed = c(0,1), `Site mechanistic`= c(0,1))) +
    theme_bw() +
    labs(x = 'Observed late-summer fraction', y = 'Predicted late-summer fraction')
  
  
  pdf('/Users/dallen@middlebury.edu/My Drive/Research/lyme/2024/pheno_manu/figures/early_v_late.pdf',width=4,height=4)  
    p2
  dev.off()

  summary_data %>%
    select(-when_early, -when_late, -larva_early, -larva_late) %>%
    pivot_wider(names_from = mod_type, values_from = late_frac) %>% 
    lm(`Site mechanistic` ~ Smoothed, data = .) %>%
    summary()
  
  summary_data %>%
    select(-when_early, -when_late, -larva_early, -larva_late) %>%
    pivot_wider(names_from = mod_type, values_from = late_frac) %>% 
    ungroup() %>%
    summarise(mean(abs(`Site mechanistic`- Smoothed)))
    
}


## Make a map of the sites
{
  require(leaflet)
  site_loc <- tibble(
    site = c('BRF', 'Chipman', 'Cyrstal', 'Foote', 'Frost', 'Gilmore', 'Gorge', 'Jackson', 'Lourie', 'Major', 'SPIN', 'Snowbowl', 'UpperChipman'),
    lat = c( 44.0314,  44.0338,  43.9453,  44.0127,  43.9638,  43.9605,  43.9724,  43.9973,  44.0707,  44.0238,  43.9600,  43.9351,  44.0243),
    lon = c(-73.0820, -73.1609, -72.9675, -73.1371, -73.0034, -72.9795, -73.0738, -73.2011, -73.2605, -73.2613, -73.0288, -72.9511, -73.1626),
    larvae = c(1,1,0,rep(1,10))
  )
  
  site_loc %>%
    mutate(larvae_color = ifelse(larvae, 'black','transparent')) %>%
    leaflet(options = leafletOptions(zoomControl = FALSE)) %>%
    addCircleMarkers(radius = 3, 
                     fillColor = ~larvae_color, 
                     color = 'black', 
                     opacity = 1, 
                     fillOpacity = 1,
                     weight = 2) %>%
    addProviderTiles(providers$Esri.NatGeoWorldMap, options = tileOptions(opacity =  1)) %>% # I also like providers$CartoDB
    setMaxBounds(-73.27,43.92,-72.94,44.08) %>%
    addMiniMap(position = 'topright',
               tiles = providers$Stadia.StamenTonerLite,
               width = 115, height = 115,
               zoomLevelFixed = 5) %>%
    addScaleBar(options = scaleBarOptions(imperial = F))
}