# Analysis Code for 
# "A mechanistic model explains variation in larval  
# tick questing phenology along an elevation gradient"
# Submitted Nov 2020

require(tidyverse)
require(bbmle)
require(cowplot)
set.seed(31416)

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

# This code fits the two-peak phenology curves to data from each elevation category
# Also gets confidence intervals around each parameter
# Input: data/drag_sampling.csv
# Output: results/phenology_fits.RData, 
# this is already included if you don't want to run all this
{
  samples <- read_csv('data/drag_sampling.csv') %>% 
    mutate(elevCat = cut(elev,c(0,200,400,1000),c('low','mid','high')))
  fitList<-list()
  fitList[['low']] <- list(); fitList[['mid']] <- list(); fitList[['high']] <- list() 

  param_names <- c('peak_e', 'tau_e', 'mu_e', 'peak_l', 'tau_l', 'mu_l', 'sigma_l', 'k')
  
  param_guess <- list()
  param_guess[['low']] <- list(peak_e=30, tau_e=160, mu_e=20, peak_l=70, tau_l=200, mu_l=20, sigma_l=0.1,k=0.2)
  param_guess[['mid']] <- list(peak_e=25, tau_e=135, mu_e=35, peak_l=5, tau_l=200, mu_l=60, sigma_l=0.1,k=0.2)
  param_guess[['high']] <- list(peak_e=7, tau_e=170, mu_e=11.5, peak_l=0.03, tau_l=203, mu_l=50, sigma_l=1.5,k=0.03) 
  
  for (which_elev in c('low', 'mid', 'high'))
  {
    data_list <- samples %>% 
      filter(elevCat == which_elev) %>%
      select(day = julian, tickNum = larva) %>%
      as.list()
    fit1 <- mle2(twoPeak, start=param_guess[[which_elev]], data=data_list,method='BFGS')
    fitList[[which_elev]][['fit']] <- fit1
    for (param in param_names) {
      if(param == 'mu_l' & which_elev == 'high')
      {
        fitList[[which_elev]][[param]] <- c(coef(fit1)['mu_l'], coef(fit1)['mu_l']-25,coef(fit1)['mu_l']+25)
      } else {
        fitList[[which_elev]][[param]] <- paramRangeTwoPeak(fit1, param, data_list)
      }  
    }
  }
 save(fitList,file = 'results/phenology_fits.RData')
}

# This code generates 500 phenology fit curves for each elev. cat.
# Uses CIs of fit parameters to do that
# This is to quantify uncertainty in the fits
# Inputs: results/phenology_fits.RData and data/drag_sampling.csv
# Output: results/many_fits.RData (also included if you don't want to run this)
{
  samples <- read_csv('data/drag_sampling.csv') %>% 
    mutate(elevCat = cut(elev,c(0,200,400,1000),c('low','mid','high')))
  load(file = 'results/phenology_fits.RData')
  pheno_smooth <- tibble(julian = numeric(), larva = numeric(), larva_frac = numeric(),fitnum = numeric(),elevCat= character())
  for (which_elev in c('low', 'mid', 'high'))
  {
    elev_data <- samples %>% filter(elevCat == which_elev)
    j<-0
    while (j <500)
    {
      peak_e <- runif(1,fitList[[which_elev]]$peak_e[2],fitList[[which_elev]]$peak_e[3])
      tau_e <- runif(1,fitList[[which_elev]]$tau_e[2],fitList[[which_elev]]$tau_e[3])
      mu_e <- runif(1,fitList[[which_elev]]$mu_e[2],fitList[[which_elev]]$mu_e[3])
      peak_l <- runif(1,fitList[[which_elev]]$peak_l[2],fitList[[which_elev]]$peak_l[3])
      tau_l <- runif(1,fitList[[which_elev]]$tau_l[2],fitList[[which_elev]]$tau_l[3])
      mu_l <- runif(1,fitList[[which_elev]]$mu_l[2],fitList[[which_elev]]$mu_l[3])
      sigma_l <- runif(1,fitList[[which_elev]]$sigma_l[2],fitList[[which_elev]]$sigma_l[3])
      myK <- coef(fitList[[which_elev]]$fit)['k']
      
      does_fit <- ifelse(twoPeak(peak_e, tau_e, mu_e, peak_l, tau_l, mu_l, sigma_l, myK, elev_data$julian, elev_data$larva) < -logLik(fitList[[which_elev]]$fit) +  qchisq(p = 0.95, df = 7)/2,1,0)
  
      twoPeakCurveTemp <- function(x) {twoPeakCurve(x,peak_e=peak_e, tau_e=tau_e, mu_e=mu_e,peak_l=peak_l, tau_l=tau_l, mu_l=mu_l, sigma_l=sigma_l)}
      if (does_fit)
      {
        j <- j + 1
        temp_pred <-tibble(
            julian =1:365, 
            larva = twoPeakCurveTemp(1:365), 
            larva_frac = larva/sum(larva), 
            fitnum = rep(j,365),
            elevCat = which_elev)
        pheno_smooth <- rbind(pheno_smooth,temp_pred)
        print(j)
      }
    }
  }
  save(pheno_smooth, file = 'results/many_fits.RData')
}

# Larval questing phenology model and parameters
# doesn't do anything on its own, 
# but defines required functions for the next section
{
  # Paramters for model, See table 1
  params_with_CI <- list(
    ovi_m = c(132,243),
    ovi_sd = c(37,63),
    ecl_m = c(429,638),
    ecl_sd = c(6,71),
    diapause = c(0.25,0.75),
    hardening = c(14,28),
    start_quest = c(5, 15),
    max_quest = c(20, 30),
    host_find = c(0.004,0.1),
    mort = c(0.007,0.014),
    overwinter_surv = c(0.1,0.8)
  )
  
  # Function to uniformly generate parameters
  rand_param <- function(p)
  {
    toreturn <- list(
      ovi_m = runif(1,p$ovi_m[1],p$ovi_m[2]),
      ovi_sd = runif(1,p$ovi_sd[1],p$ovi_sd[2]),
      ecl_m = runif(1,p$ecl_m[1],p$ecl_m[2]),
      ecl_sd = runif(1,p$ecl_sd[1],p$ecl_sd[2]),
      diapause = runif(1,p$diapause[1],p$diapause[2]),
      hardening = round(runif(1,p$hardening[1],p$hardening[2])),
      max_quest = runif(1,p$max_quest[1],p$max_quest[2]),
      start_quest = runif(1,p$start_quest[1],p$start_quest[2]),
      host_find = runif(1,p$host_find[1],p$host_find[2]),
      mort = runif(1,p$mort[1],p$mort[2]),
      overwinter_surv = runif(1,p$overwinter_surv[1],p$overwinter_surv[2])
    )
    toreturn
  }
  
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
      # add number becoming active on day i to alive pool
      active <- active + start_act[i] 
      # fraction questing on day i 
      questing[i] <- active * f_quest[i]
      # remove ticks from active cohort if they find a host or die
      active <- (1 - params$host_find * f_quest[i] - params$mort) * active
    }
    # fraction of the cohort which survives overwinters
    # overwinter survival times three groups: active at end of questing season; those starting activity after questing season; those that entered diapause earlier
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
}

# This code generates 500 runs of the phenology model for each
# elevation category. It requires the code block above
# Input: data/leaf_litter_temp.csv
# Output: result/many_model_runs.RData
{
  leaf_litter_temp <- read_csv(file = 'data/leaf_litter_temp.csv')
  
  model_pred <- tibble(
    julian = numeric(),
    larva_frac = numeric(),
    elevCat = character(),
    fitnum = numeric())
  
  for (j in 1:500) 
  {
    for (which_elev in c('low', 'mid', 'high'))
    {
      
      elev_cat_temp <- leaf_litter_temp %>%
        filter(elevCat == which_elev) %>%
        group_by(jday) %>%
        summarise(tmean = mean(tmean))
      
      param <- rand_param(params_with_CI) 
      temp_pred <- tibble(
        julian = 1:365,
        larva_frac = larval_quest(elev_cat_temp$tmean[1:365], param)$qst_norm,
        elevCat = which_elev,
        fitnum = j )
      model_pred<- rbind(model_pred,temp_pred)
    }
  }
  
  save(model_pred, file = 'results/many_model_runs.RData')
}

# Make Figure 1, comparing fit and predicted phenology curves
# Fig 1A shows phenology curves fit to data
# Fig 1B shows those fit curves compared to model predictions
# Inputs: data/drag_sampling.csv, results/many_fits.RData, results/many_model_runs.RData
# Output: Figure 1 (pheno_curves.pdf)
{
  load(file = 'results/many_fits.RData') 
  word_to_num <- c('low' = '<200 m', 'mid' = '200-400 m', 'high' = '>400 m')
  pheno_smooth <- pheno_smooth %>% 
    mutate(elevCat = recode(elevCat,!!!word_to_num),
           elevCat = factor(elevCat, levels = c('<200 m', '200-400 m', '>400 m')))
  samples <- read_csv('data/drag_sampling.csv') %>% 
    mutate(elevCat = cut(elev,c(0,200,400,1000),c('<200 m','200-400 m','>400 m')),
           elevCat = factor(elevCat, levels = c('<200 m', '200-400 m', '>400 m')),
           fitnum = 1)
  
  ylab <- expression(paste('Larvae (per 200 ',m^2,')'))
  # this is Figure 1A
  fit_pheno_plot <- samples %>%
    ggplot(aes(julian,larva,group = fitnum)) +
    geom_point(cex=0.25) +
    facet_wrap(~elevCat) +
    theme_classic() +
    theme(strip.background = element_rect(color='transparent'),
          strip.text = element_text(color='black',size = 10),
          axis.text = element_text(color='black',size = 10),
          plot.margin = unit(c(0,0,-0.35,0), 'cm'),
          axis.title = element_text(size = 10)) +
    geom_path(data=pheno_smooth,aes(julian,larva,group=fitnum), alpha = 0.015) +
    scale_x_continuous(limits = c(105,305),
                       breaks =c(121,  182, 244, 305),
                       labels=c('', '','', '')) +
    scale_y_continuous(limits = c(0,200)) +
    labs(x='',y=ylab)

  load(file = 'results/many_model_runs.RData') 
  
  pheno_smooth <- pheno_smooth %>%
    mutate(type = 'observed') %>%
    select(-larva)

  model_pred <- model_pred %>%
    mutate(type = 'modeled') %>%
    rbind(pheno_smooth) %>%
    mutate(elevCat = recode(elevCat,!!!word_to_num),
           elevCat = factor(elevCat, levels = c('<200 m', '200-400 m', '>400 m'))) %>%
    arrange(fitnum,type,julian) 
  
  mean_vals <- model_pred %>%
    group_by(julian, elevCat, type) %>%
    summarise(larva_frac = mean(larva_frac)) %>%
    mutate(fitnum = 1) %>%
    arrange(type,julian)

  # Make Figure 1B
  mod_v_obs_plot <- model_pred %>%
    ggplot(aes(julian,larva_frac, group = fitnum, color=type)) +
    geom_path(alpha = 0.015) +
    geom_path(data = mean_vals,lwd=0.75) +
    facet_wrap(~elevCat) +
    theme_classic() +
    theme(strip.background = element_rect(color='transparent'),
          strip.text = element_blank(),
          axis.text = element_text(color='black',size = 10),
          axis.title = element_text(size = 10),
          plot.margin = unit(c(-0.35,0,0,0), 'cm'),
          legend.position = 'none') +
    scale_color_manual(values = c( '#d7191c', '#2b83ba')) +
    scale_x_continuous(limits = c(105,305),
                       breaks =c(121,  182, 244, 305),
                       labels=c('May 1', 'Jul 1','Sep 1', 'Nov 1')) +
    labs(x='',y='Fraction questing')
 
  # Combine subfigures into Figure 1
  pdf('figures/pheno_curves.phdf',width=7.25,height=6)
    plot_grid(fit_pheno_plot, 
              mod_v_obs_plot,
              ncol = 1, 
              labels = c('A', 'B'), 
              align = 'v')
  dev.off()
}

# Make Figure 2, comparing day and height of peaks for 
# fit phenology to model predicted phenology
# Fig 2A shows comparison of day of early and late peaks
# Fig 2B shows comparison of height of early and late peaks
# Inputs: results/many_fits.RData, results/many_model_runs.RData
# Output: Figure 2 (peak_comp.pdf)
{
  load(file = 'results/many_fits.RData') 
  load(file = 'results/many_model_runs.RData') 

  word_to_num <- c('low' = '<200 m', 'mid' = '200-400 m', 'high' = '>400 m')
  
  pheno_smooth <- pheno_smooth %>% 
    mutate(type = 'observed',
           elevCat = recode(elevCat,!!!word_to_num)) %>% 
    select(-larva) 
  model_pred <- model_pred %>% 
    mutate(type = 'modelled',
           elevCat = recode(elevCat,!!!word_to_num))
  
  # for each fit or model find the day when the late questing starts 
  find_break <- pheno_smooth %>%
    rbind(model_pred) %>%
    group_by(fitnum, elevCat, type) %>%
    mutate(delta = larva_frac - lag(larva_frac, default = 0),
           mid_point = ifelse(delta > 0 & lag(delta) < 0, 1, 0) ) %>%
    filter(mid_point == 1) %>%
    select(fitnum, elevCat, type, e_to_l_day = julian)
           
  # summarize each fit or model run to find timing and height of both peaks 
  peak_comp <- pheno_smooth %>%
    rbind(model_pred) %>%
    left_join(find_break, by = c('fitnum', 'elevCat', 'type')) %>%
    mutate(e_or_l = ifelse(julian < e_to_l_day | is.na(e_to_l_day),'early', 'late')) %>%
    group_by(elevCat, fitnum, e_or_l, type) %>%
    summarise(peak = max(larva_frac), when = julian[which.max(larva_frac)]) %>%
    group_by(elevCat, e_or_l, type) %>%
    mutate(e_or_l = factor(e_or_l,levels=c('late','early'))) 
  
  # some fits/runs don't have a late peak need to go back and add that in as 0
  # this first makes all combos
  peak_comp_0 <- expand_grid(
    fitnum = 1:500,
    elevCat = c('<200 m', '200-400 m', '>400 m'),
    type = c('observed', 'modelled'),
    e_or_l = c('early', 'late')
  )
  
  # now fill in all combos with values, missing peak heights become 0
  peak_comp <- peak_comp_0 %>%
    left_join(peak_comp, by=c('fitnum','elevCat','type', 'e_or_l') ) %>%
    mutate(peak = ifelse(is.na(peak), 0, peak),
           elevCat = factor(elevCat, levels = c('<200 m', '200-400 m', '>400 m')) )
  
  peak_comp_e <- peak_comp %>% filter(e_or_l == 'early')
  peak_comp_l <- peak_comp %>% filter(e_or_l == 'late')
    
  # This makes subfigure 2A
  when_plot <- peak_comp_l %>%
    ggplot(aes(x = when, color = type, fill = type, group= type, y = ..scaled..)) +
    geom_density(alpha = 0.5) +
    geom_density(data = peak_comp_e,aes(y=-..scaled..,),alpha=0.5) +
    facet_wrap(~elevCat) +
     scale_color_manual(values = c( '#d7191c', '#2b83ba')) +
    scale_fill_manual(values = c( '#d7191c', '#2b83ba')) +
    coord_flip() + 
    theme_classic() +
    theme(strip.background = element_rect(color='transparent'),
          strip.text =  element_text(color='black',size = 10),
          axis.text = element_text(color='black',size = 10),
          axis.title = element_text(size = 10),
          legend.position = 'none') +
    scale_x_continuous(limits = c(105,350),
                       breaks =c(121,  182, 244, 305),
                       labels=c('May 1', 'Jul 1','Sep 1', 'Nov 1')) +
    scale_y_continuous(breaks = NULL) +
    geom_hline(yintercept = 0) +
    labs(y='',x='Day of peak')
    
  # This makes subfigure 2B
  peak_plot <- peak_comp_l %>%
    ggplot(aes(x = peak, color = type, fill = type, group= type, y = ..scaled..)) +
    geom_density(alpha = 0.5) +
    geom_density(data = peak_comp_e,aes(y=-..scaled..,),alpha=0.5) +
    facet_wrap(~elevCat) +
    scale_color_manual(values = c( '#d7191c', '#2b83ba')) +
    scale_fill_manual(values = c('#d7191c', '#2b83ba')) +
    coord_flip() + 
    theme_classic() +
    theme(strip.background = element_rect(color='transparent'),
          strip.text = element_blank(),
          axis.text = element_text(color='black',size = 10),
          axis.title = element_text(size = 10),
          legend.position = 'none') +
    scale_y_continuous(breaks = c(-0.5,0.5), labels = c('Early', 'Late')) +
    geom_hline(yintercept = 0) +
    labs(y='',x='Larval fraction at peak')
    
  # Make figure 2  
  pdf('figures/peak_comp.pdf',width=3.5,height=6)
    plot_grid(when_plot, 
              peak_plot,
              ncol = 1, 
              labels = c('A', 'B'), 
              align = 'v')
  dev.off()
  
}

# Makes Figure A2, which shows sampling results by site
# Input: data/drag_sampling.csv
{
  samples <- read_csv('data/drag_sampling.csv')
  
  # filter out the two sites which never had any larval ticks
  samplesMod <- samples %>%
    filter(!(site %in% c('Snowbowl', 'Crystal'))) %>%
    mutate(elevText = paste(elev, ' m'))    
  
  # dummy data so plot looks nice
  newRow <- tibble(
    site = NA,
    elev = NA,
    date = NA,
    julian = NA,
    larva = NA,
    elevText = ''
  )
  
  samplesMod <- rbind(samplesMod, newRow)
  tempLevels <- levels(as.factor(samplesMod$elevText))
  samplesMod <- samplesMod %>% rbind(samplesMod,newRow) %>%
    mutate(elevText = factor(elevText,tempLevels[c(2:8,1,9:12)]))
  
  
  ylab <- expression(paste('Larvae (per 200 ',m^2,')'))
  
  p1 <- samplesMod %>% 
    filter(elev < 200) %>%
    ggplot(aes(julian, larva)) +
    geom_point(cex = 0.5) +
    facet_wrap(~elevText, nrow = 1) +
    stat_smooth(se = F) +
    coord_cartesian(ylim = c(0,125)) +
    theme_classic() +
    theme(axis.text = element_text(color='black',size = 10),
          axis.title = element_text(size = 10),
          plot.margin = unit(c(0,0,-0.15,0), 'cm'),
          axis.text.x = element_blank()) +
    labs(x='',y='') +
    scale_x_continuous(limits = c(105,305),
                       breaks =c(121,  182, 244),
                       labels = c('','',''))
  
  p2 <- samplesMod %>% 
    filter((elev > 200 & elev < 400) | elevText == '') %>%
    ggplot(aes(julian, larva)) +
    geom_point(cex = 0.5) +
    facet_wrap(~elevText, nrow = 1) +
    stat_smooth(se = F) +
    coord_cartesian(ylim = c(0,60)) +
    theme_classic() +
    theme(axis.text = element_text(color='black',size = 10),
          axis.title = element_text(size = 10),
          plot.margin = unit(c(-0.15,0,-0.15,0), 'cm'),
          axis.text.x = element_blank()) +
    labs(x='',y=ylab) +
    scale_x_continuous(limits = c(105,305),
                       breaks =c(121,  182, 244),
                       labels = c('','',''))
  
  p3 <- samplesMod %>% 
    filter(elev > 400) %>%
    ggplot(aes(julian, larva)) +
    geom_point(cex = 0.5) +
    facet_wrap(~elevText,  nrow = 1) +
    stat_smooth(se = F) +
    coord_cartesian(ylim = c(0,30)) +
    theme_classic() +
    theme(axis.text = element_text(color='black',size = 10),
          plot.margin = unit(c(-0.15,0,0,0), 'cm'),
          axis.title = element_text(size = 10)) +
    labs(x='',y='') +
    scale_x_continuous(limits = c(105,305),
                       breaks =c(121,  182, 244),
                       labels=c('May 1', 'Jul 1','Sep 1'))
  

  
  pdf('figures/pheno_bysite.pdf',width=6,height=5)  
    plot_grid(p1,p2,p3, align = 'v', ncol = 1)
  dev.off()
  
}
