## David Allen
## Code for larval tick phenology manuscript
## April 2020

require(tidyverse)
require(bbmle)

# required phenology functions
{
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

# fit curves to phenology patterns for three elevation classes
{
  samples <- read_csv('data/drag_samplingwith2020.csv') %>% mutate(elevCat = cut(elev,c(0,200,410,1000),c('low','mid','high')))
  fitList<-list()
  fitList[['low']] <- list(); fitList[['mid']] <- list(); fitList[['high']] <- list() 

  param_of_interest <- c('peak_e', 'tau_e', 'mu_e', 'peak_l', 'tau_l', 'mu_l', 'sigma_l', 'k')
  julianSt <- 120
  julianEnd <- 305
  
  paramGuess <- list()
  paramGuess[['low']] <- list(peak_e=30, tau_e=160, mu_e=20, peak_l=70, tau_l=200, mu_l=20, sigma_l=0.1,k=0.2)
  paramGuess[['mid']] <- list(peak_e=25, tau_e=135, mu_e=35, peak_l=5, tau_l=200, mu_l=60, sigma_l=0.1,k=0.2)
  paramGuess[['high']] <- list(peak_e=7, tau_e=170, mu_e=11.5, peak_l=0.03, tau_l=203, mu_l=50, sigma_l=1.5,k=0.03) 
  
  for (which_elev in c('low', 'mid', 'high'))
  {
    subSetData <- samples %>% filter(elevCat == which_elev)
    dataList <- with(subSetData,list(day = julian, tickNum = larva))
    fit1 <- mle2(twoPeak, start=paramGuess[[which_elev]], data=dataList,method='BFGS')
    # with(subSetData,plot(julian,larva))
    # curve(twoPeakCurve,from=julianSt,to=julianEnd,add=TRUE)
    fitList[[which_elev]][['fit']] <- fit1
    for (param in param_of_interest) {
      if(param == 'mu_l' & which_elev == 'high')
      {
        fitList[[which_elev]][[param]] <- c(coef(fit1)['mu_l'], coef(fit1)['mu_l']-25,coef(fit1)['mu_l']+25)
      } else {
        fitList[[which_elev]][[param]] <- paramRangeTwoPeak(fit1, param, dataList)
      }  
    }
  }
 save(fitList,file = 'data/phenology_fits.RData')
}


# get range of phenology fits
{
  samples <- read_csv('data/drag_samplingwith2020.csv') %>% mutate(elevCat = cut(elev,c(0,200,410,1000),c('low','mid','high')))
  load(file = 'data/phenology_fits.RData')
  pred <- tibble(julian = numeric(), larva = numeric(), larva_frac = numeric(),fitnum = numeric(),elevCat= character())
  for (which_elev in c('low', 'mid', 'high'))
  {
    subSetData <- samples %>% filter(elevCat == which_elev)
    j<-0
    while (j <100)
    {
      peak_e <- runif(1,fitList[[which_elev]]$peak_e[2],fitList[[which_elev]]$peak_e[3])
      tau_e <- runif(1,fitList[[which_elev]]$tau_e[2],fitList[[which_elev]]$tau_e[3])
      mu_e <- runif(1,fitList[[which_elev]]$mu_e[2],fitList[[which_elev]]$mu_e[3])
      peak_l <- runif(1,fitList[[which_elev]]$peak_l[2],fitList[[which_elev]]$peak_l[3])
      tau_l <- runif(1,fitList[[which_elev]]$tau_l[2],fitList[[which_elev]]$tau_l[3])
      mu_l <- runif(1,fitList[[which_elev]]$mu_l[2],fitList[[which_elev]]$mu_l[3])
      sigma_l <- runif(1,fitList[[which_elev]]$sigma_l[2],fitList[[which_elev]]$sigma_l[3])
      myK <- coef(fitList[[which_elev]]$fit)['k']
      
      doesFit <- ifelse(twoPeak(peak_e, tau_e, mu_e, peak_l, tau_l, mu_l, sigma_l, myK, subSetData$julian, subSetData$larva) < -logLik(fitList[[which_elev]]$fit) +  2,1,0)
  
      twoPeakCurveTemp <- function(x) {twoPeakCurve(x,peak_e=peak_e, tau_e=tau_e, mu_e=mu_e,peak_l=peak_l, tau_l=tau_l, mu_l=mu_l, sigma_l=sigma_l)}
      if (doesFit)
      {
        j <- j + 1
        temp_pred <-tibble(
            julian =1:365, 
            larva = twoPeakCurveTemp(1:365), 
            larva_frac = larva/sum(larva), 
            fitnum = rep(j,365),
            elevCat = which_elev)
        pred <- rbind(pred,temp_pred)
        print(j)
      }
      
    }
  }
  pheno_smooth <- pred  %>% mutate(elevCat = factor(elevCat,levels=c('low','mid','high')))
  save(pheno_smooth, file = 'data/smoothed_pheno.RData')
  
}


# make elevation versus fraction early summer larva plot
{
  samples <- read_csv('data/drag_sampling.csv') %>% mutate(elevCat = cut(elev,c(0,200,410,1000),c('low','mid','high')))
  late <- samples %>%
    filter(julian > 212) %>%
    group_by(site) %>%
    summarise(fall_l = mean(larva),
              elev = unique(elev))
  early <- samples %>%
    filter(julian <= 212) %>%
    group_by(site) %>%
    summarise(spring_l = mean(larva))
  
  larva_comparision <- late %>%
    inner_join(early,by='site') %>%
    mutate(early_frac = spring_l/(spring_l+fall_l) )
  
  larva_comparision %>%
    filter(is.finite(early_frac)) %>%
    lm(early_frac ~ elev, data = .) %>%
    summary()
  
  pdf('figures/elev_l_frac.pdf',width=3.14,height=2.5)
    larva_comparision %>% 
      ggplot(aes(elev,early_frac)) +
      geom_point() +
      theme_classic() +
      theme(axis.text = element_text(color='black', size =10),
            axis.title = element_text(size = 10)) +
      labs(x = 'Elevation (m)', y = 'Early summer fraction') +
      stat_smooth(se=F,method='lm',color='black') + 
      coord_cartesian(xlim=c(120,600))
  dev.off()
  
}

# make phenology plot with fit curves
{
  load(file = 'data/smoothed_pheno.RData') 
  samples <- read_csv('data/drag_samplingwith2020.csv') %>% 
    mutate(elevCat = cut(elev,c(0,200,410,1000),c('low','mid','high')),
           fitnum = 1)
   
  labv2 <- tibble(elevCat = c('low','mid','high'),
                  julian=rep(212,3),
                  larva=rep(200,3),
                  lab = c('<200 m', '200 - 400 m', '>400 m'),
                  fitnum = 1) %>%
    mutate(elevCat = factor(elevCat,levels=c('low','mid','high')))
  ylab <- expression(paste('Larvae (per 200 ',m^2,')'))
  pdf('figures/l_pheno_fit_CI.pdf',width=6.68,height=3)
    samples %>%
      ggplot(aes(julian,larva,group = fitnum)) +
      geom_point(cex=0.25) +
      facet_wrap(~elevCat) +
      theme_classic() +
      theme(strip.background = element_rect(color='transparent'),
            strip.text = element_blank(),
            axis.text = element_text(color='black',size = 10),
            axis.title = element_text(size = 10)) +
      geom_path(data=pheno_smooth,aes(julian,larva,group=fitnum), alpha = 0.05) +
      scale_x_continuous(limits = c(100,300),
                         #                   breaks =c(121, 152, 182, 213, 244, 274),
                         #                  labels=c('May 1', 'Jun 1', 'Jul 1', 'Aug 1','Sep 1', 'Oct 1')) +
                         breaks =c(121,  182, 244),
                         labels=c('May 1', 'Jul 1','Sep 1')) +
      scale_y_continuous(limits = c(0,200)) +
      labs(x='',y=ylab) +
      geom_text(data=labv2,aes(label=lab),size=4)
  dev.off()
  
}

# make plot comparing parameters
{
  load(file = 'data/phenology_fits.RData') # see above for code which fit these functions
  # compare parameters
  paramcomp <- tibble(
    elev = c('<200 m','<200 m','200-400 m','200-400 m','>400 m','>400 m'),
    when = rep(c('Early','Late'),3),
    param = c(fitList$low$peak_e[1],fitList$low$peak_l[1],fitList$mid$peak_e[1],fitList$mid$peak_l[1],fitList$high$peak_e[1],fitList$high$peak_l[1]),
    upper_param = c(fitList$low$peak_e[3],fitList$low$peak_l[3],fitList$mid$peak_e[3],fitList$mid$peak_l[3],fitList$high$peak_e[3],fitList$high$peak_l[3]),
    lower_param =c(fitList$low$peak_e[2],fitList$low$peak_l[2],fitList$mid$peak_e[2],fitList$mid$peak_l[2],fitList$high$peak_e[2],fitList$high$peak_l[2]),
  ) %>% 
    mutate(elev= factor(elev,levels = c('<200 m','200-400 m','>400 m')))
  
  pdf('figures/peak_param_comp.pdf',width=3.14,height=3)
    ylab <- expression(paste('Peak (larvae per 200 ',m^2,')'))
    paramcomp %>%
      ggplot(aes(when,param)) +
      geom_point() +
      facet_wrap(~elev) +
      geom_errorbar(aes(ymin=lower_param,ymax=upper_param), width=0.5) +
      theme_classic() +
      theme(strip.background = element_rect(color='transparent'),
            axis.text = element_text(color='black',size = 10),
            axis.title = element_text(size = 10),
            strip.text = element_text(size = 10)) +
      labs(x='',y=ylab) +
      scale_y_log10(breaks=c(1e-2,1e-1,1,1e1,1e2),
                    labels = c(.001,0.1,1,10,100))    
  dev.off()
}

# larval development and quest model
{
  
  params_with_CI <- list(
    ovi_m = c(131.8,330.5),
    ovi_sd = c(57.6,69.8),
    ecl_m = c(363.4,469.2),
    ecl_sd = c(63.2,80.5),
    hardening = c(14,28),
    quest25 = c(0.49, 1),
    start_quest = c(7.0, 10),
    host_find = c(0.01,0.03),
    mort = c(0.006,0.011),
    diapause = c(0.25,0.75),
    overwinter_surv = c(0.104,0.44)
  )
  
  rand_param <- function(p)
  {
    toreturn <- list(
      ovi_m = runif(1,p$ovi_m[1],p$ovi_m[2]),
      ovi_sd = runif(1,p$ovi_sd[1],p$ovi_sd[2]),
      ecl_m = runif(1,p$ecl_m[1],p$ecl_m[2]),
      ecl_sd = runif(1,p$ecl_sd[1],p$ecl_sd[2]),
      hardening = round(runif(1,p$hardening[1],p$hardening[2])),
      quest25 = runif(1,p$quest25[1],p$quest25[2]),
      start_quest = runif(1,p$start_quest[1],p$start_quest[2]),
      host_find = runif(1,p$host_find[1],p$host_find[2]),
      mort = runif(1,p$mort[1],p$mort[2]),
      diapause = runif(1,p$diapause[1],p$diapause[2]),
      overwinter_surv = runif(1,p$overwinter_surv[1],p$overwinter_surv[2])
    )
    toreturn
  }
  
  # parameters for model
  params <- list(
    ovi_m = 243.2, # in degree days base 6 C Rand et al. 2004
    ovi_sd = 63.1, # in degree days base 6 C Rand et al. 2004
    ecl_m = 428.5, # in degree days base 11 C Rand et al. 2004
    ecl_sd = 70.8, # in degree days base 11 C Rand et al. 2004
    hardening = 21, # in days Daniels et al. 1996
    start_quest = 10, # in C Ogden et al. 2005
    max_quest = 25, # in C Ogden et al. 2005
    quest_slope1 = 0.067, # in 1/C Odgen et al. 2005
    quest_slope2 = -0.067, # in 1/C Odgen et al. 2005
    host_find = 0.0207, # in 1/day Odgen et al. 2005
    mort = 0.006, # in 1/day Odgen et al. 2005
    diapause = 0.25, # Ogden et al. 2018
    overwinter_surv = 0.44 # Lindsay et al. 1995
  )
  
  # takes in daily average temp at leaf litter and parameters
  # outputs fraction of larvae oviposited, ecolosed, questing on each day
  larval_quest <- function(tmean, params)
  {
    tot_day <- length(tmean)
    dd6 <- cumsum(ifelse(tmean > 6, tmean - 6, 0))
    # fraction of cohort oviposited on each day
    ovi = pnorm(dd6, params$ovi_m, params$ovi_sd) - pnorm(lag(dd6, n = 1, default = 0), params$ovi_m, params$ovi_sd)             
    
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
    act_prediapause <- lag(ecl, params$hardening, 0) # hardening time
    # all ticks eclosing before solstice become active, fraction of htem ecolsing afterwards enter diapause (Ogden et al. 2018)
    diapause_frac_daily <- c(rep(0, sum_sol + params$hardening), rep(params$diapause, tot_day - (sum_sol + params$hardening)))
    # fraction of cohort starting activity on each day
    start_act <- (1 - diapause_frac_daily) * act_prediapause
    # fraction of cohort that enters diapause after ecolosion
    diapause_frac_total <- sum(diapause_frac_daily*act_prediapause)
    
    # this section calculates the fraction of cohort questing on each day
    # sums to greater than 1 since a tick can quest on more than one day
    questing <- rep(0, tot_day)
    days_to_winter <- min(which(tmean<params$start_quest)[which(tmean<params$start_quest) > sum_sol]) #here by 'winter' i just mean when temp falls below that which larvae quest
    active <- 0 # tracks number of active ticks
    
    # fraction of ticks active on a given day that will quest based on temp
    quest_slope <- params$quest25/(25-params$start_quest)
    f_quest <- ifelse(tmean < 25, 
                      quest_slope*(tmean-25) + params$quest25,
                      -quest_slope*(tmean-25) + params$quest25 )
    f_quest <- ifelse(f_quest>1,1,ifelse(f_quest<0,0,f_quest))
    
    for (i in 1:days_to_winter)
    {
      active <- active + start_act[i] # add number becoming active on day i to alive pool
      # fraction question on day i 
      questing[i] <- active * f_quest[i]
      # remove ticks from active cohort if they find a host or die
      active <- (1 - params$host_find * f_quest[i] - params$mort) * active
    }
    # fraction of the cohort which survives overwinters
    # overwinter survival times three groups: active at end of questing season; those starting activity after questing season; those that entered diapause earlier
    larvae_overwinter <- params$overwinter_surv * (active + sum(start_act[(days_to_winter+1):tot_day]) + diapause_frac_total)
    
    # add overwintered larvae to questing vector
    for (i in 1:days_to_winter)
    {
      questing[i] <- questing[i] + larvae_overwinter * f_quest[i]
      larvae_overwinter <- (1 - params$host_find * f_quest[i] - params$mort) * larvae_overwinter
    }
    
    output_tibble <- tibble(
      ovi = ovi, # fraction of cohort oviposited each day
      ecl = ecl, # fraction of cohort eclosed each day
      qst = questing, # fraction of cohort questing each day (sums > 1)
      qst_norm = questing/sum(questing) # fraction of questing larvae on a given day, for comparison to real data
    )
    
    return(output_tibble)
  }
}

# make plot comparing model predictions to smoothed observed larval phenology
{
  load(file = 'data/smoothed_pheno.RData') 
  load(file = 'data/processed_prism.RData')
  
  model_pred <- tibble(
    julian = numeric(),
    larva_frac = numeric(),
    elevCat = character(),
    fitnum = numeric())
  
  for (j in 1:100) 
  {
    for (which_elev in c('low', 'mid', 'high'))
    {
      
      if (which_elev == 'low') temp_climate <- siteClimate %>% filter(site == 'Foote')
      if (which_elev == 'mid') temp_climate <- siteClimate %>% filter(site == 'Gorge')
      if (which_elev == 'high') temp_climate <- siteClimate %>% filter(site == 'Snowbowl')
      
      param <- rand_param(params_with_CI) 
      temp_pred <- tibble(
        julian = 1:365,
        larva_frac = larval_quest(temp_climate$tmean[1:365], param)$qst_norm,
        elevCat = which_elev,
        fitnum = j )
      model_pred<- rbind(model_pred,temp_pred)
    }
  }
  
  
  model_pred <- model_pred %>%
    mutate(elevCat = factor(elevCat,levels=c('low','mid','high')))
  
  labv3 <- tibble(elevCat = c('low','mid','high'),
                  day=rep(175,3),
                  larva_frac=rep(0.035,3),
                  lab = c('<200 m', '200 - 400 m', '>400 m')) %>%
    mutate(elevCat = factor(elevCat,levels=c('low','mid','high')))
  
  pdf('figures/pheno_v_mod_CI.pdf',width=6.68,height=3)
    pheno_smooth %>%
      ggplot(aes(julian,larva_frac, group = fitnum)) +
      geom_path(alpha = 0.1) +
      facet_wrap(~elevCat) +
      theme_classic() +
      theme(strip.background = element_rect(color='transparent'),
            strip.text = element_blank(),
            axis.text = element_text(color='black',size = 10),
            axis.title = element_text(size = 10),
            legend.position = c(0.5,0.5)) +
      scale_x_continuous(limits = c(100,300),
                         #                   breaks =c(121, 152, 182, 213, 244, 274),
                         #                  labels=c('May 1', 'Jun 1', 'Jul 1', 'Aug 1','Sep 1', 'Oct 1')) +
                         breaks =c(121,  182, 244),
                         labels=c('May 1', 'Jul 1','Sep 1')) +
      labs(x='',y='Fraction questing') +
      #geom_text(data=labv3,aes(label=lab),cex=4) +
      geom_path(data=model_pred,color='red',alpha = 0.1)
  dev.off()

}




# make plot comparing model predictions to smoothed observed larval phenology by site
{
  load(file = 'data/phenology_fits_bysite.RData') 
  load(file = 'data/processed_prism.RData')
  
  
  mysites <- names(fitList)
  fit1 <- fitList[[ mysites[1] ]]$fit
  smooth <- twoPeakCurve(1:365)
  fit_larvae <- tibble(
    day = 1:365,
    larva = smooth,
    larva_frac = smooth/sum(smooth),
    site = rep(mysites[1],365)
  )
  
  for (s in mysites[2:length(mysites)])
  {
    fit1 <- fitList[[ s ]]$fit
    smooth <- twoPeakCurve(1:365)
    temp_larvae <- tibble(
      day = 1:365,
      larva = smooth,
      larva_frac = smooth/sum(smooth),
      site = rep(s,365)
    )
    fit_larvae <- rbind(fit_larvae,temp_larvae)
  }
  
  fit_larvae <- fit_larvae %>%
    mutate(site = factor(site, levels = c('Foote','Chipman','Snake','Gorge','Chipman2','BRF','SPIN', 'Frost','Gilmore')))    
  
  
  # paramters for model
  params <- list(
    ovi_m = 243.2, # in degree days base 6 C Rand et al. 2004
    ovi_sd = 63.1, # in degree days base 6 C Rand et al. 2004
    ecl_m = 428.5, # in degree days base 11 C Rand et al. 2004
    ecl_sd = 70.8, # in degree days base 11 C Rand et al. 2004
    hardening = 21, # in days Daniels et al. 1996
    start_quest = 10, # in C Ogden et al. 2005
    max_quest = 25, # in C Ogden et al. 2005
    quest_slope1 = 0.067, # in 1/C Odgen et al. 2005
    quest_slope2 = -0.067, # in 1/C Odgen et al. 2005
    host_find = 0.0207, # in 1/day Odgen et al. 2005
    mort = 0.006, # in 1/day Odgen et al. 2005
    diapause = 0.25, # Ogden et al. 2018
    overwinter_surv = 0.44 # Lindsay et al. 1995
  )
  
  pred_larvae <- tibble(
    day = 1:365,
    larva_frac = larval_quest(siteClimate[[ mysites[1] ]]$tmean[1:365], params)$qst_norm,
    site = rep(mysites[1],365)
  )
  
  for (s in mysites[2:length(mysites)])
  {
    temp_pred <- tibble(
      day = 1:365,
      larva_frac = larval_quest(siteClimate[[ s ]]$tmean[1:365], params)$qst_norm,
      site = rep(s,365)
    )
    pred_larvae <- rbind(pred_larvae, temp_pred)
  }
  
  pred_larvae <-
    pred_larvae %>%
    mutate(site = factor(site, levels = c('Foote','Chipman','Snake','Gorge','Chipman2','BRF','SPIN', 'Frost','Gilmore')))    
  
  labv3 <- tibble(elevCat = c('low','mid','high'),
                  day=rep(175,3),
                  larva_frac=rep(0.035,3),
                  lab = c('<200 m', '200 - 400 m', '>400 m')) %>%
    mutate(elevCat = factor(elevCat,levels=c('low','mid','high')))
  
  pdf('figures/pheno_v_mod_bysite.pdf',width=6,height=5)
  fit_larvae %>%
    ggplot(aes(day,larva_frac)) +
    geom_path() +
    facet_wrap(~site) +
    theme_classic() +
    theme(strip.background = element_rect(color='transparent'),
          strip.text = element_blank(),
          axis.text = element_text(color='black',size = 10),
          axis.title = element_text(size = 10),
          legend.position = c(0.5,0.5)) +
    scale_x_continuous(limits = c(100,300),
                       #                   breaks =c(121, 152, 182, 213, 244, 274),
                       #                  labels=c('May 1', 'Jun 1', 'Jul 1', 'Aug 1','Sep 1', 'Oct 1')) +
                       breaks =c(121,  182, 244),
                       labels=c('May 1', 'Jul 1','Sep 1')) +
    labs(x='',y='Fraction questing') +
   # geom_text(data=labv3,aes(label=lab),cex=4) +
    geom_path(data=pred_larvae,lty=2)
  dev.off()
  
}





# fit curves to phenology patterns for each site
{
  samples <- read_csv('data/drag_samplingwith2020.csv')
  
  samplesMod <- samples %>%
    mutate(site = ifelse(site == 'BRF2','BRF',site),
           site = ifelse(site %in% c('Lourie','Major'),'Snake',site)) %>%
    filter(!(site %in% c('Snowbowl', 'Crystal'))) %>%
    mutate(site = factor(site, levels = c('Foote','Chipman','Snake','Gorge','Chipman2','BRF','SPIN', 'Frost','Gilmore')))    
  
  samplesMod %>%
    ggplot(aes(julian,larva)) +
    geom_point() +
    facet_wrap(~site) +
    coord_cartesian(ylim=c(0,75)) +
    stat_smooth(se=F)
  
  fitList<-list()
  mysites <- samplesMod %>% pull(site) %>% unique() 
  mysites2 <- mysites[!mysites %in% c('Snowbowl', 'Crystal')]
  
  for (s in mysites2) fitList[[s]] <- list()
  
  param_of_interest <- c('peak_e', 'tau_e', 'mu_e', 'peak_l', 'tau_l', 'mu_l', 'sigma_l', 'k')
  julianSt <- 120
  julianEnd <- 305
  
  paramGuess <- list()
  paramGuess[['low']] <- list(peak_e=30, tau_e=150, mu_e=20, peak_l=70, tau_l=200, mu_l=20, sigma_l=0.1,k=0.2)
  paramGuess[['mid']] <- list(peak_e=5, tau_e=170, mu_e=10, peak_l=5, tau_l=195, mu_l=50, sigma_l=0.1,k=0.2)
  paramGuess[['high']] <- list(peak_e=7, tau_e=170, mu_e=11.5, peak_l=0.03, tau_l=203, mu_l=50, sigma_l=1.5,k=0.03) 
  do_conf_int <- TRUE  
  
  for (s in mysites2)
  {
    subSetData <- samplesMod %>% filter(site == s)
    dataList <- with(subSetData,list(day = julian, tickNum = larva))
    which_elev <- ifelse(s %in% c('Foote', 'Snake', 'Chipman'), 'low', ifelse(s %in% c('Frost','Gilmore','SPIN'),'high','mid'))
    fit1 <- mle2(twoPeak, start=paramGuess[[which_elev]], data=dataList,method='BFGS')
    # with(subSetData,plot(julian,larva))
    # curve(twoPeakCurve,from=julianSt,to=julianEnd,add=TRUE)
    fitList[[s]][['fit']] <- fit1
    if (do_conf_int)
    {
      for (param in param_of_interest) {
        if (which_elev == 'high' & param == 'mu_l')
        {
          fitList[[s]][[param]] <- c(coef(fit1)['mu_l'], coef(fit1)['mu_l']-25,coef(fit1)['mu_l']+25)
        } else {
          fitList[[s]][[param]] <- paramRangeTwoPeak(fit1, param, dataList)
          print(param)
        }
      }
    }
    print(s)
  }
  
  s <- 'Foote'
  subSetData <- samplesMod %>% filter(site == s)
  fit1 <- fitList[[s]][['fit']]
  pred <- tibble(julian =julianSt:julianEnd, larva = twoPeakCurve(julianSt:julianEnd))
  
  ggplot(subSetData,aes(julian,larva)) +
    geom_point() +
    stat_smooth(se=F,span = 0.25) +
    geom_path(data = pred, color = 'red',lwd=1.5)
  
  save(fitList,file = 'data/phenology_fits_bysite.RData')
}



# make phenology plot with fit curves by site
{
  load(file = 'data/phenology_fits_bysite.RData') # see above for code which fit these functions
  samples <- read_csv('data/drag_samplingwith2020.csv') 
  samplesMod <- samples %>%
    mutate(site = ifelse(site == 'BRF2','BRF',site),
           site = ifelse(site %in% c('Lourie','Major'),'Snake',site)) %>%
    filter(!(site %in% c('Snowbowl', 'Crystal'))) %>%
    mutate(site = factor(site, levels = c('Foote','Chipman','Snake','Gorge','Chipman2','BRF','SPIN', 'Frost','Gilmore')))    
  
  
  
  mysites <- names(fitList)
  
  fit1 <- fitList[[ mysites[1] ]]$fit
  smooth <- twoPeakCurve(1:365)
  
  fit_larvae <- tibble(
    day = 1:365,
    larva = smooth,
    larva_frac = smooth/sum(smooth),
    site = rep(mysites[1],365)
  )
  
  for (s in mysites[2:length(mysites)])
  {
    fit1 <- fitList[[ s ]]$fit
    smooth <- twoPeakCurve(1:365)
    
    temp_larvae <- tibble(
      day = 1:365,
      larva = smooth,
      larva_frac = smooth/sum(smooth),
      site = rep(s,365)
    )
    fit_larvae <- rbind(fit_larvae,temp_larvae)
    
  }
  
  
  fit_larvae <- fit_larvae %>%
    mutate(site = factor(site, levels = c('Foote','Chipman','Snake','Gorge','Chipman2','BRF','SPIN', 'Frost','Gilmore')))    
  
  labv2 <- tibble(elevCat = c('low','mid','high'),
                  julian=rep(212,3),
                  larva=rep(100,3),
                  lab = c('<200 m', '200 - 400 m', '>400 m')) %>%
    mutate(elevCat = factor(elevCat,levels=c('low','mid','high')))
  
  ylab <- expression(paste('Larvae (per 200 ',m^2,')'))
  pdf('figures/l_pheno_fit_bysite.pdf',width=6,height=5)
  samplesMod %>%
    ggplot(aes(julian,larva)) +
    geom_point(cex=0.25) +
    facet_wrap(~site) +
    theme_classic() +
    theme(strip.background = element_rect(color='transparent'),
          strip.text = element_blank(),
          axis.text = element_text(color='black',size = 10),
          axis.title = element_text(size = 10)) +
    geom_line(data=fit_larvae,aes(day,larva)) +
    scale_x_continuous(limits = c(100,300),
                       #                   breaks =c(121, 152, 182, 213, 244, 274),
                       #                  labels=c('May 1', 'Jun 1', 'Jul 1', 'Aug 1','Sep 1', 'Oct 1')) +
                       breaks =c(121,  182, 244),
                       labels=c('May 1', 'Jul 1','Sep 1')) +
    scale_y_continuous(limits = c(0,125)) +
    #geom_text(data=labv2,aes(label=lab),size=4) +
    labs(x='',y=ylab)
  dev.off()
  
}


