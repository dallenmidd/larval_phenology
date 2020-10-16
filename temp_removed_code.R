# some removed code




### fit phenology curve for each site -- going to remove I think
{
  
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
