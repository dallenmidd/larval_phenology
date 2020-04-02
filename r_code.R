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
  samples <- read_csv('data/drag_sampling.csv') %>% mutate(elevCat = cut(elev,c(0,200,410,1000),c('low','mid','high')))
  fitList<-list()
  fitList[['low']] <- list(); fitList[['mid']] <- list(); fitList[['high']] <- list() 
  fitList[['low']][['nymph']] <- list(); fitList[['low']][['larva']] <- list(); fitList[['mid']][['nymph']] <- list(); fitList[['mid']][['larva']] <- list(); fitList[['high']][['nymph']] <- list(); fitList[['high']][['larva']] <- list() 
  
  param_of_interest <- c('peak_e', 'peak_l')
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
    fitList[[which_elev]][['larva']][['fit']] <- fit1
    for (param in param_of_interest) fitList[[which_elev]][['larva']][[param]] <- paramRangeTwoPeak(fit1, param, dataList)
  }
 save(fitList,file = 'data/phenology_fits.RData')
}

# make elevation versus fraction early summer larva plot
{
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
  
  pdf('figures/elev_l_frac.pdf',width=6,height=3)
    larva_comparision %>% 
      ggplot(aes(elev,early_frac)) +
      geom_point() +
      theme_classic() +
      theme(axis.text = element_text(color='black', size =14),
            axis.title = element_text(size = 14)) +
      labs(x = 'Elevation (m)', y = 'Early summer fraction') +
      stat_smooth(se=F,method='lm',color='black') + 
      coord_cartesian(xlim=c(120,600))
  dev.off()
  
}

# make plot phenology plot with fit curves
{
  load(file = 'data/phenology_fits.RData') # see above for code which fit these functions
  samples <- read_csv('data/drag_sampling.csv') %>% mutate(elevCat = cut(elev,c(0,200,410,1000),c('low','mid','high')))
 
   fit1 <- fitList$low$larva$fit
  pred1 <- twoPeakCurve(1:365)
  fit1 <- fitList$mid$larva$fit
  pred2 <- twoPeakCurve(1:365)
  fit1 <- fitList$high$larva$fit
  pred3 <- twoPeakCurve(1:365)
  
  fit_larvae <- tibble(day = rep(1:365,3),
                        larva = c(pred1,pred2,pred3), 
                        larva_frac = c(pred1/sum(pred1),pred2/sum(pred2),pred3/sum(pred3)),
                        elevCat = c(rep('low',365),rep('mid',365),rep('high',365))) %>%
    mutate(elevCat = factor(elevCat,levels=c('low','mid','high')))
  
  labv2 <- tibble(elevCat = c('low','mid','high'),
                  julian=rep(212,3),
                  larva=rep(100,3),
                  lab = c('<200 m', '200 - 400 m', '>400 m')) %>%
    mutate(elevCat = factor(elevCat,levels=c('low','mid','high')))
  ylab <- expression(paste('Larvae (per 200 ',m^2,')'))
  pdf('figures/l_pheno_fit.pdf',width=7,height=3)
    samples %>%
      ggplot(aes(julian,larva)) +
      geom_point(cex=0.25) +
      facet_wrap(~elevCat) +
      theme_classic() +
      theme(strip.background = element_rect(color='transparent'),
            strip.text = element_blank(),
            axis.text = element_text(color='black')) +
      geom_line(data=fit_larvae,aes(day,larva)) +
      scale_x_continuous(limits = c(100,300),
                         #                   breaks =c(121, 152, 182, 213, 244, 274),
                         #                  labels=c('May 1', 'Jun 1', 'Jul 1', 'Aug 1','Sep 1', 'Oct 1')) +
                         breaks =c(121,  182, 244),
                         labels=c('May 1', 'Jul 1','Sep 1')) +
      scale_y_continuous(limits = c(0,100)) +
      labs(x='',y=ylab) +
      geom_text(data=labv2,aes(label=lab),cex=5)
  dev.off()
  
}



