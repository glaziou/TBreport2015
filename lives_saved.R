#------------------------------------------------------
# Lives saved
# PG, 22/08/2014
#------------------------------------------------------
setwd('../tbreport2015') 
rm(list=ls())
load('Rdata/global.Rdata')
load('Rdata/regional.Rdata')
source('fun.R')

m <- 1e5
M <- 1e6

lsaved <- function(dta){
  # CFRs untreated
  # HIV negative not on TB treatment  0.43 (0.28 - 0.53), assume ~beta
  cfrn <- 0.43
  cfrn.sd <- (0.53 - 0.28)/4
  # HIV positive not on ART, not on TB treatment  0.78 (0.65 - 0.94), assume ~beta
  cfrp <- 0.78
  cfrp.sd <- (0.94 - 0.65)/4
  
  # lives saved HIV-neg, HIV-pos and total by row
  for (i in 1:dim(dta)[1]){ # will optimize this loop later
    # counterfactual (rates)
    cfn <- betaop(dta$inc.nh[i]/m, cfrn, dta$inc.nh.sd[i]/m, cfrn.sd, op='*', dist=T, nsim=M)
    cfp <- betaop(dta$inc.h[i]/m, cfrp, dta$inc.h.sd[i]/m, cfrp.sd, op='*', dist=T, nsim=M)
    
    # factual (rates)
    parn <- get.beta(dta$mort.nh[i]/m, dta$mort.nh.sd[i]/m)
    fn <- rbeta(M, parn[1], parn[2])
    if (dta$mort.h[i]>0){
      parp <- get.beta(dta$mort.h[i]/m, dta$mort.h.sd[i]/m)
      fp <- rbeta(M, parp[1], parp[2])
    } else fp <- 0
    
    # difference counterfactual minus factual (absolute numbers)
    savedn <- pmax((cfn - fn), 0) * dta$e.pop.num[i]
    savedp <- pmax((cfp - fp), 0) * dta$e.pop.num[i]
    saved <- savedn + savedp
    
    dta$savedn[i] <- mean(savedn) 
    dta$savedn.sd[i] <- sd(savedn) 
    dta$savedn.lo[i] <- quantile(savedn, probs=0.025)
    dta$savedn.hi[i] <- quantile(savedn, probs=0.975) 
    
    dta$savedp[i] <- mean(savedp)
    dta$savedp.sd[i] <- sd(savedp) 
    dta$savedp.lo[i] <- quantile(savedp, probs=0.025) 
    dta$savedp.hi[i] <- quantile(savedp, probs=0.975)
    
    dta$saved[i] <- mean(saved) 
    dta$saved.sd[i] <- sd(saved) 
    dta$saved.lo[i] <- quantile(saved, probs=0.025) 
    dta$saved.hi[i] <- quantile(saved, probs=0.975) 
  }
  
  return(dta)
}



clsaved <- function(dta){
  # cumulative lives saved
  savedn <- sum(dta$savedn)
  savedn.sd <- sqrt(sum(dta$savedn.sd^2))
  savedn.lo <- savedn + qnorm(0.025) * savedn.sd
  savedn.hi <- savedn + qnorm(0.975) * savedn.sd
  
  savedp <- sum(dta$savedp)
  savedp.sd <- sqrt(sum(dta$savedp.sd^2))
  savedp.lo <- savedp + qnorm(0.025) * savedp.sd
  savedp.hi <- savedp + qnorm(0.975) * savedp.sd
  
  saved <- sum(dta$saved)
  saved.sd <- sqrt(sum(dta$saved.sd^2))
  saved.lo <- saved + qnorm(0.025) * saved.sd
  saved.hi <- saved + qnorm(0.975) * saved.sd
  dta <- data.table(savedn, savedn.lo, savedn.hi, savedn.sd,
                       savedp, savedp.lo, savedp.hi, savedp.sd,
                       saved, saved.lo, saved.hi, saved.sd)
  return(dta)
}


# vectorized lohi
vlohi <- Vectorize(lohi, c('ev','se'))

# vectorized betaop
vbetaop <- Vectorize(betaop, c('ev1','ev2','se1','se2'))





#-------------------------------------------------
# globally (2000-2013)
#-------------------------------------------------
set.seed(1234)
global.lsaved <- lsaved(global[year>1999])
global.lsaved <- global.lsaved[, list(year=as.character(year), savedn, savedn.lo, savedn.hi, savedn.sd,
                                      savedp, savedp.lo, savedp.hi, savedp.sd,
                                      saved, saved.lo, saved.hi, saved.sd)]
global.clsaved <- clsaved(global.lsaved)
global.clsaved$year <- 'Cumulative'

saved.global <- rbind(global.lsaved, global.clsaved, use.names=TRUE)

saved.global.print <- saved.global[, list(year, 
                                          saved.hivneg=signif(savedn/M, 3), 
                                          saved.hivneg.ui=paste("(", as.character(signif(savedn.lo/M, 3)),"-", 
                                                              as.character(signif(savedn.hi/M, 3)), ")", sep=''),
                                          saved.hivpos=signif(savedp/M, 3), 
                                          saved.hivpos.ui=paste("(", as.character(signif(savedp.lo/M, 3)),"-", 
                                                         as.character(signif(savedp.hi/M, 3)), ")", sep=''),
                                          saved=signif(saved/M, 3), 
                                          saved.ui=paste("(", as.character(signif(saved.lo/M, 3)),"-", 
                                                         as.character(signif(saved.hi/M, 3)), ")", sep='')
                                          )]
(saved.global.print)

write.csv(saved.global.print, file='output/globalSaved.csv', row.names=FALSE)


p <- qplot(as.numeric(year), savedn/M, data=saved.global[1:14, ], geom='line', colour=I('purple')) +
  geom_ribbon(aes(as.numeric(year), ymax=savedn.hi/M, ymin=savedn.lo/M), fill=I('purple'), alpha=I(0.4)) +
  geom_line(aes(as.numeric(year), savedp/M), colour=I('red')) +
  geom_ribbon(aes(as.numeric(year), ymax=savedp.hi/M, ymin=savedp.lo/M), fill=I('red'), alpha=I(0.4)) +
  xlab('') + ylab('Million') + expand_limits(y=0) +
  theme_bw(base_size=28)
print(p)


#-------------------------------------------------
# by WHO region (2000-2013)
#-------------------------------------------------
set.seed(101)
regional.lsaved <- lsaved(regional[year>1999])
regional.lsaved <- regional.lsaved[, list(year=as.character(year), region=as.character(g.whoregion), 
                                          savedn, savedn.lo, savedn.hi, savedn.sd,
                                          savedp, savedp.lo, savedp.hi, savedp.sd,
                                          saved, saved.lo, saved.hi, saved.sd)]
regional.clsaved <- ddply(regional.lsaved, .(region), function(x)clsaved(x))
regional.clsaved$year <- 'Cumulative'

saved.regional <- rbind(regional.lsaved, regional.clsaved, use.names=TRUE)
saved.regional.print <- saved.regional[, list(year, region, 
                                          saved.hivneg=signif(savedn/M, 3), 
                                          saved.hivneg.ui=paste("(", as.character(signif(savedn.lo/M, 3)),"-", 
                                                                as.character(signif(savedn.hi/M, 3)), ")", sep=''),
                                          saved.hivpos=signif(savedp/M, 3), 
                                          saved.hivpos.ui=paste("(", as.character(signif(savedp.lo/M, 3)),"-", 
                                                                as.character(signif(savedp.hi/M, 3)), ")", sep=''),
                                          saved=signif(saved/M, 3), 
                                          saved.ui=paste("(", as.character(signif(saved.lo/M, 3)),"-", 
                                                         as.character(signif(saved.hi/M, 3)), ")", sep='')
)]
(saved.regional.print)

write.csv(saved.regional, file='output/regionalSaved.csv', row.names=FALSE)


# for global tb report 2015
tab1 <- saved.regional.print[year=='Cumulative'][, year:=NULL]
tab2 <- saved.global.print[.N][, year:=NULL]
tab2$region <- "Global"
tab.saved <- rbind(tab1, tab2, use.names=TRUE)
(tab.saved)
write.csv(tab.saved, file='tab/livesSaved.csv', row.names=FALSE)





#-------------------------------------------------
# country requests
#-------------------------------------------------
# BGD since 1985
load('Rdata/est.Rdata')
cty.lsaved <- function(iso='BGD', start=1995, csv=FALSE){
  th <- 1000  
  bgd <- lsaved(est[iso3==iso & year>=start])
  bgd.cum <- clsaved(bgd)
  bgd.saved <- bgd[, list(year=as.character(year), savedn, savedn.lo, savedn.hi, savedn.sd,
                          savedp, savedp.lo, savedp.hi, savedp.sd,
                          saved, saved.lo, saved.hi, saved.sd)]
  bgd.cum$year <- 'Cumulative'
  bgd.saved <- rbind(bgd.saved, bgd.cum, use.names=T)
  bgd.saved <- bgd.saved[, list(year, 
                                saved.hivneg=signif(savedn/th, 3), 
                                saved.hivneg.ui=paste("(", as.character(signif(savedn.lo/th, 3)),"-", 
                                                      as.character(signif(savedn.hi/th, 3)), ")", sep=''),
                                saved.hivpos=signif(savedp/th, 3), 
                                saved.hivpos.ui=paste("(", as.character(signif(savedp.lo/th, 3)),"-", 
                                                      as.character(signif(savedp.hi/th, 3)), ")", sep=''),
                                saved=signif(saved/th, 3), 
                                saved.ui=paste("(", as.character(signif(saved.lo/th, 3)),"-", 
                                               as.character(signif(saved.hi/th, 3)), ")", sep='')
  )]
  
  if (csv) write.csv(bgd.saved, file=paste('output/',iso,'_livesSaved.csv', sep=''), row.names=FALSE)
  return(bgd.saved)
}

# cty.lsaved(csv=T)



# all countries - request from Mehran, Aug 2015
load('Rdata/est.Rdata')
yr <- 2014
lst <- est[year==yr & inc.nh>0 & inc.h>0 & mort.nh>0, iso3] # current countries with valid data

#---- temp code
sel <- is.na(est$inc.h)
table(sel)
unique(est$iso3[sel])
est$inc.h[sel] <- 1e-5
est$inc.h.sd[sel] <- 1e-5

sel <- is.na(est$inc.nh)
table(sel)
unique(est$iso3[sel])
est$inc.nh[sel] <- 1e-5
est$inc.nh.sd[sel] <- 1e-5

sel <- is.na(est$mort.h)
table(sel)
sel <- est$mort.h==0 
table(sel)
est$mort.h[sel] <- 1e-6
sel <- is.na(est$mort.nh)
table(sel)
sel <- est$mort.nh==0 
table(sel)
est$mort.nh[sel] <- 1e-6
#---- end temp code

# out <- lsaved(est[iso3 %in% lst % year>1999])

# th <- 1e3
# out <- within(out, {
#                    saved.hivneg <- signif(savedn/th, 3) 
#                    saved.hivneg.lo <- signif(savedn.lo/th, 3)
#                    saved.hivneg.hi <- signif(savedn.hi/th, 3)
#                    saved.hivpos <- signif(savedp/th, 3) 
#                    saved.hivpos.lo <- signif(savedp.lo/th, 3) 
#                    saved.hivpos.hi <- signif(savedp.hi/th, 3) 
#                    saved <- signif(saved/th, 3) 
#                    saved.lo <- signif(saved.lo/th, 3) 
#                    saved.hi <- signif(saved.hi/th, 3)

#   }

# cty.clsaved <- ddply(out, .(iso3), function(x)clsaved(x))
# cty.clsaved$year <- 'Cumulative'

# out2 <- rbind(out, cty.clsaved, use.names=TRUE)
# setkey(out2, iso3)

# out3 <- out2[, .(iso3, year, saved.hivneg, saved.hivneg.lo, saved.hivneg.hi, 
#                saved.hivpos, saved.hivpos.lo, saved.hivpos.hi, 
#                saved, saved.lo, saved.hi)]

# write.csv(out3, file='output/lsaved.csv', row.names=F)

# load('Rdata/est.Rdata')

cty.lsaved2 <- function(iso='BGD', start=2000){
  th <- 1000  
  bgd <- lsaved(est[iso3==iso & year>=start])
  bgd.cum <- clsaved(bgd)
  bgd.saved <- bgd[, list(year=as.character(year), savedn, savedn.lo, savedn.hi, savedn.sd,
                          savedp, savedp.lo, savedp.hi, savedp.sd,
                          saved, saved.lo, saved.hi, saved.sd)]
  bgd.cum$year <- 'Cumulative'

  bgd.saved <- rbind(bgd.saved, bgd.cum, use.names=T)
  bgd.saved <- bgd.saved[, list(iso3=iso, year, 
                                saved.hivneg = signif(savedn/th, 3), 
                                saved.hivneg.lo = signif(savedn.lo/th, 3),
                                saved.hivneg.hi = signif(savedn.hi/th, 3),
                                saved.hivpos = signif(savedp/th, 3), 
                                saved.hivpos.lo = signif(savedp.lo/th, 3), 
                                saved.hivpos.hi = signif(savedp.hi/th, 3), 
                                saved = signif(saved/th, 3), 
                                saved.lo = signif(saved.lo/th, 3), 
                                saved.hi = signif(saved.hi/th, 3))]
  
  return(bgd.saved)
}


# lst <- est[year==yr & inc.nh>0 & inc.h>0 & mort.nh>0, iso3] # current countries with valid data
est[iso3 %in% lst, table(year)]  

tmp <- cty.lsaved2(iso=lst[1])
save(tmp, file=paste('misc/AFG.Rdata', sep='')) 

# fail the loop below
fail <- c('MNE', 'SRB')


for (i in setdiff(lst[-1], fail)){ # time to grap a coffee ... or two
  print(i)
  sdl <- est$iso3==i & est$year>1999
  styr <- min(est$year[sdl])
  tmp2 <- cty.lsaved2(iso=i, start=styr)
  save(tmp2, file=paste('misc/',i , '.Rdata', sep='')) 
  tmp <- rbind(tmp, tmp2)
}

# add failed cases
otmp <- tmp
for (i in fail) {
  tmp2 <- cty.lsaved2(iso=i, start=2005)
  tmp <- rbind(tmp, tmp2)
}

setkey(tmp, iso3)

write.csv(tmp, file='output/lsaved.csv', row.names=FALSE)

load('Rdata/est.Rdata')


