##################################################
# Author:   Philippe Glaziou
#
# GLOBAL TB REPORT figs
##################################################
library(ggplot2)
library(grid)
library(gridExtra)
library(ggthemes)



#--------------------------------------------
# Load datasets -----------------------------
#--------------------------------------------
load ('Rdata/tb.Rdata')
load('Rdata/est.Rdata')
load('Rdata/global.Rdata')
load('Rdata/regional.Rdata')
load('Rdata/hbc.Rdata')
load('Rdata/regional_ff.Rdata')
load('Rdata/global_ff.Rdata')
load('Rdata/hbc_ff.Rdata')

m <- 100000
source('fun.R') # we need multiplot()



# function to generate low and high bounds assuming Beta distribution
lohi <- function(ev, se){
  par <- get.beta(ev, se)
  lo <- qbeta(0.025, par[1], par[2])
  hi <- qbeta(0.975, par[1], par[2])
  return(c(lo=lo, hi=hi))
}

# vectorized 
vlohi <- Vectorize(lohi, c('ev','se'))




#----------------------------------------------------------
# Global rates of incidence, prevalence and mortality -----
#----------------------------------------------------------
# global plot
pdf(width=8, height=4, file='fig/fig2_6_global.pdf')
p1 <- qplot(year, inc, data=global, geom='line', colour=I('blue')) +
    geom_ribbon(aes(year, ymin=inc.lo, ymax=inc.hi), fill=I('blue'), alpha=0.4) +
    geom_line(aes(year, inc.h), colour=I('red')) +
    geom_line(aes(year, newinc)) + 
    geom_ribbon(aes(year, ymin=inc.h.lo, ymax=inc.h.hi), 
      fill=I('red'), alpha=0.4) +
    ylab('Rate per 100,000 population/year') + xlab('') +
    expand_limits(y=0) +
    theme_bw(base_size=10) 

p2 <- qplot(year, mort.nh, data=global.ff, geom='line', colour=I('blue'), linetype=forecast) +
    geom_ribbon(aes(year, ymin=mort.nh.lo, ymax=mort.nh.hi), fill=I('blue'), alpha=0.4) +
    geom_hline(aes(yintercept=mort.nh[1] / 2), linetype=2) +
    ylab('Mortality Rate per 100,000 population/year') + xlab('') +
    expand_limits(y=0) +
    theme_bw(base_size=10) + theme(legend.position='none') 

p3 <- qplot(year, prev, data=global.ff, geom='line', colour=I('blue'), linetype=forecast) +
    geom_ribbon(aes(year, ymin=prev.lo, ymax=prev.hi), fill=I('blue'), alpha=0.4) +
    geom_hline(aes(yintercept=prev[1] / 2), linetype=2) +
    ylab('Prevalence rate per 100,000 population/year') + xlab('') +
    expand_limits(y=0) +
    theme_bw(base_size=10) + theme(legend.position='none') 
 
multiplot(p1, p3, p2, cols=3)
dev.off()



pdf(width=8, height=8, file='fig/fig4_1_globalNotificationRates.pdf')
p1 <- qplot(year, inc, data=global, geom='line', colour=I('darkgreen')) +
  geom_ribbon(aes(year, ymin=inc.lo, ymax=inc.hi), fill=I('darkgreen'), alpha=0.4) +
  geom_line(aes(year, newinc)) +   ylab('Rate per 100,000 population/year') + xlab('') +
  expand_limits(y=0) +
  theme_bw(base_size=18) 
print(p1)
dev.off()




#----------------------------------------------------
# Global numbers (incidence and mortality) ----------
#----------------------------------------------------
library(gridExtra)
pdf(width=8, height=4, file='fig/fig2_2_global_num.pdf', 
    title = paste0("Figure 2.2 Estimated absolute numbers of TB cases and deaths (in millions per year), 1990–", max(global$year)))
mil <- 1e6
p1 <- qplot(year, inc.num/mil, data=global, geom='line', colour=I('black')) +
  geom_ribbon(aes(year, ymin=inc.lo.num/mil, ymax=inc.hi.num/mil), fill=I('grey'), alpha=0.8) +
  geom_line(aes(year, inc.h.num/mil), colour=I('red')) +
  # geom_line(aes(year, c.newinc/mil)) + 
  geom_ribbon(aes(year, ymin=inc.h.lo.num/mil, ymax=inc.h.hi.num/mil), 
              fill=I('red'), alpha=0.4) +
                ylab('Millions') + xlab('') +
                expand_limits(y=0) +
                theme_bw(base_size=11) +
  ggtitle("TB incidence")

p2 <- qplot(year, mort.nh.num/mil, data=global, geom='line', colour=I('blue')) +
  geom_ribbon(aes(year, ymin=mort.nh.lo.num/mil, ymax=mort.nh.hi.num/mil), fill=I('blue'), alpha=0.3) +
  geom_line(aes(year,mort.h.num/mil), colour=I('red')) +
  geom_ribbon(aes(year, ymin=mort.h.lo.num/mil, ymax=mort.h.hi.num/mil), fill=I('red'), alpha=0.3) +
  # geom_line(aes(year, mort.num/mil), colour=I('black')) +
  # geom_ribbon(aes(year, ymin=mort.lo.num/mil, ymax=mort.hi.num/mil), fill=I('grey'), alpha=0.8) +
  ylab('Millions') + xlab('') +
  expand_limits(y=0) +
  theme_bw(base_size=11) +
  theme(legend.position='none') +
  ggtitle("TB deaths")

grid.arrange(p1, p2, ncol=2, nrow=1, 
             main=paste0("Figure 2.2 Estimated absolute numbers of TB cases and deaths (in millions per year), 1990–", max(global$year)))
# multiplot(p1, p2, cols=2)
dev.off()









#----------------------------------------------------
# Regional notification rate plus inc ---------------
#----------------------------------------------------
# chapter 2
pdf(width=8, height=8, file='fig/fig2_7_incidenceByRegion.pdf')
p <- qplot(year, inc, data=as.data.frame(regional), geom='line', colour=I('blue')) +
  geom_ribbon(aes(year, ymin=inc.lo, ymax=inc.hi), fill=I('blue'), alpha=0.4) +
  geom_line(aes(year, inc.h), colour=I('red')) +
  geom_line(aes(year, newinc)) + 
  geom_ribbon(aes(year, ymin=inc.h.lo, ymax=inc.h.hi), 
              fill=I('red'), alpha=0.4) +
  facet_wrap(~g.whoregion, scales='free_y') +
  xlab('') + ylab('Rate per 100,000 population/year') +
  expand_limits(y=0) +
  theme_bw(base_size=12) 
print(p)
dev.off()



# chapter 3
pdf(width=8, height=8, file='fig/fig4_2_notifByRegion.pdf')
p <- qplot(year, inc, data=as.data.frame(regional), geom='line', colour=I('darkgreen')) +
    geom_ribbon(aes(year, ymin=inc.lo, ymax=inc.hi), fill=I('darkgreen'), alpha=0.4) +
    geom_line(aes(year, newinc)) + 
    facet_wrap(~g.whoregion, scales='free_y') +
    xlab('') + ylab('Rate per 100,000 population/year') +
    expand_limits(y=0) +
    theme_bw(base_size=12) 
print(p)
dev.off()









#------------------------------------------------------
# Prevalence and Mortality by WHO region --------------
#------------------------------------------------------
detach(package:dplyr)
library(plyr); library(dplyr)
regional.ff2 <- ddply(as.data.frame(regional.ff), .(g.whoregion), transform,
  target.prev=prev[1]/2, target.mort=mort.nh[1]/2, target.morth=mort.h[15]/2)

pdf(width=8, height=8, file='fig/fig2_10_prevalenceByRegion.pdf')
p <- qplot(year, prev, data=regional.ff2, geom='line', colour=I('blue')) +
    geom_ribbon(aes(year, ymin=prev.lo, ymax=prev.hi), fill=I('blue'), alpha=0.4) +
    geom_hline(aes(yintercept=target.prev), linetype=2) +
    facet_wrap(~g.whoregion, scales='free_y') +
    xlab('') + ylab('Rate per 100,000 population/year') +
    expand_limits(y=0) +
    theme_bw(base_size=12) +
    theme(legend.position='none')
print(p)
dev.off()

pdf(width=8, height=8, file='fig/fig2_13_mortalityByRegion.pdf')
p <- qplot(year, mort.nh, data=regional.ff2, geom='line', colour=I('blue')) +
    geom_ribbon(aes(year, ymin=mort.nh.lo, ymax=mort.nh.hi), fill=I('blue'), alpha=0.4) +
    geom_hline(aes(yintercept=target.mort), linetype=2) +
    facet_wrap(~g.whoregion, scales='free_y') +
    xlab('') + ylab('Rate per 100,000 population/year') +
    expand_limits(y=0) +
    theme_bw(base_size=12) +
    theme(legend.position='none')
print(p)
dev.off()


pdf(width=8, height=8, file='fig/fig2_13b_mortalityHIVposByRegion.pdf')
p <- qplot(year, mort.h, data=subset(regional.ff2, year<2013), geom='line', colour=I('blue')) +
  geom_ribbon(aes(year, ymin=mort.h.lo, ymax=mort.h.hi), fill=I('blue'), alpha=0.4) +
  geom_hline(aes(yintercept=target.morth), linetype=2) +
  facet_wrap(~g.whoregion, scales='free_y') +
  xlab('') + ylab('Rate per 100,000 population/year') +
  expand_limits(y=0) +
  theme_bw(base_size=12) +
  theme(legend.position='none')
print(p)
dev.off()






#----------------------------------------------------
# AFR trends in AIDS-TB mortality -------------------
#----------------------------------------------------
pdf(width=8, height=8, file='fig/fig_afr_mortality.pdf')
p <- qplot(year, mort.h, data=regional[g.whoregion=='AFR'], geom='line', colour=I('red')) +
  geom_ribbon(aes(year, ymin=mort.h.lo, ymax=mort.h.hi), fill=I('red'), alpha=0.4) +
  geom_line(aes(year, mort), colour=I('grey')) +
  geom_ribbon(aes(year, ymin=mort.lo, ymax=mort.hi), fill=I('grey'), alpha=0.7) +  
  geom_line(aes(year, mort.nh), colour=I('blue')) +
  geom_ribbon(aes(year, ymin=mort.nh.lo, ymax=mort.nh.hi), fill=I('blue'), alpha=0.2) +  
  xlab('') + ylab('Rate per 100,000 population/year') +
  expand_limits(y=0) +
  theme(legend.position='none') +
  theme_bw(base_size=18)
print(p)
dev.off()






#----------------------------------------------------
# Incidence in 22 HBCs ------------------------------
#----------------------------------------------------
hest <- subset(est, g.hbc22=='high' & !is.na(g.hbc22))
levels(hest$country)[match('Democratic Republic of the Congo', levels(hest$country))] <- 'DR Congo'
levels(hest$country)[match('United Republic of Tanzania', levels(hest$country))] <- 'UR Tanzania'

hest <- within(hest, {
  inc.h.lo <- vlohi(inc.h/m, inc.h.sd/m)[1, ]*m
  inc.h.hi <- vlohi(inc.h/m, inc.h.sd/m)[2, ]*m
})

#pdf(width=8, height=8, file='fig/fig2_8_incidence_22hbc.pdf')
p1 <- qplot(year, inc, data=as.data.frame(hest), geom='line', colour=I('blue')) +
    geom_ribbon(aes(year, ymin=inc.lo, ymax=inc.hi), fill=I('blue'), alpha=0.4) +
    geom_line(aes(year, inc.h), colour=I('red')) +
    geom_line(aes(year, newinc)) + 
    geom_ribbon(aes(year, ymin=inc.h.lo, ymax=inc.h.hi), 
      fill=I('red'), alpha=0.4) +
    facet_wrap(~country, scales='free_y') +
    xlab('') + ylab('Rate per 100,000 population/year') +
    expand_limits(y=0) +
    theme_bw(base_size=8) 
print(facetAdjust(p1))
ggsave('fig/fig2_8_incidence_22hbc.pdf', width=10, height=8)
#dev.off()

# notifs and incidence (chapter 3)
#pdf(width=8, height=8, file='fig/fig3_3_incidenceRates_22hbc.pdf')
p2 <- qplot(year, inc, data=as.data.frame(hest), geom='line', colour=I('blue')) +
  geom_ribbon(aes(year, ymin=inc.lo, ymax=inc.hi), fill=I('blue'), alpha=0.4) +
  geom_line(aes(year, newinc)) + 
  facet_wrap(~country, scales='free_y') +
  xlab('') + ylab('Rate per 100,000 population/year') +
  expand_limits(y=0) +
  theme_bw(base_size=8) 
print(facetAdjust(p2))
ggsave('fig/fig3_3_incidenceRates_22hbc.pdf', width=10, height=8)
#dev.off()






#----------------------------------------------------
# Incidence in top 10 countries ---------------------
#----------------------------------------------------
(top1 <- head(est[year==yr][order(inc.num, decreasing=T), list(country, inc.num)], 10))
levels(top1$country)[match('Democratic Republic of the Congo', levels(top1$country))] <- 'DR Congo'

pdf(width=8, height=8, file='fig/fig2_3_inc_top10inc.pdf')
p1 <- qplot(year, inc, data=est[country %in% as.character(top1$country)], geom='line', colour=I('blue')) +
  geom_ribbon(aes(year, ymin=inc.lo, ymax=inc.hi), fill=I('blue'), alpha=0.4) +
  geom_line(aes(year, newinc)) + 
  facet_wrap(~country, scales='free_y') +
  xlab('') + ylab('Rate per 100,000 population/year') +
  expand_limits(y=0) +
  theme_bw(base_size=10) 
print(p1)
#print(facetAdjust(p2))
dev.off()








#------------------------------------------------------
# Prevalence and Mortality in 22 HBCs -----------------
#------------------------------------------------------
hbc.ff2 <- ddply(as.data.frame(hbc.ff), .(iso3), transform,
  target.prev=prev[1]/2, target.mort=mort.nh[1]/2)

dim(hbc.ff2)
load('Rdata/cty.Rdata')
hbc.ff2 <- merge(hbc.ff2, cty[, list(iso3, country)], by='iso3', all.x=TRUE, all.y=FALSE)
dim(hbc.ff2)

levels(hbc.ff2$country)[match('Democratic Republic of the Congo', levels(hbc.ff2$country))] <- 'DR Congo'
levels(hbc.ff2$country)[match('United Republic of Tanzania', levels(hbc.ff2$country))] <- 'UR Tanzania'


#pdf(width=8, height=8, file='fig/fig2_14_prevalence_22hbc.pdf')
p1 <- qplot(year, prev, data=hbc.ff2, geom='line', colour=I('blue')) +
    geom_ribbon(aes(year, ymin=prev.lo, ymax=prev.hi), fill=I('blue'), alpha=0.4) +
    geom_hline(aes(yintercept=target.prev), linetype=2) +
    facet_wrap(~country, scales='free_y') +
    xlab('') + ylab('Rate per 100,000 population/year') +
    expand_limits(y=0) +
    theme_bw(base_size=8) +
    theme(legend.position='none')
print.facetAdjust(facetAdjust(p1))
ggsave(file='fig/fig2_14_prevalence_22hbc.pdf', width=10, height=8)
#dev.off()




#pdf(width=8, height=8, file='fig/fig2_14_mortalityRates_hbc22.pdf')
hbc.ff3 <- merge(hbc.ff2, est[, list(iso3,year,vr.raw)], by=c('iso3','year'), all.x=TRUE, all.y=FALSE)
p2 <- qplot(year, mort.nh, data=hbc.ff3, geom='line', colour=I('blue'), linetype=forecast) +
    geom_ribbon(aes(year, ymin=mort.nh.lo, ymax=mort.nh.hi), fill=I('blue'), alpha=0.4) +
    geom_hline(aes(yintercept=target.mort), linetype=2) +
    geom_point(aes(year, vr.raw), shape=I(4)) +
    facet_wrap(~country, scales='free_y') +
    xlab('') + ylab('Rate per 100,000 population/year') +
    expand_limits(y=0) +
    theme_bw(base_size=8) +
    theme(legend.position='none')
print(facetAdjust(p2))
ggsave(file='fig/fig2_14_mortalityRates_hbc22.pdf', width=10, height=8)
#dev.off()
#ggsave(file='fig/fig2_14b_mortality_22hbc.pdf', width=8, height=8)


#pdf(width=8, height=8, file='fig/fig2_14c_mortalityHIVpos_hbc22.pdf')
p2 <- qplot(year, mort.h, data=subset(est, g.hbc22=='high'), geom='line', colour=I('blue')) +
  geom_ribbon(aes(year, ymin=mort.h.lo, ymax=mort.h.hi), fill=I('blue'), alpha=0.4) +
  facet_wrap(~country, scales='free_y') +
  xlab('') + ylab('Rate per 100,000 population/year') +
  expand_limits(y=0) +
  theme_bw(base_size=8)
print(facetAdjust(p2))
ggsave(file='fig/fig2_14c_mortalityHIVpos_22hbc.pdf', width=10, height=8)
#dev.off()



# Indonesia plots ---------------------------------------------------
pdf(width=8, height=5, file='fig/IDNbox.pdf')
p1 <- qplot(year, inc, data=est['IDN'], geom='line', colour=I('blue')) +
  geom_ribbon(aes(year, ymin=inc.lo, ymax=inc.hi), fill=I('blue'), alpha=0.4) +
  geom_line(aes(year, newinc)) + 
  xlab('') + ylab('Rate per 100,000 population/year') +
  ggtitle('Incidence') +
  expand_limits(y=0) +
  theme_bw(base_size=10) 
p2 <- qplot(year, mort.nh, data=subset(hbc.ff3, iso3=='IDN'), 
            geom='line', colour=I('blue')) +
    geom_ribbon(aes(year, ymin=mort.nh.lo, ymax=mort.nh.hi), fill=I('blue'), alpha=0.4) +
    geom_hline(aes(yintercept=target.mort), linetype=2) +
    geom_point(aes(year, vr.raw), shape=I(4)) +
    xlab('') + ylab('Rate per 100,000 population/year') +
  ggtitle('Mortality (HIV-neg)') +
    expand_limits(y=0) +
    theme_bw(base_size=10) +
    theme(legend.position='none')
p3 <- qplot(year, prev, data=subset(hbc.ff2, iso3=='IDN'), 
            geom='line', colour=I('blue')) +
    geom_ribbon(aes(year, ymin=prev.lo, ymax=prev.hi), fill=I('blue'), alpha=0.4) +
    geom_hline(aes(yintercept=target.prev), linetype=2) +
    xlab('') + ylab('Rate per 100,000 population') +
    ggtitle('Prevalence') +
    expand_limits(y=0) +
    theme_bw(base_size=10) +
    theme(legend.position='none')
multiplot(p1, p3, p2, cols=3)
dev.off()




# India plots ---------------------------------------------------
pdf(width=8, height=5, file='fig/INDbox.pdf')
p1 <- qplot(year, inc, data=est['IND'], geom='line', colour=I('blue')) +
  geom_ribbon(aes(year, ymin=inc.lo, ymax=inc.hi), fill=I('blue'), alpha=0.4) +
  geom_line(aes(year, newinc)) + 
  xlab('') + ylab('Rate per 100,000 population/year') +
  ggtitle('Incidence') +
  expand_limits(y=0) +
  theme_bw(base_size=10) 
p2 <- qplot(year, mort.nh, data=subset(hbc.ff3, iso3=='IND'), 
            geom='line', colour=I('blue')) +
    geom_ribbon(aes(year, ymin=mort.nh.lo, ymax=mort.nh.hi), fill=I('blue'), alpha=0.4) +
    geom_hline(aes(yintercept=target.mort), linetype=2) +
    geom_point(aes(year, vr.raw), shape=I(4)) +
    xlab('') + ylab('Rate per 100,000 population/year') +
  ggtitle('Mortality (HIV-neg)') +
    expand_limits(y=0) +
    theme_bw(base_size=10) +
    theme(legend.position='none')
p3 <- qplot(year, prev, data=subset(hbc.ff2, iso3=='IND'), 
            geom='line', colour=I('blue')) +
    geom_ribbon(aes(year, ymin=prev.lo, ymax=prev.hi), fill=I('blue'), alpha=0.4) +
    geom_hline(aes(yintercept=target.prev), linetype=2) +
    xlab('') + ylab('Rate per 100,000 population') +
    ggtitle('Prevalence') +
    expand_limits(y=0) +
    theme_bw(base_size=10) +
    theme(legend.position='none')
multiplot(p1, p3, p2, cols=3)
dev.off()





#--------------------------------------------------------------------
# top 10 causes of deaths
#--------------------------------------------------------------------
library(grid)

# 2012 top 10 
cod <- data.table(cause=c('Ischaemic heart disease','Stroke','COPD','Lower respiratory infections',
                          'Tracheal, bronchus, lung cancer','HIV/AIDS','TB','Diarrheal diseases',
                          'Diabetes mellitus','Road injury'),
                  ncd=c(T,T,T,F,F,F,F,F,T,T),
#                  n=c(7.4, 6.7, 3.1, 3.1, 1.6, 1.5, 1.3, 1.5, 1.5, 1.3),
#                  tbhiv=c(rep(0, 5), rep(0.32, 2), rep(F,3)))
                  n=c(7.4, 6.7, 3.1, 3.1, 1.6, 0.98, 1.18, 1.5, 1.5, 1.3),
                  tbhiv=c(rep(0, 5), rep(0.422545, 2), rep(F,3)))
cod[, total := n + tbhiv]

# 2014 TB and HIV only
cod2 <- data.table(cause=c('HIV/AIDS','TB'),
                  n=c(0.794683, 1.113),
                  tbhiv=rep(0.387010, 2))
cod2[, total := n + tbhiv]

p <- qplot(reorder(cause, total), total, data=cod, geom='bar', stat='identity', 
           fill=I('grey50'), colour=I('black'),
           main='Top 10 causes of deaths, with TB/HIV deaths shown in grey.\n') +
    geom_bar(aes(cause, n), stat='identity', fill=I('white'), colour=I('black')) +
    xlab('') + ylab('Millions (2012)') +
    scale_y_continuous(breaks=0:7) +
    coord_flip()

inset <- qplot(reorder(cause, total), total, data=cod2, geom='bar', stat='identity', 
           fill=I('grey50'), colour=I('black')) +
    geom_bar(aes(cause, n), stat='identity', fill=I('white'), colour=I('black')) +
    xlab('') + ylab('(2014)') +
    coord_flip()

vp <- viewport(width=0.3, height=0.3, x=0.8, y=0.25)

pdf(file='fig/top10C0D.pdf', width=8, height=6)
print(p)
#print(inset, vp = vp)
dev.off()


p2 <- qplot(reorder(cause, total), total, data=cod2, geom='bar', stat='identity', 
           fill=I('grey50'), colour=I('black')) +
    geom_bar(aes(cause, n), stat='identity', fill=I('white'), colour=I('black')) +
    xlab('') + ylab('Millions (2014)') +
    coord_flip()

#pdf(file='fig/TBHIVdeaths2014.pdf', width=8, height=6)
print(p2)
#print(inset, vp = vp)
#dev.off()


p2 <- qplot(reorder(cause, total), total, data=cod[ncd==F], geom='bar', stat='identity', 
           fill=I('grey50'), colour=I('black'),
           main='Top 5 communicable causes of deaths, TB/HIV deaths in grey.\n') +
    geom_bar(aes(cause, n), stat='identity', fill=I('white'), colour=I('black')) +
    xlab('') + ylab('Millions (2012)') +
    scale_y_continuous(breaks=0:7) +
    coord_flip()

pdf(file='fig/topCommunicableC0D.pdf', width=8, height=6)
print(p2)
print(inset, vp = vp)
dev.off()

# 2012 top 10 infectious
codb <- data.table(cause=c('Respiratory infections','Diarrhoeal diseases','Malaria','Meningitis, encephalitis',
                          'Vaccine preventable','HIV/AIDS','TB','Hepatitis B and C',
                          'Neglected tropical diseases','STIs, excluding HIV'),
                  n=c(3.61, 1.5, 0.62, 0.47, 0.27, 0.98, 1.18, 0.188, 0.171, 0.084),
                  tbhiv=c(rep(0, 5), rep(0.422545, 2), rep(F,3)))
codb[, total := n + tbhiv]

pb <- qplot(reorder(cause, total), total, data=codb, geom='bar', stat='identity', 
           fill=I('grey50'), colour=I('black'),
           main='Top 10 infectious causes of deaths, TB/HIV deaths in grey.\n') +
    geom_bar(aes(cause, n), stat='identity', fill=I('white'), colour=I('black')) +
    xlab('') + ylab('Millions (2012)') +
    scale_y_continuous(breaks=0:7) +
    coord_flip()

pdf(file='fig/top10infectiousC0D.pdf', width=8, height=6)
print(pb)
print(inset, vp = vp)
dev.off()




#----------------------------------------------------
# Maps
#----------------------------------------------------
library(whomap)

# epi regions
dta <- merge(est[year==2014,.(iso3)], cty[, .(iso3, g.est)], by='iso3')
dta$g.est[dta$iso3=='SSD'] <- 'AFRhigh'
dta$var <- factor(dta$g.est, labels=c('Africa, high-HIV', 'Africa, low-HIV',
                                      'Central Europe','East Europe',
                                      'High-income','East Mediterranean',
                                      'Latin America','South East Asia', 'West Pacific'))

whomap(X=dta) + scale_fill_brewer("Epi regions", palette='Set1')
ggsave(file='fig/epiregions.pdf', width=10, height=8)


# prev
yr <- 2014
dta <- est[year==2014, list(iso3, g.hbc22, newinc, snewinc, inc, inc.num, tbhiv, mort, mort.nh, prev, 
                            mort.num, source.inc, source.mort, source.prev)]
dta$var <- dta$source.prev
dta$var[is.na(dta$var)] <- 'indirect'
dta$var[dta$iso3 %in% c('TUR', 'PRK', 'KAZ', 'SSD')] <- 'indirect'
whomap(X=dta, Z=scale_fill_manual("survey", values=c('grey','darkgreen')), legendpos='none')
ggsave(file='fig/PScountries.pdf', width=10, height=8)

# VR countries 
dta$var <- dta$source.mort=='imputed' 
whomap(X=dta, Z=scale_fill_manual("VR/Survey", values=c('grey','darkgreen')), legendpos='none')

ggsave(file='fig/fig2_11_map_VRcountries.pdf', width=10, height=8)
write.csv(dta, file='tab/vrCountries.csv', row.names=FALSE)


# source Incidence (2014)
dta <- est[year==yr, list(iso3, g.hbc22, newinc, snewinc, inc, inc.num, tbhiv, mort, mort.nh, prev, 
                            mort.num, source.inc, source.mort, source.prev)]
#dta$var <- dta$source.inc %in% c("Capture-recapture","High income","Survey")
dta$var <- dta$source.inc
dta$var[dta$iso3 %in% c('CHN', 'GMB', 'IDN', 'MMR', 'PAK', 'PHL', 'RWA', 'VNM')] <- 'Prevalence survey'
dta$var[dta$var=='Survey'] <- 'Prevalence survey'
dta$var[dta$var %ni% c('Prevalence survey', 'Capture-recapture', 'High income')] <- 'Case notifications'
dta$var[dta$iso3 %in% c('FRA', 'RUS')]  <- 'High income'
dta$var[dta$iso3 %in% c('EGY', 'NLD')]  <- 'Capture-recapture'
dta$var <- factor(dta$var, levels=c('Case notifications','Prevalence survey',
                                    'High income','Capture-recapture'),
                  labels=c('Case notifications,\nexpert opinion', 'Prevalence survey',
                           'Case notifications,\nstandard adjustment','Capture-recapture'))

p <- whomap(X=dta) + scale_fill_brewer('Main method', palette='Set1') 

title.grob <- textGrob(
    label = '\nFigure 2.2. Main method used to estimate TB incidence\n
              In the first method, case notification data are combined with expert opinion about 
              case detection gaps (under-reporting and under-diagnosis), and trends are estimated 
              using either mortality data, surveys of the annual risk of infection or exponential 
              interpolation using estimates of case detection gaps for three years. For all 
              high-income countries except the Netherlands and the United Kingdom, notifications 
              are adjusted by a standard amount to account for case detection gaps. 
              For further details about all four methods, see text.',
    x = unit(0, "lines"), 
    y = unit(0, "lines"),
    hjust = 0, vjust = 0,
    gp = gpar(fontsize = 16))

p1 <- arrangeGrob(p, main = title.grob)

pdf(file='fig/fig2_1_map_incmethod.pdf', width=10, height=8)
(p1)
dev.off()

pdf(file='fig/fig2_1_map_incmethod_clean.pdf', width=10, height=8)
(p)
dev.off()

write.csv(dta[, list(iso3, source.inc=var)], file='tab/inc_method.csv', row.names=F)

table(dta$var)
print(dta[, .(sum(inc.num)), by=var][, prop:=V1/sum(V1)])




# country consultations
load('Rdata/cty.Rdata')
emr <- cty$iso3[cty$g.whoregion=='EMR']
emr <- setdiff(emr, c('SSD','LBY','PSE'))
sea <- cty$iso3[cty$g.whoregion=='SEA']
sea <- setdiff(sea, c('PRK'))
wpr <- c('KHM','CHN','LAO','MYS','MNG','PNG','PHL','VNM')
amr <- c('BOL','BRA','COL','DOM','ECU','SLV','GTM','GUY','HTI','HND','MEX','NIC','PRY','PER','URY','VEN')
afr <- c('BWA','BFA','BDI','CIV','COD','ETH','GHA','KEN','MWI','MLI','MOZ','NAM','NGA','RWA','UGA','ZMB','ZWE')
eur <- cty$iso3[cty$g.est=='EEUR' & !is.na(cty$g.est)]
eur <- setdiff(eur, 'ARM')
eur <- c(eur, 'ALB','BGR','EST','LVA','LTU','MNE','MKD','RUS','SRB','TUR')
tot <- c(amr, emr, sea, afr, wpr, eur)
dta$var <- dta$iso3 %in% tot
whomap(X=dta, Z=scale_fill_manual("Surveillance", values=c('grey','darkgreen')), legendpos='none')

ggsave(file='fig/fig2_5_map_ctyconsult.pdf', width=10, height=8)


# incidence
dta$var <- cut(dta$inc, c(0, 25, 50, 100, 200, 300, Inf), 
               c('0-24.9', '25-49.9', '50-99', '100-199.9', '200-299', '300+'), 
               right=F, ordered_result=T)
whomap(X=dta, Z=scale_fill_brewer("Incidence\nper 100,000", palette="YlGnBu", type="seq"))
ggsave(file='fig/fig2_5_map_incidenceRates.pdf', width=10, height=8)



# mortality HIV-
dta$var <- cut(dta$mort.nh, c(0, 1, 2, 5, 10, 20, 40, Inf), 
               c('0-0.9', '1-1.9', '2-4.9', '5-9.9', '10-19','20-39','40+'), 
               right=F, ordered_result=T)
whomap(X=dta, Z=scale_fill_brewer("Mortality\nper 100,000", palette="YlGnBu", type="seq"))
ggsave(file='fig/fig2_12_map_mortalityRatesHIVneg.pdf', width=10, height=8)



# prevalence
# dta$var <- cut(dta$prev, c(0, 5, 10, 20, 50, 100, 200, Inf), 
#                c('0-4.9', '5-9.9', '10-19', '20-49', '50-99', '100-199','200+'), 
#                right=F, ordered_result=T)
# whomap(X=dta, Z=scale_fill_brewer("Prevalence\nper 100,000", palette="YlGnBu", type="seq"))
# ggsave(file='fig/fig2_4_map_prevalence.pdf')


# HIV prevalence in TB
dta$var <- cut(dta$tbhiv*100, c(0, 1, 2, 5, 10, 20, 50, Inf), 
               c('0-0.9', '1-1.9', '2-4.9', '5-9.9', '10-19', '20-49','50+'), 
               right=F, ordered_result=T)
whomap(X=dta, Z=scale_fill_brewer("HIV Prevalence\namong TB (%)", palette="YlGnBu", type="seq"))
ggsave(file='fig/fig2_4_map_HIVamongTB.pdf', width=10, height=8)



# Prevalence surveys
# Prevalence survey countries (2013)
pr.sub <- 'IND'
pr.one <- c('MMR','PHL','BGD','ERI','ETH','VNM','PAK','LAO','NGA','GMB','RWA','GHA','SDN')
pr.ongoing <- c('TZA','THA','MNG','ZWE')
pr.late <- c('KEN','ZAF','MWI','UGA','ZMB','BGD','PRK','NEP','MOZ')
pr.torepeat <- c('VNM','BGD','MMR','ERI','IDN')
pr.repeated <- c('KOR','CHN','KHM','PHL')

dta$var <- NULL
dta$var <- 'Not planned'
dta$var[dta$iso3 %in% c(pr.ongoing, pr.late)] <- 'Ongoing/Planned'
dta$var[dta$iso3 %in% pr.sub] <- 'Sub-national'
dta$var[dta$iso3 %in% pr.one] <- 'One completed'
dta$var[dta$iso3 %in% pr.repeated] <- 'Repeated'
dta$var[dta$iso3 %in% pr.torepeat] <- 'Repeat planned'
dta$var <- factor(dta$var, 
                  levels=c('Not planned','Sub-national','Ongoing/Planned','One completed','Repeat planned','Repeated'))
(table(dta$var))

whomap(X=dta, Z=scale_fill_manual("Prevalence survey", values=c('grey','lightblue','blue','yellow4','brown','green4')), legendpos=c(0.14, 0.3))
ggsave(file='fig/fig2_9_map_prevalence_surveys.pdf', width=10, height=8)




# change in prevalence estimates ------------------------------------
library(reshape2)
setwd('../tbreport2015')
post <- read.csv('input/prev_prepost.csv')
load('Rdata/ws.Rdata') 
vlohi <- Vectorize(lohi, c('ev', 'se'))
post2 <- ws[!is.na(uprev), .(iso3, country=as.character(country), year, prev=uprev, prev.se=uprev.sd,
                             old.prev=prev, old.prev.se=prev.sd)]
post3 <- rbind(post, post2)
m <- 1e5
post <- within(post3, {
  prev.hi <- vlohi(prev/m, prev.se/m)[2, ] * 1e3
  prev.lo <- vlohi(prev/m, prev.se/m)[1, ] * 1e3
  old.prev.hi <- vlohi(old.prev/m, old.prev.se/m)[2, ] * 1e3
  old.prev.lo <- vlohi(old.prev/m, old.prev.se/m)[1, ] * 1e3
})
post <- within(post, {
  prev <- prev/100
  old.prev <- old.prev/100
})
(post)

# difference
post <- within(post, {
               change <- prev - old.prev
               change.se <- sqrt(prev.se^2 + old.prev.se^2)
})
post <- within(post, {
               change.lo <- change - 1.96*change.se
               change.hi <- change + 1.96*change.se
})

# p1 <- qplot(change, country, data=post, geom='point', size=I(1.1), colour=I('red'), fill=I('red')) +
#  geom_segment(aes(x=change.lo, xend=change.hi, y=country, yend=country), colour=I('red')) +
#  xlab('Difference survey - prior (per 1000)') + ylab('')
# print(p1)

load('Rdata/cty.Rdata')
post2 <- merge(post, cty[, .(iso3, g.whoregion)], by='iso3', all.x=T, all.y=F)
post2$region <- 'Africa'
post2$region[post2$g.whoregion %in% c('WPR','SEA')]  <- 'Asia'
post2$region[post2$iso3 %in% c('Sudan')]  <- 'Africa'
post2$region[post2$iso3 %in% c('PAK')]  <- 'Asia'
post2$region <- factor(post2$region, levels=c('Asia', 'Africa'))

p2 <- qplot(old.prev, reorder(country, change), data=post2[iso3 %ni% c('PHL', 'VNM')], geom='point', size=I(3), colour=I('lightblue'), 
            fill=I('lightblue')) +
  geom_segment(aes(x=old.prev.lo, xend=old.prev.hi, y=country, yend=country), colour=I('lightblue'),
               size=I(3)) +
  geom_point(aes(prev, country), colour=I('red'), size=I(1)) +
  geom_segment(aes(x=prev.lo, xend=prev.hi, y=country, yend=country), 
               colour=I('red'), size=I(1)) +
  facet_wrap(~region, nrow=1, scales='free_y') +
  scale_x_log10(breaks=c(.25, 0.5, 1, 2, 5, 10)) +
  xlab('Prevalence per 1000 population (log scale)') + ylab('') +
  theme_bw(base_size=18)

title.grob <- textGrob(
    label = '\nFigure B2.2.1: Estimates of TB prevalence (all ages, all forms of TB) for 17 
    countries, before (in blue) and after (in red) results from national prevalence surveys 
    became available. Panels are ordered according to the size of the before-after difference. 
    The wide uncertainty interval of the post-survey estimate for the United Republic of 
    Tanzania is because laboratory challenges meant that it was only possible to directly 
    estimate the prevalence of smear-positive (as opposed to bacteriologically confirmed) TB.\n',
    x = unit(0, "lines"), 
    y = unit(0, "lines"),
    hjust = 0, vjust = 0,
    gp = gpar(fontsize = 16))

p2b <- arrangeGrob(p2, main = title.grob)

pdf(width=10, height=8, file='fig/prepostsurvey.pdf')
print(p2b)
dev.off()

write.csv(post2, file='output/before_after.csv', row.names=F)
 

# change in incidence pre vs post survey
post3 <- post2[year>=2012]
load('Rdata/est.Rdata')
load('../gtb2015/Rdata/old.Rdata')
post3b <- merge(post3, est[,.(iso3,year,inc,inc.lo,inc.hi)], by=c('iso3','year'), all.x=T, all.y=F)
lcty <- post3b$iso3[post3b$year==2013]
inc2013 <- old[iso3 %in% lcty & year==2013, .(iso3,year,oinc=inc,oinc.lo=inc.lo,oinc.hi=inc.hi)]

load('../gtb2014/Rdata/old.Rdata')
lcty <- post3b$iso3[post3b$year==2012]
inc2012 <- old[iso3 %in% lcty & year==2012, .(iso3,year,oinc=inc,oinc.lo=inc.lo,oinc.hi=inc.hi)]
addinc <- rbind(inc2012, inc2013)
post3c <- merge(post3b, addinc, by=c('iso3','year'))
post3c$change <- post3c$inc - post3c$oinc

sel <- post3c$iso3=='IDN'
sel2 <- old$iso3=='IDN' & old$year==2012
post3c$oinc[sel] <- old$inc[sel2]
post3c$oinc.lo[sel] <- old$inc.lo[sel2]
post3c$oinc.hi[sel] <- old$inc.hi[sel2]

p3 <- qplot(oinc, reorder(country, change), data=post3c, geom='point', size=I(3), colour=I('lightblue'), 
            fill=I('lightblue')) +
  geom_segment(aes(x=oinc.lo, xend=oinc.hi, y=country, yend=country), colour=I('lightblue'),
               size=I(3)) +
  geom_point(aes(inc, country), colour=I('red'), size=I(1)) +
  geom_segment(aes(x=inc.lo, xend=inc.hi, y=country, yend=country), 
               colour=I('red'), size=I(1)) +
  facet_wrap(~region, nrow=1, scales='free_y') +
  scale_x_log10(breaks=c(50,100,200,500)) +
  xlab('Incidence per 100,000 population per year (log scale)') + ylab('') +
  theme_bw(base_size=18)
pdf(width=10, height=6, file='fig/prepostsurveyInc.pdf')
print(p3)
dev.off()

