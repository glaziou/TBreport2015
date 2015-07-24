##################################################
# Author:   Philippe Glaziou
#
# GLOBAL TB REPORT figs
##################################################
library(ggplot2)
library(grid)
library(gridExtra)




#--------------------------------------------
# Load datasets 
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




#----------------------------------------------------
# Global rates of incidence, prevalence and mortality
#----------------------------------------------------
# global plot
pdf(width=8, height=4, file='tbreport/figs/fig2_6_global.pdf')
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



pdf(width=8, height=4, file='tbreport/figs/fig4_1_globalNotificationRates.pdf')
qplot(year, inc, data=global, geom='line', colour=I('darkgreen')) +
  geom_ribbon(aes(year, ymin=inc.lo, ymax=inc.hi), fill=I('darkgreen'), alpha=0.4) +
  geom_line(aes(year, newinc)) +   ylab('Rate per 100,000 population/year') + xlab('') +
  expand_limits(y=0) +
  theme_bw(base_size=10) 
dev.off()




#----------------------------------------------------
# Global numbers (incidence and mortality)
#----------------------------------------------------
pdf(width=8, height=4, file='tbreport/figs/fig2_2_global_num.pdf')
mil <- 1e6
p1 <- qplot(year, inc.num/mil, data=global, geom='line', colour=I('black')) +
  geom_ribbon(aes(year, ymin=inc.lo.num/mil, ymax=inc.hi.num/mil), fill=I('grey'), alpha=0.8) +
  geom_line(aes(year, inc.h.num/mil), colour=I('red')) +
  geom_line(aes(year, c.newinc/mil)) + 
  geom_ribbon(aes(year, ymin=inc.h.lo.num/mil, ymax=inc.h.hi.num/mil), 
              fill=I('red'), alpha=0.4) +
                ylab('Cases/year (millions)') + xlab('') +
                expand_limits(y=0) +
                theme_bw(base_size=11) 

p2 <- qplot(year, mort.nh.num/mil, data=global, geom='line', colour=I('blue')) +
  geom_ribbon(aes(year, ymin=mort.nh.lo.num/mil, ymax=mort.nh.hi.num/mil), fill=I('blue'), alpha=0.3) +
  geom_line(aes(year,mort.h.num/mil), colour=I('red')) +
  geom_ribbon(aes(year, ymin=mort.h.lo.num/mil, ymax=mort.h.hi.num/mil), fill=I('red'), alpha=0.3) +
  geom_line(aes(year, mort.num/mil), colour=I('black')) +
  geom_ribbon(aes(year, ymin=mort.lo.num/mil, ymax=mort.hi.num/mil), fill=I('grey'), alpha=0.8) +
  ylab('TB deaths/year (millions)') + xlab('') +
  expand_limits(y=0) +
  theme_bw(base_size=11) +
  theme(legend.position='none')


multiplot(p1, p2, cols=2)
dev.off()









#----------------------------------------------------
# Regional notification rate plus inc
#----------------------------------------------------
# chapter 2
pdf(width=8, height=8, file='tbreport/figs/fig2_7_incidenceByRegion.pdf')
qplot(year, inc, data=as.data.frame(regional), geom='line', colour=I('blue')) +
  geom_ribbon(aes(year, ymin=inc.lo, ymax=inc.hi), fill=I('blue'), alpha=0.4) +
  geom_line(aes(year, inc.h), colour=I('red')) +
  geom_line(aes(year, newinc)) + 
  geom_ribbon(aes(year, ymin=inc.h.lo, ymax=inc.h.hi), 
              fill=I('red'), alpha=0.4) +
  facet_wrap(~g.whoregion, scales='free_y') +
  xlab('') + ylab('Rate per 100,000 population/year') +
  expand_limits(y=0) +
  theme_bw(base_size=10) 
dev.off()



# chapter 3
pdf(width=8, height=8, file='tbreport/figs/fig4_2_notifByRegion.pdf')
qplot(year, inc, data=as.data.frame(regional), geom='line', colour=I('darkgreen')) +
    geom_ribbon(aes(year, ymin=inc.lo, ymax=inc.hi), fill=I('darkgreen'), alpha=0.4) +
    geom_line(aes(year, newinc)) + 
    facet_wrap(~g.whoregion, scales='free_y') +
    xlab('') + ylab('Rate per 100,000 population/year') +
    expand_limits(y=0) +
    theme_bw(base_size=10) 
dev.off()









#------------------------------------------------------
# Prevalence and Mortality by WHO region
#------------------------------------------------------
regional.ff2 <- ddply(as.data.frame(regional.ff), .(g.whoregion), transform,
  target.prev=prev[1]/2, target.mort=mort.nh[1]/2, target.morth=mort.h[15]/2)

  
pdf(width=8, height=8, file='tbreport/figs/fig2_10_prevalenceByRegion.pdf')
qplot(year, prev, data=regional.ff2, geom='line', colour=I('blue'), linetype=forecast) +
    geom_ribbon(aes(year, ymin=prev.lo, ymax=prev.hi), fill=I('blue'), alpha=0.4) +
    geom_hline(aes(yintercept=target.prev), linetype=2) +
    facet_wrap(~g.whoregion, scales='free_y') +
    xlab('') + ylab('Rate per 100,000 population/year') +
    expand_limits(y=0) +
    theme_bw(base_size=10) +
    theme(legend.position='none')
dev.off()

pdf(width=8, height=8, file='tbreport/figs/fig2_13_mortalityByRegion.pdf')
qplot(year, mort.nh, data=regional.ff2, geom='line', colour=I('blue'), linetype=forecast) +
    geom_ribbon(aes(year, ymin=mort.nh.lo, ymax=mort.nh.hi), fill=I('blue'), alpha=0.4) +
    geom_hline(aes(yintercept=target.mort), linetype=2) +
    facet_wrap(~g.whoregion, scales='free_y') +
    xlab('') + ylab('Rate per 100,000 population/year') +
    expand_limits(y=0) +
    theme_bw(base_size=10) +
    theme(legend.position='none')
dev.off()


pdf(width=8, height=8, file='tbreport/figs/fig2_13b_mortalityHIVposByRegion.pdf')
qplot(year, mort.h, data=subset(regional.ff2, year<2013), geom='line', colour=I('blue')) +
  geom_ribbon(aes(year, ymin=mort.h.lo, ymax=mort.h.hi), fill=I('blue'), alpha=0.4) +
  geom_hline(aes(yintercept=target.morth), linetype=2) +
  facet_wrap(~g.whoregion, scales='free_y') +
  xlab('') + ylab('Rate per 100,000 population/year') +
  expand_limits(y=0) +
  theme_bw(base_size=10) +
  theme(legend.position='none')
dev.off()






#----------------------------------------------------
# AFR trends in AIDS-TB mortality
#----------------------------------------------------
pdf(width=8, height=8, file='tbreport/figs/fig_afr_mortality.pdf')
qplot(year, mort.h, data=regional[g.whoregion=='AFR'], geom='line', colour=I('red')) +
  geom_ribbon(aes(year, ymin=mort.h.lo, ymax=mort.h.hi), fill=I('red'), alpha=0.4) +
  geom_line(aes(year, mort), colour=I('grey')) +
  geom_ribbon(aes(year, ymin=mort.lo, ymax=mort.hi), fill=I('grey'), alpha=0.7) +  
  geom_line(aes(year, mort.nh), colour=I('blue')) +
  geom_ribbon(aes(year, ymin=mort.nh.lo, ymax=mort.nh.hi), fill=I('blue'), alpha=0.2) +  
  xlab('') + ylab('Rate per 100,000 population/year') +
  expand_limits(y=0) +
  theme(legend.position='none')
dev.off()






#----------------------------------------------------
# Incidence in 22 HBCs
#----------------------------------------------------
hest <- subset(est, g.hbc22=='high' & !is.na(g.hbc22))
levels(hest$country)[match('Democratic Republic of the Congo', levels(hest$country))] <- 'DR Congo'
levels(hest$country)[match('United Republic of Tanzania', levels(hest$country))] <- 'UR Tanzania'

hest <- within(hest, {
  inc.h.lo <- vlohi(inc.h/m, inc.h.se/m)[1, ]*m
  inc.h.hi <- vlohi(inc.h/m, inc.h.se/m)[2, ]*m
})

#pdf(width=8, height=8, file='tbreport/figs/fig2_8_incidence_22hbc.pdf')
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
#dev.off()

# notifs and incidence (chapter 3)
#pdf(width=8, height=8, file='tbreport/figs/fig3_3_incidenceRates_22hbc.pdf')
p2 <- qplot(year, inc, data=as.data.frame(hest), geom='line', colour=I('blue')) +
  geom_ribbon(aes(year, ymin=inc.lo, ymax=inc.hi), fill=I('blue'), alpha=0.4) +
  geom_line(aes(year, newinc)) + 
  facet_wrap(~country, scales='free_y') +
  xlab('') + ylab('Rate per 100,000 population/year') +
  expand_limits(y=0) +
  theme_bw(base_size=8) 
print(facetAdjust(p2))
#dev.off()






#----------------------------------------------------
# Incidence in top 10 countries
#----------------------------------------------------
(top1 <- head(est[year==yr][order(inc.num, decreasing=T), list(country, inc.num)], 10))
levels(top1$country)[match('Democratic Republic of the Congo', levels(top1$country))] <- 'DR Congo'

pdf(width=8, height=8, file='tbreport/figs/fig2_3_inc_top10inc.pdf')
p1 <- qplot(year, inc, data=est[country %in% as.character(top1$country)], geom='line', colour=I('blue')) +
  geom_ribbon(aes(year, ymin=inc.lo, ymax=inc.hi), fill=I('blue'), alpha=0.4) +
  geom_line(aes(year, newinc)) + 
  facet_wrap(~country, scales='free_y') +
  xlab('') + ylab('Rate per 100,000 population/year') +
  expand_limits(y=0) +
  theme_bw(base_size=8) 
print(p1)
#print(facetAdjust(p2))
dev.off()








#------------------------------------------------------
# Prevalence and Mortality in 22 HBCs
#------------------------------------------------------
hbc.ff2 <- ddply(as.data.frame(hbc.ff), .(iso3), transform,
  target.prev=prev[1]/2, target.mort=mort.nh[1]/2)

dim(hbc.ff2)
load('Rdata/cty.Rdata')
hbc.ff2 <- merge(hbc.ff2, cty[, list(iso3, country)], by='iso3', all.x=TRUE, all.y=FALSE)
dim(hbc.ff2)

levels(hbc.ff2$country)[match('Democratic Republic of the Congo', levels(hbc.ff2$country))] <- 'DR Congo'
levels(hbc.ff2$country)[match('United Republic of Tanzania', levels(hbc.ff2$country))] <- 'UR Tanzania'


#pdf(width=8, height=8, file='tbreport/figs/fig2_14_prevalence_22hbc.pdf')
p1 <- qplot(year, prev, data=hbc.ff2, geom='line', colour=I('blue'), linetype=forecast) +
    geom_ribbon(aes(year, ymin=prev.lo, ymax=prev.hi), fill=I('blue'), alpha=0.4) +
    geom_hline(aes(yintercept=target.prev), linetype=2) +
    facet_wrap(~country, scales='free_y') +
    xlab('') + ylab('Rate per 100,000 population/year') +
    expand_limits(y=0) +
    theme_bw(base_size=8) +
    theme(legend.position='none')
print.facetAdjust(facetAdjust(p1))
#ggsave(file='tbreport/figs/fig2_14_prevalence_22hbc.pdf', width=8, height=8)
#dev.off()




#pdf(width=8, height=8, file='tbreport/figs/fig2_14_mortalityRates_hbc22.pdf')
hbc.ff3 <- merge(hbc.ff2, est[, list(iso3,year,vr.tbrate.raw)], by=c('iso3','year'), all.x=TRUE, all.y=FALSE)
p2 <- qplot(year, mort.nh, data=hbc.ff3, geom='line', colour=I('blue'), linetype=forecast) +
    geom_ribbon(aes(year, ymin=mort.nh.lo, ymax=mort.nh.hi), fill=I('blue'), alpha=0.4) +
    geom_hline(aes(yintercept=target.mort), linetype=2) +
    geom_point(aes(year, vr.tbrate.raw), shape=I(4)) +
    facet_wrap(~country, scales='free_y') +
    xlab('') + ylab('Rate per 100,000 population/year') +
    expand_limits(y=0) +
    theme_bw(base_size=8) +
    theme(legend.position='none')
print(facetAdjust(p2))
#dev.off()
#ggsave(file='tbreport/figs/fig2_14b_mortality_22hbc.pdf', width=8, height=8)


#pdf(width=8, height=8, file='tbreport/figs/fig2_14c_mortalityHIVpos_hbc22.pdf')
p2 <- qplot(year, mort.h, data=subset(est, g.hbc22=='high'), geom='line', colour=I('blue')) +
  geom_ribbon(aes(year, ymin=mort.h.lo, ymax=mort.h.hi), fill=I('blue'), alpha=0.4) +
  facet_wrap(~country, scales='free_y') +
  xlab('') + ylab('Rate per 100,000 population/year') +
  expand_limits(y=0) +
  theme_bw(base_size=8)
print(facetAdjust(p2))
#ggsave(file='tbreport/figs/fig2_14c_mortalityHIVpos_22hbc.pdf', width=8, height=8)
#dev.off()








#----------------------------------------------------
# Maps
#----------------------------------------------------
library(whomap)
yr <- 2013
dta <- est[year==2013, list(iso3, g.hbc22, newinc, snewinc, inc, inc.num, tbhiv, mort, mort.nh, prev, 
                            mort.num, source.inc, source.mort, source.prev)]
# VR countries (2012)
dta$var <- dta$source.mort %in% c("VR","VR imputed")
whomap(X=dta, Z=scale_fill_manual("VR/Survey", values=c('grey','darkgreen')), legendpos='none')

ggsave(file='tbreport/figs/fig2_11_map_VRcountries.pdf')
write.csv(dta, file='output/vrCountries.csv', row.names=FALSE)

# Source TBHIV (2012)
# dta$var <- dta$source.tbhiv
# whomap(X=dta, Z=scale_fill_brewer("TBHIV", palette="Set1"))
# 
# ggsave(file='tbreport/figs/map_sourceTBHIV.pdf')


# source Incidence (2012)
dta$var <- dta$source.inc %in% c("Capture-recapture","High income","Survey")
whomap(X=dta, Z=scale_fill_manual("Surveillance", values=c('grey','darkgreen')), legendpos='none')

ggsave(file='tbreport/figs/fig2_1_map_incmethod.pdf')
write.csv(dta[, list(iso3, consult=var)], file='output/inc_method.csv', row.names=F)



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

ggsave(file='tbreport/figs/fig2_5_map_ctyconsult.pdf')


# incidence
dta$var <- cut(dta$inc, c(0, 25, 50, 100, 200, 300, Inf), 
               c('0-24.9', '25-49.9', '50-99', '100-199.9', '200-299', '300+'), 
               right=F, ordered_result=T)
whomap(X=dta, Z=scale_fill_brewer("Incidence\nper 100,000", palette="YlGnBu", type="seq"))
ggsave(file='tbreport/figs/fig2_5_map_incidenceRates.pdf')



# mortality HIV-
dta$var <- cut(dta$mort.nh, c(0, 1, 2, 5, 10, 20, 40, Inf), 
               c('0-0.9', '1-1.9', '2-4.9', '5-9.9', '10-19','20-39','40+'), 
               right=F, ordered_result=T)
whomap(X=dta, Z=scale_fill_brewer("Mortality\nper 100,000", palette="YlGnBu", type="seq"))
ggsave(file='tbreport/figs/fig2_12_map_mortalityRatesHIVneg.pdf')



# prevalence
# dta$var <- cut(dta$prev, c(0, 5, 10, 20, 50, 100, 200, Inf), 
#                c('0-4.9', '5-9.9', '10-19', '20-49', '50-99', '100-199','200+'), 
#                right=F, ordered_result=T)
# whomap(X=dta, Z=scale_fill_brewer("Prevalence\nper 100,000", palette="YlGnBu", type="seq"))
# ggsave(file='tbreport/figs/fig2_4_map_prevalence.pdf')


# HIV prevalence in TB
dta$var <- cut(dta$tbhiv*100, c(0, 1, 2, 5, 10, 20, 50, Inf), 
               c('0-0.9', '1-1.9', '2-4.9', '5-9.9', '10-19', '20-49','50+'), 
               right=F, ordered_result=T)
whomap(X=dta, Z=scale_fill_brewer("HIV Prevalence\namong TB (%)", palette="YlGnBu", type="seq"))
ggsave(file='tbreport/figs/fig2_4_map_HIVamongTB.pdf')



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
ggsave(file='tbreport/figs/fig2_9_map_prevalence_surveys.pdf')



# MN ratio (2011)
# dta$var <- dta$mort / dta$sm.newinc
# dta$var <- cut(dta$var, c(0, 0.05, 0.1, 0.2, 0.3, 0.5, Inf), 
#                c('0-4.9', '5-9.9', '10-19', '20-29', '30-49', '50+'), 
#                right=F, ordered_result=T)
# whomap(X=dta, Z=scale_fill_brewer("M:N Ratio (%)", palette="YlGnBu", type="seq"))
# 
# ggsave(file='tbreport/figs/map_MNratio.pdf')



























