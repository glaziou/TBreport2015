##################################################
# Author:   Philippe Glaziou
#
# Global TB Report various calcs
##################################################
library(foreign)
library(data.table)

rm(list=ls())
source('fun.R')

m <- 100000











#--------------------------------------------
# Global IRR
#--------------------------------------------
load('Rdata/global.Rdata')
load('Rdata/est.Rdata')
phiv <- est[, list(hiv=sum(hiv*e.pop.num, na.rm=T)/sum(e.pop.num, na.rm=T), hiv.sd=sqrt(sum(hiv.sd^2, na.rm=T))), by=year]
(phiv)
global <- merge(global, phiv, by='year')
set.seed(123)
for (i in 1:dim(global)[1]){
  print(paste(round(i/dim(global)[1] * 100, 1), '%',.sdp=''))
  inc.h <- betaop(global$inc.h[i]/m, global$hiv[i], global$inc.h.sd[i]/m, global$hiv.sd[i]/m, op='/', nsim=m)
  inc.nh <-betaop(global$inc.nh[i]/m, (1 - global$hiv[i]), global$inc.h.sd[i]/m, global$hiv.sd[i]/m, op='/', nsim=m)
  
  ratio <- betaop(inc.h[1], inc.nh[1], inc.h[2], inc.nh[2], op='/', nsim=m)
  global$irr[i] <- ratio[1]          # TB incidence rate ratio HIV+/HIV-
  global$irr.sd[i] <- ratio[2]
  global$irr.lo[i] <- ratio[3]
  global$irr.hi[i] <- ratio[4]
  global$inc.hpos[i] <- inc.h[1]     # TB incidence rate among HIV+
  global$inc.hpos.sd[i] <- inc.h[2]
  global$inc.hpos.lo[i] <- inc.h[3]
  global$inc.hpos.hi[i] <- inc.h[4]
  global$inc.hneg[i] <- inc.nh[1]    # TB incidence rate among HIV-
  global$inc.hneg.sd[i] <- inc.nh[2]
  global$inc.hneg.lo[i] <- inc.nh[3]
  global$inc.hneg.hi[i] <- inc.nh[4]
}






#--------------------------------------------
# Regional IRR
#--------------------------------------------
load('Rdata/regional.Rdata')
load('Rdata/est.Rdata')
phiv <- est[, list(hiv=sum(hiv*e.pop.num, na.rm=T)/sum(e.pop.num, na.rm=T), hiv.sd=sqrt(sum(hiv.sd^2, na.rm=T))), by=list(g.whoregion, year)]
(phiv)
regional <- merge(regional, phiv, by=c('g.whoregion', 'year'))
set.seed(105)
for (i in 1:dim(regional)[1]){
  print(paste(round(i/dim(regional)[1] * 100, 1), '%',.sdp=''))
  inc.h <- betaop(regional$inc.h[i]/m, regional$hiv[i], regional$inc.h.sd[i]/m, regional$hiv.sd[i]/m, op='/', nsim=m)
  inc.nh <-betaop(regional$inc.nh[i]/m, (1 - regional$hiv[i]), regional$inc.h.sd[i]/m, regional$hiv.sd[i]/m, op='/', nsim=m)
  
  ratio <- betaop(inc.h[1], inc.nh[1], inc.h[2], inc.nh[2], op='/', nsim=m)
  regional$irr[i] <- ratio[1]          # TB incidence rate ratio HIV+/HIV-
  regional$irr.sd[i] <- ratio[2]
  regional$irr.lo[i] <- ratio[3]
  regional$irr.hi[i] <- ratio[4]
  regional$inc.hpos[i] <- inc.h[1]     # TB incidence rate among HIV+
  regional$inc.hpos.sd[i] <- inc.h[2]
  regional$inc.hpos.lo[i] <- inc.h[3]
  regional$inc.hpos.hi[i] <- inc.h[4]
  regional$inc.hneg[i] <- inc.nh[1]    # TB incidence rate among HIV-
  regional$inc.hneg.sd[i] <- inc.nh[2]
  regional$inc.hneg.lo[i] <- inc.nh[3]
  regional$inc.hneg.hi[i] <- inc.nh[4]
}


pdf(width=8, height=8, file='tbreport/figs/fig_chapter2_incidenceAmongHIV_byRegion.pdf')
qplot(year, inc.hpos*100, data=regional, geom='line', colour=I('red')) + 
  xlab('') + ylab('Incidence in HIV-positive (%/year)') +
  geom_ribbon(aes(year, ymin=inc.hpos.lo*100, ymax=inc.hpos.hi*100), fill=I('red'), alpha=I(.4)) +
  facet_wrap(~g.whoregion)
dev.off()


pdf(width=8, height=8, file='tbreport/figs/fig_chapter2_IRR_byRegion.pdf')
qplot(year, irr, data=regional, geom='line', colour=I('red')) + 
  xlab('') + ylab('Incidence in HIV-positive (%/year)') +
  geom_ribbon(aes(year, ymin=irr.lo, ymax=irr.hi), fill=I('red'), alpha=I(.4)) +
  facet_wrap(~g.whoregion)
dev.off()

# combining both at global level

p1 <- qplot(year, inc.hpos*100, data=global, geom='line', colour=I('red')) + 
  xlab('') + ylab('Incidence in HIV-positive (%/year)') +
  geom_ribbon(aes(year, ymin=inc.hpos.lo*100, ymax=inc.hpos.hi*100), fill=I('red'), alpha=I(.4)) 

p2 <- qplot(year, irr, data=global, geom='line', colour=I('red')) + 
  xlab('') + ylab('TB incidence Rate Ratio (HIV+/HIV-)') +
  geom_ribbon(aes(year, ymin=irr.lo, ymax=irr.hi), fill=I('red'), alpha=I(.4)) 

multiplot(p1, p2, cols=2)




#--------------------------------------------
# PAF - Amy Bloom question on resurgence
#--------------------------------------------
paf <- function(pe, rr) pe * (rr - 1) / (pe * (rr - 1) + 1)

load('Rdata/global.Rdata')
load('Rdata/unaids.Rdata')
hiv <- unaids[, list(hiv=sum(hiv=hiv1549*pop15/100, na.rm=TRUE)), by=year]
pe <- global[, list(pop=e.pop.num, tb=inc.num, tbinc=inc, hiv.tb=as.vector(unlist(divXY(inc.h, inc, inc.se^2, inc.h.se^2)[1]))), by=year]
res <- merge(pe, hiv, by='year')
res$hivprev <- res$hiv/res$pop
res$hivattr.tb <- paf(res$hivprev, 20)
res$nonattr.tb <- res$tb * (1 - res$hivattr.tb)
res$nonattr.tb.inc <- res$nonattr.tb/res$pop*100000
(res)

qplot(year, tbinc, data=res, geom='line') + 
  geom_line(aes(year, nonattr.tb.inc), colour=I('blue')) +
#  expand_limits(y=0) +
  xlab('') + ylab('Rate per 100,000 / year')






#--------------------------------------------
# Annex 1 numbers
#--------------------------------------------
# pop with VR data
table(est$source.mort[est$year==2012])
length(unique(est$iso3[est$source.mort %in% c('VR','VR imputed')])) # 125
table(est$year)

unique(est$iso3[est$source.mort %in% c('VR','VR imputed') & est$g.hbc22=='high'])

load('Rdata/vr.Rdata')
table(vr$keep.vr)
table(vr$year)
# data points per country
vr <- data.table(vr); setkey(vr, iso3, year)
(res <- vr[, list(points=sum(keep.vr)), by=iso3])
mean(res$points); sd(res$points)

length(unique(est$iso3[est$source.mort %in% c('Indirect')]))

(res <- est[year==2012, list(tbdeaths=sum(mort.nh.num, na.rm=T)), by=source.mort])
(res$tbdeaths/sum(res$tbdeaths))











### FROM previous year's report_misc.R

#--------------------------------------------
# TBHIV outcomes 2010
#--------------------------------------------
load('Rdata/tx.Rdata')
select <- tx$year==2010

# death rate, pooled cases
(tx[select][, list(cohort=sum(rowSums(cbind(hiv.new.sp.coh, hiv.new.snep.coh, hiv.ret.coh), na.rm=T)),
                   death.rate=sum(rowSums(cbind(hiv.new.sp.died, hiv.new.snep.died, hiv.ret.died), na.rm=T)) / 
                     sum(rowSums(cbind(hiv.new.sp.coh, hiv.new.snep.coh, hiv.ret.coh), na.rm=T)), 
                   ncountries=sum(!is.na(hiv.new.sp.died) & !is.na(hiv.new.snep.died) & !is.na(hiv.ret.died))), 
            by=g.whoregion])

(tx[select][, list(cohort=sum(rowSums(cbind(hiv.new.sp.coh, hiv.new.snep.coh, hiv.ret.coh), na.rm=T)),
                   death.rate=sum(rowSums(cbind(hiv.new.sp.died, hiv.new.snep.died, hiv.ret.died), na.rm=T)) / 
                     sum(rowSums(cbind(hiv.new.sp.coh, hiv.new.snep.coh, hiv.ret.coh), na.rm=T)), 
                   ncountries=sum(!is.na(hiv.new.sp.died) & !is.na(hiv.new.snep.died) & !is.na(hiv.ret.died))), 
            by=year])






#--------------------------------------------
# TB-AIDS deaths by sex
#--------------------------------------------
load('Rdata/est.Rdata')
load('Rdata/cty.Rdata')
mf <- read.csv('input/HIV/aids_deaths_sex.csv')
dim(mf)
tofix <- setdiff(unique(as.character(mf$country)), cty$country) # countries in unaids and not in cty
# !!!!!
# run this block line by line 
# checking the LHS and RHS returned values before running the complete line
# !!!!!
levels(mf$country)[levels(mf$country) %in% tofix][1] <- iso('Bolivia')[[2]]
levels(mf$country)[levels(mf$country) %in% tofix][1] <- iso('Ivoire')[[2]]
levels(mf$country)[levels(mf$country) %in% tofix][1] <- iso('Democratic People')[[2]]
levels(mf$country)[levels(mf$country) %in% tofix][1] <- iso('Lao')[[2]]
levels(mf$country)[levels(mf$country) %in% tofix][1] <- iso('Tome')[[2]]

levels(mf$country)[levels(mf$country) %in% tofix][1] <- 'Sudan'
levels(mf$country)[levels(mf$country) %in% tofix][1] <- 'South Sudan'

levels(mf$country)[levels(mf$country) %in% tofix][1] <- iso('Kingdom')[[2]]
levels(mf$country)[levels(mf$country) %in% tofix][1] <- iso('Venezuela')[[2]]
levels(mf$country)[levels(mf$country) %in% tofix] # should now be empty
setdiff(unique(as.character(mf$country)), cty$country) # and this as well
# end !!!!!

mf2 <- merge(cty[, list(iso3,country,g.whoregion)], mf, by='country', all.x=FALSE, all.y=TRUE)
mf <- mf2
mf <- data.table(mf); setkey(mf, iso3, year)
dim(mf)
mf2 <- merge(est[year==2011, list(iso3,e.pop.num,mort.h,mort.h.se,mort.h.lo,mort.h.hi,
                                 mort.h.num,mort.h.lo.num,mort.h.hi.num)], mf, by='iso3', all.x=TRUE)
dim(mf2)
mf <- mf2
sum(mf$mort.h.num[is.na(mf$sexratio)], na.rm=TRUE) # AIDS-TB deaths with no sex disaggregation, n=1490
sum(mf$mort.h.num, na.rm=TRUE) # out of 429193



# impute missing sexratio
r <- weighted.mean(mf$sexratio, weights=mf$aids_deathsM + mf$aids_deathsF, na.rm=TRUE) # r = 4.232
mf$imp.sexratio <- is.na(mf$sexratio)
mf$sexratio[is.na(mf$sexratio)] <- r


# Calculate number of female AIDS deaths, regional and global M/F
# r = m/f
# t = m + f = rf + f = f(r+1)
# f = t / (r+1)

mf$fem <- mf$mort.h.num / (1 + mf$sexratio)
mf$fem.lo <- mf$mort.h.lo.num / (1 + mf$sexratio)
mf$fem.hi <- mf$mort.h.hi.num / (1 + mf$sexratio)
mf$mort.h.f <- mf$mort.h / (1 + mf$sexratio)
mf$mort.h.f.se <- mf$mort.h.se / (1 + mf$sexratio)
mf$mort.h.m <- mf$mort.h / (1 + 1/mf$sexratio)
mf$mort.h.m.se <- mf$mort.h.se / (1 + 1/mf$sexratio)

sum(mf$fem, na.rm=TRUE) # 210,000
m <- 100000
(globalf <- add.rv(mf$mort.h.f/m, r.sd=mf$mort.h.f.se/m, weights=mf$e.pop.num)) # 190,000 - 230,000
(globalm <- add.rv(mf$mort.h.m/m, r.sd=mf$mort.h.m.se/m, weights=mf$e.pop.num)) 
(globalt <- add.rv(mf$mort.h/m, r.sd=mf$mort.h.se/m, weights=mf$e.pop.num))
(global.ratio <- divXY(globalm$r, globalf$r, globalm$r.sd^2, globalf$r.sd^2)) # 1.04 (SD 0.054)

wr <- as.character(unique(mf$g.whoregion[!is.na(mf$g.whoregion)]))

region <- data.frame(g.whoregion=wr)
region$male <- region$female <- region$sexratio <- NA
region$male.lo <- region$female.lo <- region$male.hi <- region$female.hi <- region$sexratio.se <- NA

for (i in wr){
  sel <- region$g.whoregion==i
  mf2 <- mf[mf$g.whoregion==i]

  W <- add.rv(mf2$mort.h.f/m, r.sd=mf2$mort.h.f.se/m, weights=mf2$e.pop.num) 
  M <- add.rv(mf2$mort.h.m/m, r.sd=mf2$mort.h.m.se/m, weights=mf2$e.pop.num) 
  tot <- add.rv(mf2$mort.h/m, r.sd=mf2$mort.h.se/m, weights=mf2$e.pop.num)
  ratio <- divXY(M$r, W$r, M$r.sd^2, W$r.sd^2)
  
  region$male[sel] <- M$r.num
  region$male.lo[sel] <- M$r.lo.num
  region$male.hi[sel] <- M$r.hi.num
  region$female[sel] <- W$r.num
  region$female.lo[sel] <- W$r.lo.num
  region$female.hi[sel] <- W$r.hi.num
  region$total[sel] <- tot$r.num
  region$total.lo[sel] <- tot$r.lo.num
  region$total.hi[sel] <- tot$r.hi.num
  region$sexratio[sel] <- ratio[[1]]
  region$sexratio.se[sel] <- sqrt(ratio[[2]])
}
rownames(region) <- region$g.whoregion
region$g.whoregion <- NULL
(region)
glob <- data.frame(male=globalm$r.num, male.lo=globalm$r.lo.num, male.hi=globalm$r.hi.num,
                   female=globalf$r.num, female.lo=globalf$r.lo.num, female.hi=globalf$r.hi.num,
                   total=globalt$r.num, total.lo=globalt$r.lo.num, total.hi=globalt$r.hi.num,
                   sexratio=global.ratio[[1]], sexratio.se=sqrt(global.ratio[[2]]))
rownames(glob) <- 'World'

out <- rbind(region,glob)
out$sexratio.lo <-out$sexratio + qnorm(0.025) * out$sexratio.se
out$sexratio.hi <-out$sexratio + qnorm(0.975) * out$sexratio.se
out2 <- out[, c('male','male.lo','male.hi','female','female.lo','female.hi',
                'total','total.lo','total.hi',
                'sexratio','sexratio.lo','sexratio.hi')]
out2 <- signif(out2,3)
(out2 <- out2[c(3,2,1,4,6,5,7),])
write.csv(out2, file='output/sexratio_tbhiv_deaths.csv', row.names=TRUE)


# Plot sex ratio versus HIV
# load('Rdata/tbhiv.Rdata')
# mf2 <- merge(mf, tbhiv[year==2011, list(iso3,hiv,tbhiv)])
# qplot(hiv, sexratio, data=subset(mf2, sexratio<5 & !is.na(sexratio)), colour=g.whoregion)


# Plot sex ratio versus WHO region
out$region <- rownames(out)
out$region <- factor(out$region, rev(c('AMR','AFR','EMR','EUR','SEA','WPR','World')))
levels(out$region) <- rev(c('Americas','Africa','Eastern Mediterranean','Europe','South East Asia','Western Pacific','World'))

qplot(region, sexratio, ymin=sexratio.lo, ymax=sexratio.hi, data=out, geom='pointrange') +
  xlab('') + ylab('Sex Ratio (M/F)') + coord_flip()
ggsave(file='tbreport/figs/sexratio_tbhiv_deaths.pdf')








#--------------------------------------------
# Mortality vs MDR
#--------------------------------------------
load('Rdata/mdr.Rdata')
load('Rdata/est.Rdata')
load('Rdata/cty.Rdata')

est2 <- ddply(as.data.frame(est[source.mort=='VR']), .(iso3), function(x) mort <- last(x$mort.nh))
dim(est2)
names(est2)[2] <- 'mort.nh'

mdr2 <- ddply(as.data.frame(mdr), .(iso3), function(x) mdr <- last(x$e.new.mdr.pcnt))
dim(mdr2)
names(mdr2)[2] <- 'mdr'

dta <- merge(mdr2, est2, by='iso3')
dim(dta)
dta <- merge(dta, cty[, list(iso3, g.income)], by='iso3', all.x=TRUE, all.y=FALSE)
dta$income <- 'Low/Middle Income'
dta$income[dta$g.income=='HIC' & !is.na(dta$g.income)] <- 'High Income'

qplot(mdr, mort.nh, data=subset(dta, mdr>0 & mort.nh>0)) +
  facet_wrap(~income, scales='free') +
  xlab('MDR in previously untreated cases (%, log scale)') + 
  ylab('Mortality rate per 100,000/year (log scale)') +
  scale_x_log10() + 
  scale_y_log10() +
  geom_smooth(method='lm')
ggsave('tbreport/figs/mdr_mortality.pdf')

summary(lm(mdr ~ mort.nh, data=subset(dta, mdr>0 & mort.nh>0 & g.income=='HIC')))
summary(lm(mdr ~ mort.nh, data=subset(dta, mdr>0 & mort.nh>0 & g.income!='HIC')))






#--------------------------------------------
# Misc
#--------------------------------------------
# numbers for Annex 1
load('Rdata/vr.Rdata')
dim(vr)
table(vr$keep.vr)
(chin <- 21+21-6)
(1919-30-chin) # VR data points, 30 outliers, 36 imputed in China and India

vr$m <- vr$keep.vr
vr$m[vr$iso3=='CHN' & vr$year<2006] <- FALSE
vr$m[vr$iso3=='CHN' & vr$year!=2005] <- FALSE
table(vr$m) # 1898 data points retained for analysis

mean(table(vr$iso3[vr$m])) # mean number of points per country
sd(table(vr$iso3[vr$m]))




#--------------------------------------------
# China MDR prevalence in general population
#--------------------------------------------
# chn <- est['CHN'][year %in% c(2000, 2010)]
# mdr00 <- cii(263, 20)
# mdr10 <- cii(241, 13)
# ab1 <- get.beta(mdr00$prob, mdr00$se)
# ab2 <- get.beta(mdr10$prob, mdr10$se)
# 
# ab3 <- get.beta(chn$prev[1] / m, (chn$prev.hi[1] - chn$prev.lo[1]) / m / 4)
# ab4 <- get.beta(chn$prev[2] / m, (chn$prev.hi[2] - chn$prev.lo[2]) / m / 4)
# 
# nsim <- 50000
# mx1 <- rbeta(nsim, ab1[1], ab1[2]) * rbeta(nsim, ab3[1], ab3[2])
# mx2 <- rbeta(nsim, ab2[1], ab2[2]) * rbeta(nsim, ab4[1], ab4[2])
# 
# mxp1 <- mx1 * chn$e.pop.num[1]
# mxp2 <- mx2 * chn$e.pop.num[2]
# 
# (res1 <- c(mean(mxp1), quantile(mxp1, probs=c(0.025, 0.975))))
# (res2 <- c(mean(mxp2), quantile(mxp2, probs=c(0.025, 0.975))))
# 
# (sum((mxp1 - mxp2) > 0) / nsim)





#--------------------------------------------
# Orphans due to TB, using an obfuscated model
#--------------------------------------------
(est[, sum(mort.h.num, na.rm=T)/sum(aids.mort*e.pop.num/m,na.rm=T), by=year])
# ok, about 22% aids deaths due to TB - not a surprise really

(hiv[, sum(aids.orphans.num, na.rm=T) / sum(total.orphans, na.rm=T), by=year])
# 32% of orphans due to aids
# let us boldly assume that o_h / o_t = m_h / (m_t * a)
# where o_h is aids orphans, o_t is tb orphans (HIV-), m_h is 
# aids mortality (30) and m_t is TB mortality excluding HIV (21),
# taken from Global Statistics Report 2010,
# a is the proportion of TB cases aged 15-55 (the others aren't supposed to have kids)
# Then o_t = (o_h * m_t * a / m_h) + k * o_h 
# where k is the 22% proportion of AIDS deaths due to TB.

nsim <- 50000

# get a
(am <- tb2[year>2004][, sum(new.sp.m1524 + new.sp.m2534 + new.sp.m3544 + new.sp.m4554 +
  new.sp.f1524 + new.sp.f2534 + new.sp.f3544 + new.sp.f4554, na.rm=TRUE)/ sum(new.sp, na.rm=T), by=year])
a <- mean(am$V1); a.sd <- sd(am$V1)
ab.a <- get.beta(a, a.sd)
sim.a <- rbeta(nsim, ab.a[1], ab.a[2]) # distribution of a

# get k
out <- est[, add.rv(aids.mort/m, aids.mort.lo/m, aids.mort.hi/m, weights=e.pop.num), by=year]
global2 <- cbind(global, out)
km <- with(subset(global2, year==2009), ratio.rates(mort.h/m, mort.h.lo/m, mort.h.hi/m,
                    r, r.lo, r.hi))
ab.km <- get.beta(km$r, (km$r.hi - km$r.lo) / 4)
sim.k <- rbeta(nsim, ab.km[1], ab.km[2])

# distribution of orphans from AIDS deaths
hiv2 <- merge(hiv, pop)
hiv2$oh <- hiv2$aids.orphans.num / hiv2$e.pop.num
hiv2$oh.lo <- hiv2$aids.orphans.lo.num / hiv2$e.pop.num
hiv2$oh.hi <- hiv2$aids.orphans.hi.num / hiv2$e.pop.num

out <- hiv2[year==2009][, add.rv(oh, oh.lo, oh.hi, weights=e.pop.num)]
ab.oh <- get.beta(out$r, (out$r.hi - out$r.lo) / 4)
sim.oh <- rbeta(nsim, ab.oh[1], ab.oh[2])
(mean(sim.oh)*pop)

# TB mortality HIV-
out <- global[year==2009][, list(mort.nh, mort.nh.lo, mort.nh.hi)]
ab.mt <- get.beta(out$mort.nh/m, (out$mort.nh.hi - out$mort.nh.lo) /m / 4)
sim.mt <- rbeta(nsim, ab.mt[1], ab.mt[2])

# AIDS mortality
out <- est[year==2009][, add.rv(aids.mort/m, aids.mort.lo/m, aids.mort.hi/m, weights=e.pop.num)]
ab.mh <- get.beta(out$r, (out$r.hi - out$r.lo) / 4)
sim.mh <- rbeta(nsim, ab.mh[1], ab.mh[2])


# get results
(hiv[year==2009][, sum(aids.orphans.num, na.rm=T)])

pop <- global$e.pop.num[global$year==2009]

sim.ot.nh <- (sim.oh * sim.mt * sim.a / sim.mh)
qot.nh <- quantile(sim.ot.nh, probs=c(0.025, 0.975), na.rm=TRUE) * pop
(ot.nh <- mean(sim.ot.nh, na.rm=TRUE) * pop) # TB orphans excl HIV
(ot.nh.lo <- qot.nh[1])
(ot.nh.hi <- qot.nh[2])
  
  
sim.ot.h <- sim.k * sim.oh
qot.h <- quantile(sim.ot.h, probs=c(0.025, 0.975)) * pop
(ot.h <- mean(sim.ot.h) * pop) # TB orphans among AIDS orphans
(ot.h.lo <- qot.h[1])
(ot.h.hi <- qot.h[2])
  
sim.ot <- sim.ot.nh + sim.ot.h
qot <- quantile(sim.ot, probs=c(0.025, 0.975)) * pop
(ot <- mean(sim.ot) * pop) # total TB orphans
(ot.lo <- qot[1])
(ot.hi <- qot[2])





#--------------------------------------------
# Lives saved
#--------------------------------------------
load('Rdata/tb.Rdata')
load('Rdata/dots.Rdata')

dot <- dots[, list(dots=sum(tot.newrel.dots, na.rm=T)), by=year]
tb[, list(dots=sum(tot.newrel, na.rm=T)), by=year][year>2007]
dot <- rbind(dot, tb[, list(dots=sum(tot.newrel, na.rm=T)), by=year][year>2007])

# global lives saved since 1995
treated = sum(dot$dots[dot$year>1994], na.rm=T)
saved = treated /3
(signif(treated, 3) / 1000000)
(signif(saved, 3) / 1000000)

# In India
dot <- dots['IND', list(dots=sum(tot.newrel.dots, na.rm=T)), by=year]
tb['IND', list(dots=sum(tot.newrel, na.rm=T)), by=year][year>2007]
dot <- rbind(dot, tb['IND', list(dots=sum(tot.newrel, na.rm=T)), by=year][year>2007])

treated.IND = sum(dot$dots[dot$year>1994], na.rm=T)
saved.IND = treated.IND / 3
(signif(treated.IND, 3) / 1000000)
(signif(saved.IND, 3) / 1000000)

ind.tsr <- c(.79, .79, .82, .84, .82, .84, .85, .87, .86, .86, .86, .86, .87, .87, .88, .88, .88)
(cured.IND <- sum(ind.tsr * dot$dots[dot$year>1994])/1000000)









#--------------------------------------------
# Exports 
#--------------------------------------------

# -> IER
write.csv(est[, list(iso3, year, e.pop.num, round(mort.nh.num), round(mort.nh.lo.num), 
                     round(mort.nh.hi.num))], file='output/ier.csv', row.names=FALSE)









#--------------------------------------------
# Successfully treated
#--------------------------------------------
tb2 <- merge(tb, tx[, list(iso3, year, c.new.tsr, c.ret.tsr, c.ret.nrel.tsr)], by=c('iso3','year'))

wmean <- function(x, w) sum(x*w, na.rm=T) / sum(w, na.rm=T)

out <- tb2[year %in% 1995:2011, list(notif.all=sum(c.notified, na.rm=T), 
                                     notif.new=sum(c.new, na.rm=T),
                                     notif.ret=sum(c.ret, na.rm=T),
#                                     nm.c.tsr=sum(!is.na(c.tsr)),
#                                     nm.c.new.tsr=sum(!is.na(c.new.tsr)),
#                                     nm.c.ret.tsr=sum(!is.na(c.ret.tsr)),
                                     cure.rate.new=wmean(x=c.new.tsr, w=c.new)/100, 
                                     cure.rate.ret=wmean(x=c.ret.tsr, w=c.ret)/100
                                     ), by=year]
out <- within(out, {
  cured <- notif.new * cure.rate.new + notif.ret * cure.rate.ret
  cure.rate.all <- cured / notif.all
})
(out)

newrow <- tb[year %in% 2012:2013, list(notif.all=sum(c.notified, na.rm=T), 
                                       notif.new=sum(c.new, na.rm=T), 
                                       notif.ret=sum(c.ret, na.rm=T),
                                       cure.rate.new=last(out$cure.rate.new), 
                                       cure.rate.ret=last(out$cure.rate.ret), 
                                       cure.rate.all=last(out$cure.rate.all)
                                       ), by=year]
newrow$notif.new[2] <- tb[year==2013, sum(c.newunk, na.rm=T)]

newrow$cured <- newrow$notif.all * newrow$cure.rate.all
res <- data.frame(rbind(out, newrow))
res[, 5:7] <- signif(res[ , 5:7], 2)
(res)
(sum(res$cured))













