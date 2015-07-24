##################################################
# Author:   Philippe Glaziou
# Updated:  16/07/2013
#
# Global TB Report EPI tables
##################################################
library(data.table)

rm(list=ls())
source('fun.R')

m <- 100000




#--------------------------------------------
# Load datasets 
#--------------------------------------------
load('Rdata/tb.Rdata')
load('Rdata/est.Rdata')
load('Rdata/global.Rdata')
load('Rdata/regional.Rdata')
load('Rdata/hbc.Rdata')





#--------------------------------------------
# Table 2.1 epi burden (absolute numbers)
#--------------------------------------------
yr <- 2013
tab1.1 <- subset(est, g.hbc22=='high' & !is.na(g.hbc22) & year==yr, 
  select=c('country','e.pop.num','mort.nh.num','mort.nh.lo.num','mort.nh.hi.num',
           'mort.h.num','mort.h.lo.num','mort.h.hi.num',
    'prev.num','prev.lo.num','prev.hi.num','inc.num','inc.lo.num','inc.hi.num',
    'inc.h.num','inc.h.lo.num','inc.h.hi.num'))

tab1.2 <- subset(hbc, year==yr, 
  select=c('e.pop.num','mort.nh.num','mort.nh.lo.num','mort.nh.hi.num',
           'mort.h.num','mort.h.lo.num','mort.h.hi.num',
    'prev.num','prev.lo.num','prev.hi.num','inc.num','inc.lo.num','inc.hi.num',
    'inc.h.num','inc.h.lo.num','inc.h.hi.num'))

tab1.3 <- subset(regional, year==yr, 
  select=c('g.whoregion','e.pop.num','mort.nh.num','mort.nh.lo.num','mort.nh.hi.num',
           'mort.h.num','mort.h.lo.num','mort.h.hi.num',
    'prev.num','prev.lo.num','prev.hi.num','inc.num','inc.lo.num','inc.hi.num',
    'inc.h.num','inc.h.lo.num','inc.h.hi.num'))

tab1.4 <- subset(global, year==yr, 
  select=c('e.pop.num','mort.nh.num','mort.nh.lo.num','mort.nh.hi.num',
           'mort.h.num','mort.h.lo.num','mort.h.hi.num',
    'prev.num','prev.lo.num','prev.hi.num','inc.num','inc.lo.num','inc.hi.num',
    'inc.h.num','inc.h.lo.num','inc.h.hi.num'))


setnames(tab1.1, 'country','rowname')
setnames(tab1.3, 'g.whoregion', 'rowname')
tab1.2 <- cbind(rowname='High-burden countries', tab1.2)
tab1.4 <- cbind(rowname='Global', tab1.4)

tab1 <- rbind(tab1.1, tab1.2, tab1.3, tab1.4)

tab1s <- as.data.frame(tab1)
#tab1s[, 3:14] <- signif(tab1s[, 3:14], 3)

tab1s$mort.nh.bounds <- paste("(", as.character(signif(tab1s$mort.nh.lo.num,3)),"-", 
                              as.character(signif(tab1s$mort.nh.hi.num, 3)), ")", sep='')

tab1s$mort.h.bounds <- paste("(", as.character(signif(tab1s$mort.h.lo.num,3)),"-", 
                              as.character(signif(tab1s$mort.h.hi.num, 3)), ")", sep='')

tab1s$prev.bounds <- paste("(", as.character(signif(tab1s$prev.lo.num,3)),"-", 
                              as.character(signif(tab1s$prev.hi.num, 3)), ")", sep='')

tab1s$inc.bounds <- paste("(", as.character(signif(tab1s$inc.lo.num,3)),"-", 
                           as.character(signif(tab1s$inc.hi.num, 3)), ")", sep='')

tab1s$inc.h.bounds <- paste("(", as.character(signif(tab1s$inc.h.lo.num,3)),"-", 
                          as.character(signif(tab1s$inc.h.hi.num, 3)), ")", sep='')

tab1s <- tab1s[, c('rowname', 'mort.nh.num','mort.nh.bounds',
                   'mort.h.num','mort.h.bounds',
                   'prev.num','prev.bounds',
                   'inc.num','inc.bounds',
                   'inc.h.num','inc.h.bounds')]

(tab1s)

write.csv(tab1s, file='tbreport/tables/tab2_1.csv', row.names=FALSE)






#--------------------------------------------
# Table 2.2 epi burden (rates)
#--------------------------------------------
tab1b.1 <- subset(est, g.hbc22=='high' & !is.na(g.hbc22) & year==yr, 
  select=c('country','e.pop.num','mort.nh','mort.nh.lo','mort.nh.hi',
    'prev','prev.lo','prev.hi','inc','inc.lo','inc.hi',
    'tbhiv','tbhiv.lo','tbhiv.hi'))

tab1b.2 <- subset(hbc, year==yr, 
  select=c('e.pop.num','mort.nh','mort.nh.lo','mort.nh.hi',
    'prev','prev.lo','prev.hi','inc','inc.lo','inc.hi',
    'tbhiv','tbhiv.lo','tbhiv.hi'))

tab1b.3 <- subset(regional, year==yr, 
  select=c('g.whoregion','e.pop.num','mort.nh','mort.nh.lo','mort.nh.hi',
    'prev','prev.lo','prev.hi','inc','inc.lo','inc.hi',
    'tbhiv','tbhiv.lo','tbhiv.hi'))

tab1b.4 <- subset(global, year==yr, 
  select=c('e.pop.num','mort.nh','mort.nh.lo','mort.nh.hi',
    'prev','prev.lo','prev.hi','inc','inc.lo','inc.hi',
    'tbhiv','tbhiv.lo','tbhiv.hi'))


setnames(tab1b.1, 'country','rowname')
setnames(tab1b.3, 'g.whoregion', 'rowname')
tab1b.2 <- cbind(rowname='High-burden countries', tab1b.2)
tab1b.4 <- cbind(rowname='Global', tab1b.4)

tab1b <- rbind(tab1b.1, tab1b.2, tab1b.3, tab1b.4)

tab1bs <- as.data.frame(tab1b)

tab1bs$mort.nh.bounds <- paste("(", as.character(signif(tab1bs$mort.nh.lo,3)),"-", 
                              as.character(signif(tab1bs$mort.nh.hi, 3)), ")", sep='')

tab1bs$prev.bounds <- paste("(", as.character(signif(tab1bs$prev.lo,3)),"-", 
                           as.character(signif(tab1bs$prev.hi, 3)), ")", sep='')

tab1bs$inc.bounds <- paste("(", as.character(signif(tab1bs$inc.lo,3)),"-", 
                          as.character(signif(tab1bs$inc.hi, 3)), ")", sep='')

tab1bs$tbhiv.bounds <- paste("(", as.character(signif(tab1bs$tbhiv.lo,3)),"-", 
                            as.character(signif(tab1bs$tbhiv.hi, 3)), ")", sep='')

tab1bs <- tab1bs[, c('rowname', 'mort.nh','mort.nh.bounds',
                   'prev','prev.bounds',
                   'inc','inc.bounds',
                   'tbhiv','tbhiv.bounds')]

(tab1bs)

write.csv(tab1bs, file='tbreport/tables/tab2_2.csv', row.names=FALSE)






#--------------------------------------------
# Top 10 countries
#--------------------------------------------
# top 10 in decreasing order of estimated incident cases
(top1 <- head(est[year==yr][order(inc.num, decreasing=T), list(country, inc.num)], 10))

# top 10 in decreasing order of risk of TB (incidence rate)
(top2 <- head(est[year==yr][order(inc, decreasing=T), list(country, inc)], 10))

# top 10 in decreasing order of estimated TB deaths (including HIV+)
(top3 <- head(est[year==yr][order(mort.num, decreasing=T), list(country, mort.num)], 10))

# top 10 in decreasing order of risk of dying from TB (mortality rate including HIV+)
(top4 <- head(est[year==yr][order(mort, decreasing=T), list(country, mort)], 10))

# top 10 in decreasing order of HIV+ TB cases
(top5 <- head(est[year==yr][order(inc.h.num, decreasing=T), list(country, mort)], 10))

# top 10 in decreasing order of MDR-TB cases
#(top6 <- head(est[year==2012][order(inc.h.num, decreasing=T), list(country, mort)], 10))


write.csv(top1, file='tbreport/tables/tab2_3_top10inc.csv', row.names=FALSE)
write.csv(top2, file='tbreport/tables/tab2_3_top10incnum.csv', row.names=FALSE)
write.csv(top3, file='tbreport/tables/tab2_3_top10deaths.csv', row.names=FALSE)
write.csv(top4, file='tbreport/tables/tab2_3_top10mort.csv', row.names=FALSE)
write.csv(top5, file='tbreport/tables/tab2_3_top10hivtb.csv', row.names=FALSE)
#write.csv(top6, file='tbreport/tables/tab2_3_top10mdr.csv', row.names=FALSE)







#--------------------------------------------
# Table 4.1 notifications
#--------------------------------------------
tab2a <- function(yr) {
	tb2 <- within(as.data.frame(tb), {
		new.pulm <- rowSums(cbind(new.labconf, ret.rel.labconf, new.clindx, ret.rel.clindx), na.rm=T)
    c.ep <- rowSums(cbind(new.ep, ret.rel.ep), na.rm=T)
	})
  tb2$new.prop.pulm <- tb2$new.pulm / (tb2$new.pulm + tb2$c.ep)

	hbc <- tb2[tb$g.hbc22 == 'high' & !is.na (tb$g.hbc22) & tb$year==yr, ]
	select <- c('country','c.newinc','new.labconf','new.clindx','ret.rel.labconf', 'ret.rel.clindx','c.ep',
			'ret.nrel', 'new.pulm', 'new.prop.pulm')
	hbc <- hbc[, select]
	hbc <- hbc[order(hbc$country), ]
	
	nr <- data.frame(country='High burden countries', c.newinc=sum(hbc$c.newinc, na.rm=T), 
                   new.labconf=sum(hbc$new.labconf, na.rm=T), 
                   new.clindx=sum(hbc$new.clindx, na.rm=T), 
                   ret.rel.labconf=sum(hbc$ret.rel.labconf, na.rm=T),
                   ret.rel.clindx=sum(hbc$ret.rel.clindx, na.rm=T),
                   c.ep=sum(hbc$c.ep, na.rm=T),
                   ret.nrel=sum(hbc$ret.nrel, na.rm=T), 
	                 new.pulm=sum(hbc$new.pulm, na.rm=T),
	                 new.prop.pulm=sum(hbc$new.pulm, na.rm=T)/(sum(hbc$new.pulm, na.rm=T) + sum(hbc$c.ep, na.rm=T)))
	hbc <- rbind(hbc, nr)
	tab2a <- hbc
	return(tab2a)
}

tab2b <- function(yr) {
  tb2 <- within(as.data.frame(tb), {
    new.pulm <- rowSums(cbind(new.labconf, ret.rel.labconf, new.clindx, ret.rel.clindx), na.rm=T)
    c.ep <- rowSums(cbind(new.ep, ret.rel.ep), na.rm=T)
  })
  tb2$new.prop.pulm <- tb2$new.pulm / (tb2$new.pulm + tb2$c.ep)
  
  hbc <- tb2[!is.na (tb$g.hbc22) & tb$year==yr, ]
  select <- c('country','g.whoregion', 'c.newinc','new.labconf','new.clindx','ret.rel.labconf', 'ret.rel.clindx','c.ep',
              'ret.nrel', 'new.pulm', 'new.prop.pulm')
  hbc <- hbc[, select]
  
  wr <- with(hbc, aggregate(hbc[, -c(1,2,11)], by=list(region=g.whoregion), sum, na.rm=T))
  wr$new.prop.pulm <- wr$new.pulm / (wr$new.pulm + wr$c.ep)
  wr <- wr[order(wr$region), ]
  
  nr <- data.frame(region='Global', c.newinc=sum(hbc$c.newinc, na.rm=T), 
                   new.labconf=sum(hbc$new.labconf, na.rm=T), 
                   new.clindx=sum(hbc$new.clindx, na.rm=T), 
                   ret.rel.labconf=sum(hbc$ret.rel.labconf, na.rm=T),
                   ret.rel.clindx=sum(hbc$ret.rel.clindx, na.rm=T),
                   c.ep=sum(hbc$c.ep, na.rm=T),
                   ret.nrel=sum(hbc$ret.nrel, na.rm=T), 
                   new.pulm=sum(hbc$new.pulm, na.rm=T),
                   new.prop.pulm=sum(hbc$new.pulm, na.rm=T)/(sum(hbc$new.pulm, na.rm=T) + sum(hbc$c.ep, na.rm=T)))
  hbc <- rbind(wr, nr)
  tab2b <- hbc
  return(tab2b)
}


tab2.1 <- tab2a(yr)
tab2.2 <- tab2b(yr)

names(tab2.1)[1] <- names(tab2.2)[1] <- 'rowname'
tab2 <- rbind(tab2.1, tab2.2)
(tab2)

write.csv(tab2, file='tbreport/tables/tab4_1.csv', row.names=FALSE)





#--------------------------------------------
# Table 4.5 CDR all forms
#--------------------------------------------
hbc <- within(as.data.frame(hbc), {
  cdr <- c.newinc / inc.num
  cdr.lo <- c.newinc / inc.hi.num
  cdr.hi <- c.newinc / inc.lo.num
})

regional$cdr <- regional$c.newinc / regional$inc.num
regional$cdr.lo <- regional$c.newinc / regional$inc.hi.num
regional$cdr.hi <- regional$c.newinc / regional$inc.lo.num

global$cdr <- global$c.newinc / global$inc.num
global$cdr.lo <- global$c.newinc / global$inc.hi.num
global$cdr.hi <- global$c.newinc / global$inc.lo.num

est$cdr <- est$newinc/est$inc
est$cdr.lo <- est$newinc/est$inc.hi
est$cdr.hi <- est$newinc/est$inc.lo


tab4.1 <- subset(as.data.frame(est), g.hbc22=='high' & !is.na(g.hbc22) & year %in% c(1995, 2000, 2005, 2010, yr), 
  select=c('country','year','cdr','cdr.lo','cdr.hi'))

tab4.2 <- subset(as.data.frame(hbc), year %in% c(1995, 2000, 2005, 2010, yr), 
  select=c('year','cdr','cdr.lo','cdr.hi'))

tab4.3 <- subset(as.data.frame(regional), year %in% c(1995, 2000, 2005, 2010, yr), 
  select=c('g.whoregion','year','cdr','cdr.lo','cdr.hi'))

tab4.4 <- subset(as.data.frame(global), year %in% c(1995, 2000, 2005, 2010, yr), 
  select=c('year','cdr','cdr.lo','cdr.hi'))


names(tab4.1)[1] <- 'rowname'
names(tab4.3)[1] <- 'rowname'
tab4.2 <- cbind(rowname='High-burden countries', tab4.2)
tab4.4 <- cbind(rowname='Global', tab4.4)

tab4.1r <- reshape(tab4.1, idvar='rowname', timevar='year', direction='wide')
tab4.2r <- reshape(tab4.2, idvar='rowname', timevar='year', direction='wide')
tab4.3r <- reshape(tab4.3, idvar='rowname', timevar='year', direction='wide')
tab4.4r <- reshape(tab4.4, idvar='rowname', timevar='year', direction='wide')


tab4 <- rbind(tab4.1r, tab4.2r, tab4.3r, tab4.4r)

tab4s <- tab4
#tab4s[, 2:16] <- signif(tab4s[, 2:16], 2) * 100

tab4s$cdr.1995.bounds <- paste("(", as.character(signif(tab4s$cdr.lo.1995,3)),"-", 
                               as.character(signif(tab4s$cdr.hi.1995, 3)), ")", sep='')

tab4s$cdr.2000.bounds <- paste("(", as.character(signif(tab4s$cdr.lo.2000,3)),"-", 
                                as.character(signif(tab4s$cdr.hi.2000, 3)), ")", sep='')

tab4s$cdr.2005.bounds <- paste("(", as.character(signif(tab4s$cdr.lo.2005,3)),"-", 
                                as.character(signif(tab4s$cdr.hi.2005, 3)), ")", sep='')

tab4s$cdr.2010.bounds <- paste("(", as.character(signif(tab4s$cdr.lo.2010,3)),"-", 
                                as.character(signif(tab4s$cdr.hi.2010, 3)), ")", sep='')

tab4s$cdr.2013.bounds <- paste("(", as.character(signif(tab4s$cdr.lo.2013,3)),"-", 
                                as.character(signif(tab4s$cdr.hi.2013, 3)), ")", sep='')


tab4s <- tab4s[, c('rowname', 'cdr.1995','cdr.1995.bounds',
                   'cdr.2000','cdr.2000.bounds',
                   'cdr.2005','cdr.2005.bounds',
                   'cdr.2010','cdr.2010.bounds',
                   'cdr.2013','cdr.2013.bounds')]

(tab4s)

write.csv(tab4s, file='tbreport/tables/tab4_5.csv', row.names=FALSE)










#--------------------------------------------
# Table Annex 1 - Sources of data
#--------------------------------------------
tab <- as.data.frame(est[year==yr][, list(g.whoregion, country, source.inc, source.mort, source.prev)][order(g.whoregion, country)])
(tab)
write.csv(tab, 'tbreport/tables/annex1_data_sources.csv', row.names=FALSE)





#--------------------------------------------
# Table Annex 1 - VR countries
#--------------------------------------------
tabvr <- est[, list(vr.points=sum(source.mort=='VR'), tot.points=sum(!is.na(year))), by=iso3]

load('Rdata/vr.Rdata')
tabvr2 <- vr[, list(dropped.points=sum(!keep.vr), 
          highest.cov=signif(max(estcov)*100,2), 
          lowest.garbage=signif(min(garbage)*100,2)), by=iso3]

tabvr3 <- merge(tabvr, tabvr2)

write.csv(tabvr3[, list(iso3, tot.points, vr.points, dropped.points, highest.cov, lowest.garbage)], file='output/sourcemort.csv', row.names=F)


# add imputation methods
tabvr4 <- est[year==yr, list(iso3, source.mort)]
tabvr5 <- merge(tabvr3, tabvr4, all=TRUE)
table(tabvr5$source.mort)

fa <- c('BEN','BFA','BDI','CMR','CAF','TCD','COM','COD','COG','CIV','GNQ','GAB','GIN',
        'MDG','MLI','NER','RWA','SEN','TGO')
eco <- c(fa, 'PNG','NAM','ETH','ERI','LAO','SSD','SDN')

#tabvr5$source.mort[tabvr5$iso3 %in% eco] <- 'Ecological'
tabvr5$source.mort[tabvr5$source.mort=='Indirect'] <- 'CFR'
tabvr5$source.mort[tabvr5$source.mort=='VR imputed'] <- 'VR'
tabvr5$source.mort[is.na(tabvr5$source.mort)] <- 'CFR'
table(tabvr5$source.mort)
tabvr5$vr.points[is.na(tabvr5$vr.points)] <- 0

write.csv(tabvr5[, list(iso3, vr.points, source.mort)], file='output/sourcemort.csv', row.names=F)












