#------------------------------------------------------
# Various functions
# PG, 30/06/2015
#------------------------------------------------------
# load deps
library('ggplot2')
library('data.table')

# load datasets
# load('Rdata/est.Rdata')
load('Rdata/cty.Rdata')
load('Rdata/tb.Rdata')
#load('Rdata/sty.Rdata')
#load('Rdata/unaids.Rdata')
load('Rdata/pop.Rdata')



# constants
m <- 1e5




#------------------------------------------------------
# get country name from iso3 code
isocty <- function(iso)cty[iso, country]


# get iso code from country name
iso <- function (ct = "China") {
  suppressWarnings({
    country <- as.character(cty$country[grep(ct, cty$country, ignore.case=TRUE)])
    iso3 <- as.character(cty$iso3[cty$country %in% country])
    out <- cbind(iso3, country)
    return(out)
  })
}



# operations on two random variates distributed beta
betaop <- function(ev1, ev2, se1, se2, op="*", distr=F, nsim=500000){
  # defaults to product for a ratio, op='/'
  # ev1, ev2 = expected values
  # se1, se2 = standard errors
  # dist = returns distribution in a vector of size nsim

  if (is.na(ev1) || is.na(ev2) || is.na(se1) || is.na(se2)) return (c(ev=NA, se=NA, lo=NA, hi=NA))
  
  stopifnot(ev1<=1 & ev2<=1)
  stopifnot(se1<=0.5 & se2<=0.5)
  stopifnot(ev1>=0 & ev2>=0)
  stopifnot(se1>=0 & se2>=0)
  stopifnot(op %in% c('*','/','+','-'))
  
  get.beta <- function(ev, sd){
    S = (ev * (1 - ev) / sd^2) - 1
    a = S * ev
    b = S * (1 - ev)
    return(c(a = a, b = b))    
  }
  
  par1 <- get.beta(ev1, se1)
  par2 <- get.beta(ev2, se2)
  
  out <- get(op)(rbeta(nsim, par1[1], par1[2]), rbeta(nsim, par2[1], par2[2]))
  if (distr) return(out)
  else return(c(ev=mean(out), 
           se=sd(out), 
           lo=quantile(out, prob=0.025, names=FALSE, na.rm=TRUE), 
           hi=quantile(out, prob=0.975, names=FALSE, na.rm=TRUE)))  
}


# error propagation - requires dplyr
pgate <- function(.data, f) {
  # errors in x and y are stored in x.se and y.se
  exprs = list(
    deparse(f[[3]]),
    sapply(all.vars(f[[3]]), function(v) {
      dfdp = deparse(D(f[[3]], v))
      sprintf('(%s.se*(%s))^2', v, dfdp)
    }) %>%
      paste(collapse='+') %>%
      sprintf('sqrt(%s)', .)
  )
  names(exprs) = c(
    deparse(f[[2]]),
    sprintf('%s.se', deparse(f[[2]]))
  )
  
  .data %>% mutate_(.dots=exprs)
}



# product of two random variables X and Y using Taylor expansion approx
prodXY <- function(X, Y, varX, varY, covXY=0){
  eXY <- X * Y + covXY
  varXY <- X^2*varY + Y^2*varX + 2*X*Y*covXY + varX*varY + covXY^2
  return(list("E(XY)"=eXY, "Var(XY)"=varXY))
}


# ratio of two random variables using Taylor expansion
divXY <- function(X, Y, varX, varY, covXY=0){
  eXY <- X/Y - covXY/Y^2 + X*varY/Y^3
  varXY <- varX/Y^2 - 2*X*covXY/Y^3 + X^2*varY/Y^4
  return(list("E(X/Y)"=eXY, "Var(X/Y)"=varXY))
}


# EV and variance of the log of a RV
# based on Taylor expansion
logX <- function(X, varX){ 
    # X = EV
    # varX = variance
    E.logX <- log(X) - varX / (2 * X^2)
    var.logX <- varX / X^2
    return(list("E[log(X)]" = E.logX, "Var[log(X)]" = var.logX))
}
    
# logit and invlogit
logit <- function (x) log (x/(1-x))
invlogit <- function (x) 1/(1+exp(-x))



# returns Beta shape and scale params using the method of moments
get.beta <- function(ev, sd){
  S = (ev * (1 - ev) / sd^2) - 1
  a = S * ev
  b = S * (1 - ev)
  return(c(a = a, b = b))
}


# returns gamma shape and scale params using the method of moments
get.gamma <- function(ev, sd){
  v = ev / sd^2
  r = ev^2 / sd^2
  return(c(r = r, v = v))
}


# generate low and high bounds assuming Beta distribution
lohi <- function(ev, se){
  par <- get.beta(ev, se)
  lo <- qbeta(0.025, par[1], par[2])
  hi <- qbeta(0.975, par[1], par[2])
  return(c(lo=lo, hi=hi))
}


# sum random variables
add.rv <- function (r, r.lo, r.hi, r.sd, weights = 1, method = "beta") 
{
  if (is.null(r) || length(r) == 0) 
    stop("Error: r must contain at least one value")
  if (sum(r < 0 & !is.na(r) & method == "beta")) 
    stop("Error: r must be positive with method 'beta'")
  if (sum(r > 1 & !is.na(r) & method == "beta")) 
    stop("Error: r must be between 0 and 1 with method 'beta'")
  if (missing(r.sd)) 
    r.sd <- (r.hi - r.lo)/4
  if (missing(r.lo) & !missing(r.sd)) 
    r.lo <- numeric()
  if (missing(r.hi) & !missing(r.sd)) 
    r.hi <- numeric()
  if (sum(r.lo < 0 & !is.na(r.lo) & method == "beta")) 
    stop("Error: r.lo must be positive with method 'beta'")
  if (sum(r.lo > 1 & !is.na(r.lo) & method == "beta")) 
    stop("Error: r.lo must be between 0 and 1 with method 'beta'")
  if (sum(r.hi < 0 & !is.na(r.hi) & method == "beta")) 
    stop("Error: r.hi must be positive with method 'beta'")
  if (sum(r.hi > 1 & !is.na(r.hi) & method == "beta")) 
    stop("Error: r.hi must be between 0 and 1 with method 'beta'")
  if (sum(r.sd > 1 & !is.na(r.sd) & method == "beta")) 
    stop("Error: sd must be between 0 and 1 with method 'beta'")
  if (sum(r[!is.na(r) & is.na(r.sd)])) 
    stop("Error: some values for r are supplied without uncertainty")
  if (sum(r.sd < 0 & !is.null(r.sd) & !is.na(r.sd))) 
    stop("Error: sd must be positive")
  if (!is.null(r.sd)) 
    v <- r.sd^2
  else v <- ((r.hi - r.lo)/4)^2
  sw <- ifelse(length(weights) > 1, sum(weights[!is.na(r)], 
                                        na.rm = TRUE), 1)
  out.m <- sum(r * weights, na.rm = TRUE)/sw
  out.v <- ifelse(length(weights) > 1, sum(v[!is.na(r)] * weights[!is.na(r)]^2, 
                                           na.rm = TRUE)/sw^2, sum(v))
  if (method == "beta") {
    S <- (out.m * (1 - out.m)/out.v) - 1
    a <- S * out.m
    b <- S * (1 - out.m)
    lo <- qbeta(0.025, a, b)
    hi <- qbeta(0.975, a, b)
  }
  else {
    lo <- qnorm(0.025, out.m, sqrt(out.v))
    hi <- qnorm(0.975, out.m, sqrt(out.v))
  }
  if (all(weights == 1)) 
    return(data.frame(r = out.m, r.lo = lo, r.hi = hi, r.sd = sqrt(out.v)))
  else return(data.frame(r = out.m, r.lo = lo, r.hi = hi, r.sd = sqrt(out.v), 
                         r.num = out.m * sw, r.lo.num = lo * sw, r.hi.num = hi * 
                           sw, e.pop.num = sw))
}


# ensemble
ensbeta <- function(xi, xi.sd){
    stopifnot(xi<1 & xi.sd<1)
    stopifnot(xi>0 & xi.sd>0)
    vget.beta <- Vectorize(get.beta, c('ev','sd'))
    w <- vget.beta(xi, xi.sd) - 1
    a <- sum(w[1, ])
    b <- sum(w[2, ])
    pw <- list(c = a+1, d = b+1)
    k <- pw$c + pw$d
    post.ev <- pw$c / k
    post.lo <- qbeta(0.025, pw$c, pw$d)
    post.hi <- qbeta(0.975, pw$c, pw$d)
    post.sd <- sqrt(pw$c * pw$d /(k^2 * (k + 1)))
    return(list(post.param = c(shape=pw$c, scale=pw$d),
                post.ev = post.ev,
                post.lo = post.lo,
                post.hi = post.hi,
                post.sd = post.sd))
}



# trends in beta-distributed variates, with autocorrelation
betatrend <- function(x, x.se, distr=F, nsim=1000){
  # x = best estimate
  # x.se = standard deviation
  stopifnot(x<=1)
  stopifnot(x.se<0.5)
  stopifnot(x>=0 & x.se>=0)
  stopifnot(length(x) == length(x.se))

  lg <- length(x)
  V <- matrix(NA, nrow=nsim, ncol=lg)

  for (i in 1:lg){
    par <- get.beta(x[i], x.se[i]) 
    V[ ,i] <- sort(rbeta(nsim, par[1], par[2])) 
    # NAs produced if par[1] < 0 | par[2] < 0
  }
  
  idx <- 1:lg
  trend <- function(x)coef(lm(log(x)~idx))[2]
  out <- apply(V, 1, trend)
 
  if (distr) return(out)
  else return(c(ev=mean(out), 
                se=sd(out), 
                lo=quantile(out, prob=0.025, names=FALSE, na.rm=TRUE), 
                hi=quantile(out, prob=0.975, names=FALSE, na.rm=TRUE)))    
}



# combine plots - source: R cookbook
multiplot <- function(..., plotlist=NULL, cols) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                       # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
}




# function to plot mortality trends - librarys est and cty data.tables
mplot <- function(iso, hiv=FALSE, ylog=TRUE){
  if (missing(iso) || !is.character(iso) || nchar(iso)!=3)
    stop("Error: iso must be a valid ISO3 country code")
  if (is.na(match('est',objects(envir=parent.frame(), all.names = TRUE))) || 
    is.na(match('cty',objects(envir=parent.frame(), all.names = TRUE))))
      stop("Error: est and cty must be loaded")
  iso <- toupper(iso)
  isocty <- function(iso)cty[iso, country]
  
  p <- qplot(year, mort.nh, data=subset(est, iso3==iso), geom='line', 
             main=paste('TB mortality in', as.character(isocty(iso)))) +
        geom_point(aes(year, vr.raw), shape=I(4), size=I(3)) +
        geom_ribbon(aes(year, ymin=mort.nh.lo, ymax=mort.nh.hi), fill=I('blue'), alpha=I(0.3)) +
        xlab('') + 
        geom_hline(y=est$mort.nh[est$iso3==iso][1]/2, linetype=2) +
        ylab('Rate per 100,000/year') +
        theme_bw(base_size=20)
  
  if (ylog) p <- p + coord_trans(y="log10") + ylab('Rate per 100,000/year (log scale)')
  
  q <- p + geom_line(aes(year, mort.h), colour=I('red')) + 
        geom_ribbon(aes(year, ymin=mort.h.lo, ymax=mort.h.hi), fill=I('red'), alpha=I(0.3))
  
  if (hiv) return(q)
  else return(p)
}




# function to plot incidence trends - librarys est and cty data.tables
iplot <- function(iso, hiv=FALSE, ylog=TRUE){
  if (missing(iso) || !is.character(iso) || nchar(iso)!=3)
    stop("Error: iso must be a valid ISO3 country code")
  if (is.na(match('est',objects(envir=parent.frame(), all.names = TRUE))) || 
    is.na(match('cty',objects(envir=parent.frame(), all.names = TRUE))))
    stop("Error: est and cty must be loaded")
  iso <- toupper(iso)
  isocty <- function(iso)cty[iso, country]
  
  p <- qplot(year, inc, data=subset(est, iso3==iso), geom='line', 
             main=paste('TB incidence in', as.character(isocty(iso)))) +
               geom_ribbon(aes(year, ymin=inc.lo, ymax=inc.hi), fill=I('blue'), alpha=I(0.3)) +
               geom_line(aes(year, newinc)) +
               xlab('') + ylab('Rate per 100,000/year') +
                theme_bw(base_size=20)
  
  q <- p + geom_line(aes(year, inc.h), colour=I('red')) + 
    geom_ribbon(aes(year, ymin=inc.h.lo, ymax=inc.h.hi), fill=I('red'), alpha=I(0.3))
  
  if (ylog) p <- p + coord_trans(y = "log10") + ylab('Rate per 100,000/year (log scale)')
  
  if (hiv) return(q)
  else return(p)
}



# function to plot prevalence trends - librarys est and cty data.tables
pplot <- function(iso, ylog=TRUE){
  if (missing(iso) || !is.character(iso) || nchar(iso)!=3)
    stop("Error: iso must be a valid ISO3 country code")
  if (is.na(match('est',objects(envir=parent.frame(), all.names = TRUE))) || 
    is.na(match('cty',objects(envir=parent.frame(), all.names = TRUE))))
    stop("Error: est and cty must be loaded")
  iso <- toupper(iso)
  isocty <- function(iso)cty[iso, country]
  
  p <- qplot(year, prev, data=subset(est, iso3==iso), geom='line', 
             main=paste('TB prevalence in', as.character(isocty(iso)))) +
               geom_ribbon(aes(year, ymin=prev.lo, ymax=prev.hi), fill=I('blue'), alpha=I(0.3)) +
               geom_hline(y=est$prev[est$iso3==iso][1]/2, linetype=2) +
               xlab('') + ylab('Rate per 100,000') +
                theme_bw(base_size=20)
  
  if (ylog) p <- p + coord_trans(y = "log10") + ylab('Rate per 100,000 pop (log scale)')
  
  return(p)
}



# function to plot tbhiv trends - librarys est and cty data.tables
hplot <- function(iso, toplot=T){
  if (missing(iso) || !is.character(iso) || nchar(iso)!=3)
    stop("Error: iso must be a valid ISO3 country code")
  if (is.na(match('est',objects(envir=parent.frame(), all.names = TRUE))) ||
    is.na(match('tb',objects(envir=parent.frame(), all.names = TRUE))) ||
    is.na(match('sty',objects(envir=parent.frame(), all.names = TRUE))) ||
    is.na(match('cty',objects(envir=parent.frame(), all.names = TRUE))))
    stop("Error: est, tb, sty and cty must be loaded")
  iso <- toupper(iso)
  isocty <- function(iso)cty[iso, country]
  
  sources <- length(table(est$source.tbhiv[est$iso3==iso]))
  
  dta <- merge(est[iso3==iso, list(iso3,year,tbhiv,tbhiv.lo,tbhiv.hi)], 
               tb[iso3==iso, list(iso3,year,tot.newrel,hivtest.f,hivtest.p,hivtest.pos.f,
                                  hivtest.pos.p,hiv.art.p,hiv.art.f)],
               by=c('iso3', 'year'), all.x=TRUE)
  dta <- merge(dta, sty[, list(iso3,year,tbhiv.surv.prev,tbhiv.surv.cil,tbhiv.surv.ciu,
                               tbhiv.sentin.prev,tbhiv.sentin.cil,tbhiv.sentin.ciu)],
               by=c('iso3', 'year'), all.x=TRUE)
  
  dta <- within(dta, {
    test.coverage <- 100 * hivtest.f / tot.newrel
    test.pos <- 100* hivtest.pos.f / hivtest.f
    art <- 100 * hiv.art.f / hivtest.pos.f
  })
  
  lastyr <- max(dta$year[dta$iso3==iso])
  dta$test.coverage[dta$year==lastyr] <- dta$hivtest.p[dta$year==lastyr] * 100 / dta$tot.newrel[dta$year==lastyr]
  dta$test.pos[dta$year==lastyr] <- dta$hivtest.pos.p[dta$year==lastyr] *100 / dta$tot.newrel[dta$year==lastyr]
  dta$art[dta$year==lastyr] <- dta$hiv.art.p[dta$year==lastyr] * 100 / dta$hivtest.pos.p[dta$year==lastyr]
  
  
  p <- qplot(year, tbhiv*100, data=dta, geom='line', 
             main=paste('HIV prevalence in TB in', as.character(isocty(iso)))) +
    scale_linetype_discrete(name = "Data Source") +
    geom_ribbon(aes(year, ymin=tbhiv.lo*100, ymax=tbhiv.hi*100), fill=I('red'), alpha=I(0.3)) +
    expand_limits(y=0) +
    xlab('') + ylab('Percent') +
    geom_point(aes(year, tbhiv.surv.prev), colour=I('blue'), shape=I(16), size=I(3)) +
    geom_linerange(aes(year, ymin=tbhiv.surv.cil, ymax=tbhiv.surv.ciu), colour=I('blue')) +         
    geom_point(aes(year, tbhiv.sentin.prev), colour=I('darkgreen'), shape=I(2), size=I(4)) +
    geom_linerange(aes(year, ymin=tbhiv.sentin.cil, ymax=tbhiv.sentin.ciu), colour=I('darkgreen'), size=I(3)) +
    geom_point(aes(year, test.pos), shape=I(4), size=I(3)) +
    theme_bw(base_size=20)
  if (toplot) return(p)
  else (return(dta))
}


# the above xplot functions can be combined with multiplot, e.g.:
# multiplot(iplot('CHN', hiv=T), mplot('CHN'), pplot('CHN'), hplot('CHN'), cols=2)


hivplot <- function(iso, start=1990, toplot=T){
  if (missing(iso) || !is.character(iso) || nchar(iso)!=3)
    stop("Error: iso must be a valid ISO3 country code")
  if (is.na(match('est',objects(envir=parent.frame(), all.names = TRUE))) ||
        is.na(match('cty',objects(envir=parent.frame(), all.names = TRUE))))
    stop("Error: est and cty must be loaded")
  iso <- toupper(iso)
  isocty <- function(iso)cty[iso, country]
  

  p <- qplot(year, hiv*100, data=est[iso3==iso & year>=start], geom='line', 
             main=paste('HIV prevalence (all ages) in', as.character(isocty(iso)))) +
    geom_ribbon(aes(year, ymin=hiv.lo*100, ymax=hiv.hi*100), fill=I('red'), alpha=I(0.3)) +
    expand_limits(y=0) +
    xlab('') + ylab('Percent') +
    theme_bw(base_size=20)
  if (toplot) return(p)
  else (return(dta))
}





# function to plot age-specific case trends (all forms) - librarys pop, tb and cty data.tables
ageplot <- function(iso, start=1990, all=TRUE, toplot=TRUE){
  if (missing(iso) || !is.character(iso) || nchar(iso)!=3)
    stop("Error: iso must be a valid ISO3 country code")
  if (is.na(match('tb',objects(envir=parent.frame(), all.names = TRUE))) || 
        is.na(match('cty',objects(envir=parent.frame(), all.names = TRUE))))
    stop("Error: tb and cty must be loaded")
  
  library('reshape')
  library('directlabels')
  
  iso <- toupper(iso)
  isocty <- function(iso)cty[iso, country]
  
  m <- 100000
  dta <- subset(tb, iso3==iso & year>=start)
  dta <- within(dta, {
    new04 <- rowSums(cbind(new.sp.m04, new.sp.f04, new.sn.m04*all, new.sn.f04*all, new.ep.m04*all, new.ep.f04*all), na.rm = TRUE)
    new0514 <- rowSums(cbind(new.sp.m514, new.sp.f514, new.sn.m514*all, new.sn.f514*all, new.ep.m514*all, new.ep.f514*all), na.rm = TRUE)
    new1524 <- rowSums(cbind(new.sp.m1524, new.sp.f1524, new.sn.m1524*all, new.sn.f1524*all, new.ep.m1524*all, new.ep.f1524*all), na.rm = TRUE)
    new2534 <- rowSums(cbind(new.sp.m2534, new.sp.f2534, new.sn.m2534*all, new.sn.f2534*all, new.ep.m2534*all, new.ep.f2534*all), na.rm = TRUE)
    new3544 <- rowSums(cbind(new.sp.m3544, new.sp.f3544, new.sn.m3544*all, new.sn.f3544*all, new.ep.m3544*all, new.ep.f3544*all), na.rm = TRUE)
    new4554 <- rowSums(cbind(new.sp.m4554, new.sp.f4554, new.sn.m4554*all, new.sn.f4554*all, new.ep.m4554*all, new.ep.f4554*all), na.rm = TRUE)
    new5564 <- rowSums(cbind(new.sp.m5564, new.sp.f5564, new.sn.m5564*all, new.sn.f5564*all, new.ep.m5564*all, new.ep.f5564*all), na.rm = TRUE)
    new65 <- rowSums(cbind(new.sp.m65, new.sp.f65, new.sn.m65*all, new.sn.f65*all, new.ep.m65*all, new.ep.f65*all), na.rm = TRUE)
    pop04 <- rowSums(cbind(e.pop.m04, e.pop.f04), na.rm=T)
    pop0514 <- rowSums(cbind(e.pop.m514, e.pop.f514), na.rm=T)
    pop1524 <- rowSums(cbind(e.pop.m1524, e.pop.f1524), na.rm=T)
    pop2534 <- rowSums(cbind(e.pop.m2534, e.pop.f2534), na.rm=T)
    pop3544 <- rowSums(cbind(e.pop.m3544, e.pop.f3544), na.rm=T)
    pop4554 <- rowSums(cbind(e.pop.m4554, e.pop.f4554), na.rm=T)
    pop5564 <- rowSums(cbind(e.pop.m5564, e.pop.f5564), na.rm=T)
    pop65 <- rowSums(cbind(e.pop.m65, e.pop.f65), na.rm=T)
    rnew0_4 <- new04 * m / pop04
    rnew5_14 <- new0514 * m / pop0514
    rnew15_24 <- new1524 * m / pop1524
    rnew25_34 <- new2534 * m / pop2534
    rnew35_44 <- new3544 * m / pop3544
    rnew45_54 <- new4554 * m / pop4554
    rnew55_64 <- new5564 * m / pop5564
    rnew65 <- new65 * m / pop65
    all.sp <- new.sp * m / e.pop.num
    all <- newinc
  })

  dta <- subset(dta, select=c(iso3, year, pop04, pop0514, pop1524, pop2534, pop3544, pop4554, pop5564, pop65,
                              new04, new0514, new1524, new2534, new3544, new4554, new5564, new65, tot.newrel, 
                              rnew0_4, rnew5_14, rnew15_24, rnew25_34, rnew35_44, rnew45_54, rnew55_64, rnew65, all, all.sp
                              ))
  dta2 <- melt(subset(dta, 
                      select=c(iso3,year,rnew0_4,rnew5_14, rnew15_24, rnew25_34, rnew35_44, rnew45_54, rnew55_64, rnew65, all)), 
               id.vars=1:2)
  dta2$variable <- gsub('rnew', '', dta2$variable)
  dta2sp <- melt(subset(dta, 
                      select=c(iso3,year,rnew0_4,rnew5_14, rnew15_24, rnew25_34, rnew35_44, rnew45_54, rnew55_64, rnew65, all.sp)), 
               id.vars=1:2)
  dta3 <- subset(dta2, !is.na(value) & value>0)
  dta3sp <- subset(dta2sp, !is.na(value) & value>0)
  
  vars <- length(table(dta3$variable))
  titsp <- ifelse(all, '(all forms)', '(smear positive)')  
  
  if (all) dta4 <- dta3
  else dta4 <- dta3sp
  
  p <- qplot(year, value, data=dta4, geom='line', colour=variable, size=variable, 
             main=paste('TB notification rates', titsp, 'by age groups in', as.character(isocty(iso)))) +
    xlab('') + 
    scale_colour_manual(values=rep('black', vars), guide='none') +
    scale_size_manual(values=c(rep(0.5, vars-1), 1.2), guide='none') +
    coord_trans(y = "log10") + ylab('Rate per 100,000 per year (log scale)') +
    theme_bw(base_size=14) 
  
  q <- direct.label(p, list('last.qp'))
  if (toplot) return(q)
  else return(dta)
}



# population pyramids - 2 versions:
# one based on plotrix::pyramid.plot, for a single pyramid
# the other based on HH::likert, for single & multiple pyramids
pyramid <- function(iso, yr=2010){
  if (missing(iso) || !is.character(iso) || nchar(iso)!=3)
    stop("Error: iso must be a valid ISO3 country code")
  if (is.na(match('pop',objects(envir=parent.frame(), all.names = TRUE))) || 
        is.na(match('cty',objects(envir=parent.frame(), all.names = TRUE))))
    stop("Error: pop and cty must be loaded")
  
  library(plotrix)
  
  iso <- toupper(iso)

  dta <- as.data.frame(pop[iso3==iso & year==yr])
  dta <- within(dta, {
    m04 <- e.pop.m04/e.pop.num
    m514 <- e.pop.m514/e.pop.num
    m1524 <- e.pop.m1524/e.pop.num
    m2534 <- e.pop.m2534/e.pop.num
    m3544 <- e.pop.m3544/e.pop.num
    m4554 <- e.pop.m4554/e.pop.num
    m5564 <- e.pop.m5564/e.pop.num
    m65 <- e.pop.m65/e.pop.num
    f04 <- e.pop.f04/e.pop.num
    f514 <- e.pop.f514/e.pop.num
    f1524 <- e.pop.f1524/e.pop.num
    f2534 <- e.pop.f2534/e.pop.num
    f3544 <- e.pop.f3544/e.pop.num
    f4554 <- e.pop.f4554/e.pop.num
    f5564 <- e.pop.f5564/e.pop.num
    f65 <- e.pop.f65/e.pop.num
  })
  country <- as.character(dta$country[1])
  mpop <- with(dta, c(m04,m514,m1524,m2534,m3544,m4554,m5564,m65) * 100)
  fpop <- with(dta, c(f04,f514,f1524,f2534,f3544,f4554,f5564,f65) * 100)
  agelab <- c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65+")
  mcol='blue'
  fcol='red'
  
  pyramid.plot(mpop, fpop, labels=agelab,
              main=paste("Population pyramid in ",country," (",yr,")", sep=""),
              lxcol=mcol, rxcol=fcol, gap=1.2, show.values=FALSE)
}

pyramids <- function(iso, yrs=c(1990,2000,2010), toplot=TRUE){
  if (missing(iso) || !is.character(iso) || nchar(iso)!=3)
    stop("Error: iso must be a valid ISO3 country code")
  if (is.na(match('pop',objects(envir=parent.frame(), all.names = TRUE))) || 
        is.na(match('cty',objects(envir=parent.frame(), all.names = TRUE))))
    stop("Error: pop and cty must be loaded")
  
  library(reshape)
  library(HH)
  
  iso <- toupper(iso)
  
  dta <- as.data.frame(pop[iso3==iso & year %in% yrs])
  dta <- within(dta, {
    m04 <- e.pop.m04/e.pop.num
    m514 <- e.pop.m514/e.pop.num
    m1524 <- e.pop.m1524/e.pop.num
    m2534 <- e.pop.m2534/e.pop.num
    m3544 <- e.pop.m3544/e.pop.num
    m4554 <- e.pop.m4554/e.pop.num
    m5564 <- e.pop.m5564/e.pop.num
    m65 <- e.pop.m65/e.pop.num
    f04 <- e.pop.f04/e.pop.num
    f514 <- e.pop.f514/e.pop.num
    f1524 <- e.pop.f1524/e.pop.num
    f2534 <- e.pop.f2534/e.pop.num
    f3544 <- e.pop.f3544/e.pop.num
    f4554 <- e.pop.f4554/e.pop.num
    f5564 <- e.pop.f5564/e.pop.num
    f65 <- e.pop.f65/e.pop.num
  })
  
  country <- as.character(dta$country[1])
  mdta <- as.data.frame(dta)[, c(2,45:52)]
  fdta <- as.data.frame(dta)[, c(2,37:44)]
  DFm <- melt(mdta, id.vars=1)  
  DFm <- DFm[order(DFm$year), ]
  DFf <- melt(fdta, id.vars=1)
  DFf <- DFf[order(DFf$year), ] 
  agelab <- c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65+")
  DFm$variable <- DFf$variable <- rep(rev(agelab), length(yrs))
  DFm$index <- DFf$index <- rep(1:8, length(yrs))
  names(DFm)[2:3] <- c('Age','Males')
  names(DFf)[2:3] <- c('Age','Females')
  
  DF <- merge(DFf, DFm, by=c('year','Age','index'))
  DFs <- melt(DF, id.vars=c('year','Age','index'))
  DFs <- DFs[order(DFs$year, DFs$variable, DFs$index), ]
  DFa <- array(c(DFs$value), dim=c(8,2,length(yrs)), dimnames=list(Age=rev(agelab), Sex=c('Females','Males'), Year=yrs))
  cond <- length(yrs)>1
  
  p <- likert(DFa, 
              main=paste("Population pyramid", ifelse(cond,'s',''),' in ',country, 
                         ifelse(cond,'',paste(' (',yrs,')',sep='')), sep=""), 
              xlab='%',
              strip.left=FALSE, strip=length(yrs)>1,
              box.width=0.95,
              auto.key=list(title=NULL),
              layout=c(length(yrs),1), between=list(x=.5)) 
  if(toplot) print(p)  
  else return(DFa)
}




# thermometers
plot.progress <- function(values, label)  {
  numOfBar <- length(values)
  plot(c(0,100), c(0,numOfBar), type='n', xlab='', ylab='', yaxt='n', mar=c(3,3,3,3))
  for(i in 1:numOfBar) {
    rect(0, 0.1+i-1, values[i], 0.9+i-1, col=rainbow(numOfBar)[i])
    text(0.5, 0.5+i-1, paste(label[i], ': ', round(values[i],2), '%', sep=''), adj=0)
  }
}





#----------------------------------------------
# incidence from prevalence
#----------------------------------------------
# function to derive prevalence from incidence
prev2inc <- function(prev, prev.se, prevk=NA, prevk.se=NA, newrel, tbhiv, tbhiv.se, nsim=50000) {
  # prev = prevalence per 100k population
  # prev.se = standard deviation of prevalence
  # prevk = prevalence of known cases (on tx)
  # prevk.se = standard deviation of prevalence of known cases (on tx)
  # newrel = case notification rate (new+relapse)
  # tbhiv = proportion HIV+ among newrel
  # nsim = number of simulation runs
  m <- 100000
  
  # durations
  Du.nh <- runif(nsim, 1, 4)     # not notified, HIV-neg
  Dn.nh <- runif(nsim, 0.2, 2)   # notified, HIV-neg
  Du.h <- runif(nsim, 0.01, 0.2) # not notified, HIV-pos
  Dn.h <- runif(nsim, 0.01, 1)   # notified, HIV-pos

  # prevalence
  get.beta <- function(ev, sd){  # beta params from moments
    S = (ev * (1 - ev) / sd^2) - 1
    a = S * ev
    b = S * (1 - ev)
    return(c(a = a, b = b))
  }
  ppar <- get.beta(prev/m, prev.se/m)
  P <- rbeta(nsim, ppar[1], ppar[2])
  
  # tbhiv
  hpar <- get.beta(tbhiv, tbhiv.se)
  H <- rbeta(nsim, hpar[1], hpar[2])
  
  # notified
  N.nh <- newrel/m * (1 - H)
  N.h <- newrel/m * H
  
  # prevalence notified
  known <- !is.na(prevk) & !is.na(prevk.se)
  if (known) {
    pnpar <- get.beta(prevk/m, prevk.se/m)
    Pn <- rbeta(nsim, pnpar[1], pnpar[2])
  } else {
    Pn <- N.nh * Dn.nh + N.h * Dn.h  # theoretical
  }
  
  # prevalence not notified
  Pu <- pmax(P - Pn, 0)
  
  # incidence not notified
  # Pu = (Iu+ d+) + (Iu- d-)
  # Iu = Iu+ + Iu- = Pu / Du
  # Du = hd+ + (1-h)d-
  Du <- H * Du.h + (1 - H) * Du.nh
  Iu <- Pu / Du
  
  # total incidence
  I <- newrel + Iu*m
  
  return(c(inc = mean(I), inc.se = sd(I)))
  
}




#----------------------------------------------
# Dataset expansion
#----------------------------------------------
expand <- function (aggregate.data, index.var = "Freq", retain.freq = FALSE) 
{
  output <- NULL
  for (i in 1:nrow(aggregate.data)) {
    if (retain.freq) {
      output <- rbind(output, 
                      aggregate.data[rep(i, aggregate.data[, 
                                                           which(names(aggregate.data) == index.var)][i]), 
                                     ])
    }
    else {
      output <- rbind(output, 
                      aggregate.data[rep(i, aggregate.data[, 
                                                           which(names(aggregate.data) == index.var)][i]), 
                                     ][, -which(names(aggregate.data) == index.var)])
    }
  }
  data.frame(output, row.names = 1:nrow(output))
}



#----------------------------------------------
# Miscellenia utilities
#----------------------------------------------

# utilities
# recode -- from car
# example:
# x<-rep(1:3,3)
# recode(x, "c(1,2)='A'; else='B'")
# [1] "A" "A" "B" "A" "A" "B" "A" "A" "B"
# recode(x, "1:2='A'; 3='B'")
# [1] "A" "A" "B" "A" "A" "B" "A" "A" "B"

recode <- function (var, recodes, as.factor.result, levels) 
{
  recode.list <- rev(strsplit(recodes, ";")[[1]])
  is.fac <- is.factor(var)
  if (missing(as.factor.result)) 
    as.factor.result <- is.fac
  if (is.fac) 
    var <- as.character(var)
  result <- var
  if (is.numeric(var)) {
    lo <- min(var, na.rm = TRUE)
    hi <- max(var, na.rm = TRUE)
  }
  for (term in recode.list) {
    if (0 < length(grep(":", term))) {
      range <- strsplit(strsplit(term, "=")[[1]][1], ":")
      low <- eval(parse(text = range[[1]][1]))
      high <- eval(parse(text = range[[1]][2]))
      target <- eval(parse(text = strsplit(term, "=")[[1]][2]))
      result[(var >= low) & (var <= high)] <- target
    }
    else if (0 < length(grep("else", term))) {
      target <- eval(parse(text = strsplit(term, "=")[[1]][2]))
      result[1:length(var)] <- target
    }
    else {
      set <- eval(parse(text = strsplit(term, "=")[[1]][1]))
      target <- eval(parse(text = strsplit(term, "=")[[1]][2]))
      for (val in set) {
        if (is.na(val)) 
          result[is.na(var)] <- target
        else result[var == val] <- target
      }
    }
  }
  if (as.factor.result) {
    result <- if (!missing(levels)) 
      factor(result, levels = levels)
    else as.factor(result)
  }
  else if (!is.numeric(result)) {
    result.valid <- na.omit(result)
    opt <- options(warn = -1)
    result.valid <- as.numeric(result.valid)
    options(opt)
    if (!any(is.na(result.valid))) 
      result <- as.numeric(result)
  }
  result
}



#-------------------------------------------------
# Compute sample size to observe at least N events
# translated from C code I wrote for sampsize utility
# http://sampsize.sourceforge.net
#-------------------------------------------------
small.size <- function(obs=10, pr=0.5, level=0.9){
  mid <- obs
  cubinom <- bot <- 0
  
  while (cubinom < level){
    cubinom <- pbeta(pr, obs, mid - obs + 1)
    mid <- mid * 2
  }
  
  
  repeat {
    if (cubinom >= level)
      top <- mid
    else 
      bot <- mid
    
    mid <- (bot + top) / 2
    cubinom <- pbeta(pr, obs, mid-obs+1)
    
    if (abs(bot - mid) <= 0.5) break  
  }
  return (ceiling(top))
}




#-------------------------------------------------
# Compute sample size based on hypergeometric distribution
#-------------------------------------------------
n.size <- function(p0=0.05, N=1000, d=1, alpha=0.05){
  s <- N
  for (n in N:1){
    m <- N - n
    k <- trunc(p0 * N)
    if(dhyper(d, n, m, k) > alpha) break
    s <- n
  }
  return(s)
}





#-------------------------------------------------
# Bayesian sample size, binomial
# translated from Fulvia Mecatti's Mathematica code
# in the Lime book
#-------------------------------------------------
bsize <- function(M=2000, m=2000, lo=10000, hi=100000, a=16, b=7504,
                  d=0.00053){
  # Vectorize rbinom() to avoid creating an inner loop 
  vrbinom <- Vectorize(rbinom, "prob")
  
  # generate M random integers between lo and hi
  N <- as.integer(runif(M, min=lo, max=hi))
  
  # outer loop 
  fnab <- numeric(M)
  for (j in 1:M){
    # generate m random values, distributed Beta(a, b)
    p <- rbeta(m, a, b)
    
    # inner loop vectorized to improve efficiency
    # generate m random values, distributed Bin(N_j, p_i)
    # stored in an M X m matrix
    x <- vrbinom(m, size=N[j], prob=p)
    k <- 1 / (x + a) + 1 / (N[j] + b - x)
    h <- (N[j] + b - x) / (x + a) + (x + a) / (N[j] + b - x)
    nab <- N[j] + a + b
    dj <- 2 / (nab * sqrt(k)) * 
      (1.96 - (1.96^3 + 3 * 1.96) * (h - 1) / (4 * nab) +
      1.96 * h / (2 * nab) +
      5 * (1.96^3 + 3 * 1.96) * (h - 2) / (18 * nab) -
      1.96 * (h - 2) / nab)
    fnab[j] <- mean(dj)
  }
  
  # OLS: 1 / f(N, a, b) = alpha1 + alpha2 * N 
  alpha <- coef(lm(I(1 / fnab^2) ~ N))
  
  size <- (1 / (4 * d^2) - alpha[1]) / alpha[2]
  names(size) <- 'N'
  
  return(ceiling(size))
}


#--------------------------------------------
# Frequentist sample size, lime book eq 9.1 and 9.2
#--------------------------------------------
ssize1 <- function(p1=0.002, p2=0.0014, m=600, k1=0.3, k2=0.1, beta=0.8){
  zb <- qnorm(beta)
  N <- (zb * sqrt(p2*(1-p2) + (m-1)*(p2*k2)^2) + 1.65*sqrt(p1*(1-p1)+(m-1)*(p1*k1)^2)) / (p1 - p2)
  N <- ceiling(N ^ 2)
  return(N)
}



ssize2 <- function(pii=0.0014, p1=0.002, m1=600, m2=600, k1=0.3, k2=0.1, N1=50000, beta=0.8){
  zb <- qnorm(beta)
  num <- (1.65 + zb)^2 * (pii * (1 - pii) + (m2-1)*pii^2*k2^2) * N1
  den <- (pii-p1)^2 * N1 - (1.65+zb)^2 * (p1*(1-p1) + (m1-1)*p1^2 * k1^2)
  N <- ceiling(num / den)
  if(N<0) N <- 'undefined'
  return(N)
}


# get k from DEFF
defftok <- function(deff=2, pii=0.002, m=600){
  k <- sqrt((deff - 1) * (1 - pii) / ((m - 1) * pii))
  return(k)
}


# get DEFF from k
ktodeff <- function(k=1, pii=0.002, m=600){
  deff <- 1 + (m - 1) * k^2 * pii / (1 - pii)
  return(deff)
}






#--------------------------------------------
# Sample size McNemar test
# Author: PG
# based on: Lehr RG. 
# Drug Information Journal 2001;35:1227-1233.
#--------------------------------------------
mcn <- function(p1, p2, alpha=0.05, power=0.8, two.tailed=TRUE){
  stopifnot(p2<p1)
  stopifnot(p1<=1 & p1>=0)
  stopifnot(p2<=1 & p2>=0)
  stopifnot(alpha>0 & alpha<1)
  stopifnot(power>0 & power<1)
  
  p <- (p1 + p2)/2
  q <- 1 - p
  d <- abs(p2 - p1)
  r <- sqrt(p2*(1-p1)/(p1*(1-p2)))
  k <- ifelse(two.tailed, 2, 1)
  m = 2 * (qnorm(1-alpha/k) + (qnorm(power)))^2
  psi <- p1 - p2
  n.lo <- m * p * q * (1 - r) / d^2
  n.hi <- m * p * q / d^2
  n.connor <- (qnorm(1-alpha/k)*sqrt(psi) + qnorm(power)*sqrt(psi-d^2))^2 / d^2
  return(list(p1=p1, p2=p2, psi=psi, r=r, n.lo=ceiling(n.lo), n.connor=ceiling(n.connor), n.hi=ceiling(n.hi)))
}








#------------------------------------------------------
# binomial CI (from from epical::ci.binomial)
#------------------------------------------------------
cii <- function (size, x, precision, alpha = 0.05) 
{
  success <- x
  if (missing(size)) {
    success1 <- success
    if (min(success, na.rm = TRUE) != 0 | max(success, na.rm = TRUE) != 
      1) {
      stop("This is not a binary vector.")
    }
    success <- length(na.omit(success1)[na.omit(success1) > 
      0])
    size <- length(na.omit(success1))
  }
  reverse <- rep(FALSE, length(success))
  reverse[success/size > 0.5] <- TRUE
  success[reverse] <- size[reverse] - success[reverse]
  if (missing(precision)) {
    precision <- success/size/10000
  }
  precision[success == 0 | success == size] <- 0.01/size[success == 
    0 | success == size]
  probab <- success/size
  success1 <- success
  success1[success > 0] <- success[success > 0] - 1
  for (i in 1:length(success)) {
    while (pbinom(success1[i], size[i], probab[i], lower.tail = FALSE) > 
      alpha/2) {
      probab[i] <- probab[i] - precision[i]
    }
  }
  estimate <- success/size
  se <- sqrt(estimate * (1 - estimate)/size)
  ll <- probab
  probab <- success/size
  for (i in 1:length(success)) {
    while (pbinom(success[i], size[i], probab[i], lower.tail = TRUE) > 
      alpha/2) {
      probab[i] <- probab[i] + precision[i]
    }
  }
  ul <- probab
  data.frame.a <- data.frame(events = success, total = size, 
                             prob = estimate, se = se, ll = ll, ul = ul)
  data.frame.a[reverse, ] <- data.frame(events = size[reverse] - 
    success[reverse], total = size[reverse], prob = 1 - 
    estimate[reverse], se = se[reverse], ll = 1 - ul[reverse], 
                                        ul = 1 - ll[reverse])
  names(data.frame.a)[5] <- paste("lower", 100 * (1 - 
    alpha), "ci", sep = "")
  names(data.frame.a)[6] <- paste("upper", 100 * (1 - 
    alpha), "ci", sep = "")
  if (nrow(data.frame.a) == 1) {
    rownames(data.frame.a) <- ""
  }
  data.frame.a
}






#------------------------------------------------------
# PPV, PPN
#------------------------------------------------------
ppv <- function(prev=0.15, se=0.95, sp=0.98) se*prev / (se * prev + (1 - sp) * (1 - prev))
ppn <- function(prev=0.15, se=0.95, sp=0.98) sp*(1-prev) / (sp*(1-prev) + (1-se)*prev)


# PAF
paf <- function(pe, rr) pe * (rr - 1) / (pe * (rr - 1) + 1)





#------------------------------------------------------
# capture recapture 3 lists, loglinear
#
# A: count in A only (not in B, not in C)
# AB: count in A and B not in C
# ABC: count in A and B and C
# deps = 'AIC','independent','AB','AC','BC','ABAC','ABBC','ACBC','saturated'
#
# updated on 11 Aug 2011
#------------------------------------------------------
capture <- function (A, B, C, AB, AC, BC, ABC, deps='saturated', level=0.05){
  library(MASS)
  
  lista <- c(0,1,1,0,0,1,0)
  listb <- c(1,0,1,0,1,0,0)
  listc <- c(1,1,0,1,0,0,0)
  
  frq <- c(A, B, C, AB, AC, BC, ABC)
  
  dta <- data.frame(frq, lista, listb, listc)
  if (deps=='AIC') fit <- stepAIC(
    glm(frq ~ lista * listb * listc, family=poisson, data=dta),
    direction = ('backward'), trace=0)
  else if (deps=='independent') 
    fit <- glm(frq ~ lista + listb + listc, family=poisson, data=dta)
  else if (deps=='AB')
    fit <- glm(frq ~ lista + listb + listc + lista:listb, family=poisson, data=dta)
  else if (deps=='AC')
    fit <- glm(frq ~ lista + listb + listc + lista:listc, family=poisson, data=dta)
  else if (deps=='BC')
    fit <- glm(frq ~ lista + listb + listc + lista:listb, family=poisson, data=dta)
  else if (deps=='ABAC')
    fit <- glm(frq ~ lista + listb + listc + lista:listb + lista:listc, family=poisson, data=dta)
  else if (deps=='ABBC')
    fit <- glm(frq ~ lista + listb + listc + lista:listb + listb:listc, family=poisson, data=dta)
  else if (deps=='ACBC')
    fit <- glm(frq ~ lista + listb + listc + lista:listc + listb:listc, family=poisson, data=dta)
  else if (deps=='saturated') 
    fit <- glm(frq ~  lista + listb + listc + lista:listb + lista:listc + listb:listc,
               family=poisson, data=dta)
  
  n <- sum(frq)
  x <- sum(coef (fit))
  se <- sqrt(sum(vcov(fit)))
  low <- x + qnorm(level/2) * se
  high <- x + qnorm(1 - level/2) * se
  
  n.point <- round(n + exp(x))
  n.low <- floor(n + exp(low))
  n.high <- ceiling(n + exp(high))
  
  return(list(fit=summary(fit), 
              estimated=c(n=n.point, lo=n.low, hi=n.high), 
              listed=c(A=sum(A,AB,AC,ABC), B=sum(B,AB,BC,ABC), C=sum(C,AC,BC,ABC)), 
              inventory=n,
              added=c(n=n.point - n, lo=n.low - n, hi=n.high - n)
  )
  )
}





#------------------------------------------------------
# util
#------------------------------------------------------
# long print of data.table x
lp <- function(x)print(x, nrow=Inf)

# get rate of change from two data points
rc <- function(v, time, corr=0.5)
  return(rate.change = (log(v[2] + corr) - log(v[1] + corr))/(time[2] - time[1]))

# get SD of rate of change from two data points
rc.sd <- function(v.sd, time)
  return(rate.change.sd = sqrt(v.sd[2]^2 + v.sd[1]^2)/(time[2] - time[1]))

# not in operator
`%ni%` <- Negate(`%in%`) 




#------------------------------------------------------
# logistic function
#------------------------------------------------------
lgst <- function(a, l, k, d, t) a + (l - a)/(1 + exp((k-t)/d)) 
# a = lower asymptote
# l = asymptote at infinite
# k = growth
# d = time to mid-growth
# t = time
#------------------------------------------------------




#------------------------------------------------------
# double logistic function 
# (http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2576736/)
#------------------------------------------------------
dlgst <- function(a, r, t0, l, d, t) exp(r*(t-t0))/(1 + exp(r*(t-t0))) * 
  (2 * a*(exp(-(d*(t-t0)))) / (1 + exp(-(d*(t-t0)))) + l)
# a = peak
# r = rate of increase
# t0 = peak time
# d = rate of decline
# t = time
#------------------------------------------------------





#------------------------------------------------------
# ggplto curve - simple wrapper
#------------------------------------------------------
gcurve <- function(from, to, fun, col='black', ...) 
  qplot(c(from, to), geom='density', ...) +
  stat_function(fun=fun, colour=I(col)) + 
  xlab('x')






#------------------------------------------------------
# Robust standard errors
#------------------------------------------------------
summaryR.lm <- function(model, type=c("hc3", "hc0", "hc1", "hc2", "hc4"), ...){
  
  if (!library(car)) stop("libraryd car package is missing.")
  
  type <- match.arg(type)
  V <- hccm(model, type=type)
  sumry <- summary(model)
  table <- coef(sumry)
  table[,2] <- sqrt(diag(V))
  table[,3] <- table[,1]/table[,2]
  table[,4] <- 2*pt(abs(table[,3]), df.residual(model), lower.tail=FALSE)
  
  sumry$coefficients <- table
  p <- nrow(table)
  hyp <- cbind(0, diag(p - 1))
  sumry$fstatistic[1] <- linearHypothesis(model, hyp,white.adjust=type)[2,"F"]
  
  return(sumry)
  cat("Note: Heteroscedasticity-consistent standard errors using adjustment", type, "\n")
}
#------------------------------------------------------





#------------------------------------------------------
# adjust x-axis labels in ggplot2
#------------------------------------------------------
facetAdjust <- function(x, pos = c("up", "down")){
  pos <- match.arg(pos)
  p <- ggplot_build(x)
  gtable <- ggplot_gtable(p); dev.off()
  dims <- apply(p$panel$layout[2:3], 2, max)
  nrow <- dims[1]
  ncol <- dims[2]
  panels <- sum(grepl("panel", names(gtable$grobs)))
  space <- ncol * nrow
  n <- space - panels
  if(panels != space){
    idx <- (space - ncol - n + 1):(space - ncol)
    gtable$grobs[paste0("axis_b",idx)] <- list(gtable$grobs[[paste0("axis_b",panels)]])
    if(pos == "down"){
      rows <- grep(paste0("axis_b\\-[", idx[1], "-", idx[n], "]"), 
                   gtable$layout$name)
      lastAxis <- grep(paste0("axis_b\\-", panels), gtable$layout$name)
      gtable$layout[rows, c("t","b")] <- gtable$layout[lastAxis, c("t")]
    }
  }
  class(gtable) <- c("facetAdjust", "gtable", "ggplot"); gtable
}
# The function for printing which differs only by few lines from ggplot2:::print.ggplot:
print.facetAdjust <- function(x, newpage = is.null(vp), vp = NULL) {
  if(newpage)
    grid.newpage()
  if(is.null(vp)){
    grid.draw(x)
  } else {
    if (is.character(vp)) 
      seekViewport(vp)
    else pushViewport(vp)
    grid.draw(x)
    upViewport()
  }
  invisible(x)
}
#------------------------------------------------------





#------------------------------------------------------
# sample from an arbitrary pdf
#------------------------------------------------------
gicdf <- function(fun, 
                  min=-3.5, 
                  max=3.5, 
                  bins=1000, 
                  nqratio=10, 
                  grouping=mean, 
                  ...) {
  # Generate an inverse CDF of an arbitrary function
  fun <- match.fun(fun)
  grouping <- match.fun(grouping)
  
  # Number of points to draw
  nq=nqratio*bins
  
  # Draw indexes
  qdraw <- seq(min, max,length.out=nq)
  
  # Calculate proportional probability of each draw
  pdraw <- fun(qdraw,...)
  
  # Rescale probability sum to 1, rescale
  pdraw <- pdraw/sum(pdraw)
  
  # Calculate the cumulative probability at each qdraw
  cpdraw <- cumsum(pdraw)
  
  # Caculate the cumulative probability at each bin
  pbin <- (1:bins)/bins
  xbin <- NA*(1:bins)
  
  for (i in 1:bins) {
    xbin[i] <- grouping(qdraw[cpdraw<pbin[i]&cpdraw>0], na.rm = TRUE)
    cpdraw[cpdraw<pbin[i]] <- 2
  }
  
  (draw.set <- list(digits=floor(log10(bins)), xbin=xbin, pbin=pbin))
}

# Draw from acdf
ricdf <- function(N, draw.set) {
  digits <- draw.set$digits
  pdraws <- ceiling(runif(N)*10^digits)/10^digits
  draw.set$xbin[match(pdraws,draw.set$pbin)]
}

# Define an arbitrary pdf f: normal pdf
# f <- function(x) dnorm(x)
# 
# # set the number to sample
# nsamples <- 10000
# 
# # Draw from rnorm
# sample0 <- cbind(draw=rnorm(nsamples), type=0)
# 
# # Draw inverted cdf distribution information for range between -4 and 4
# # using mean approach
# sample1 <- cbind(draw=ricdf(nsamples, gicdf(f,min=-4,max=4)), type=1)
#------------------------------------------------------





#------------------------------------------------------
# reshape data.table
# from: https://gist.github.com/mrdwab/11380733
# source_gist(11380733) 
#------------------------------------------------------
cSplit <- function(indt, splitCols, sep = ",", direction = "wide", 
                   makeEqual = NULL, fixed = TRUE) {
  ## requires data.table >= 1.8.11
  require(data.table)
  if (!is.data.table(indt)) setDT(indt)
  if (any(!vapply(indt[, splitCols, with = FALSE],
                  is.character, logical(1L)))) {
    indt[, eval(splitCols) := lapply(.SD, as.character),
         .SDcols = splitCols]
  }
  X <- lapply(indt[, splitCols, with = FALSE], function(x) {
    strsplit(x, split = sep, fixed = fixed)
  })
  if (direction == "long") {
    if (is.null(makeEqual)) {
      IV <- function(x,y) if (identical(x,y)) TRUE else FALSE
      makeEqual <- ifelse(Reduce(IV, rapply(X, length, how = "list")),
                          FALSE, TRUE)
    }
  } else if (direction == "wide") {
    if (!is.null(makeEqual)) {
      if (!isTRUE(makeEqual)) {
        message("makeEqual specified as FALSE but set to TRUE")
        makeEqual <- TRUE
      }
      makeEqual <- TRUE
    } else {
      makeEqual <- TRUE
    }
  }
  if (isTRUE(makeEqual)) {
    SetUp <- lapply(seq_along(X), function(y) {
      A <- vapply(X[[y]], length, 1L)
      list(Mat = cbind(rep(seq_along(A), A), sequence(A)),
           Val = unlist(X[[y]]))
    })    
    Ncol <- max(unlist(lapply(SetUp, function(y) y[["Mat"]][, 2]), 
                       use.names = FALSE))
    X <- lapply(seq_along(SetUp), function(y) {
      M <- matrix(NA_character_, nrow = nrow(indt), ncol = Ncol)
      M[SetUp[[y]][["Mat"]]] <- SetUp[[y]][["Val"]]
      M
    })
    if (direction == "wide") {
      X <- lapply(seq_along(X), function(x) {
        colnames(X[[x]]) <- paste(splitCols[x], 
                                  sequence(ncol(X[[x]])), 
                                  sep = "_")
        X[[x]]
      })
      cbind(indt, do.call(cbind, X))[, eval(splitCols) := NULL][]
    } else {
      indt <- indt[rep(sequence(nrow(indt)), each = Ncol)]
      X <- lapply(X, function(y) as.vector(t(y)))
      indt[, eval(splitCols) := lapply(X, unlist, use.names = FALSE)][]
    }  
  } else {
    Rep <- vapply(X[[1]], length, integer(1L))
    indt <- indt[rep(sequence(nrow(indt)), Rep)]
    indt[, eval(splitCols) := lapply(X, unlist, use.names = FALSE)][]
  }
}






#------------------------------------------------------
# bootstrap df with tidy model objects
# from: https://github.com/hadley/dplyr/issues/269
#------------------------------------------------------
bootstrap <- function(df, m) {
  n <- nrow(df)
  
  attr(df, "indices") <- replicate(m, sample(n, replace = TRUE), 
                                   simplify = FALSE)
  attr(df, "drop") <- TRUE
  attr(df, "group_sizes") <- rep(n, m)
  attr(df, "biggest_group_size") <- n
  attr(df, "labels") <- data.frame(replicate = 1:m)
  attr(df, "vars") <- list(quote(replicate)) # Change
  class(df) <- c("grouped_bootstrap_df", "grouped_df", "tbl_df", "tbl", "data.frame")
  
  df
}

#library(dplyr)
#bootnls <- bootstrap(mtcars, 100) %>% do(tidy(nls(mpg ~ k / wt + b, ., start=list(k=1, b=0)))) 
#bootnls %>% group_by(term) %>% summarize(low=quantile(estimate, 0.025), high=quantile(estimate, 0.975))







#------------------------------------------------------
# lives saved
#------------------------------------------------------
cty.lsaved <- function(iso='CHN', start=1995, csv=FALSE){
  m <- 1e5
  M <- 1e6
  
  lsaved <- function(dta){
    # CFRs untreated
    # HIV negative not on TB treatment  0.43 (0.28 - 0.53), assume ~beta
    cfrn <- 0.43
    cfrn.se <- (0.53 - 0.28)/4
    # HIV positive not on ART, not on TB treatment  0.78 (0.65 - 0.94), assume ~beta
    cfrp <- 0.78
    cfrp.se <- (0.94 - 0.65)/4
    
    # lives saved HIV-neg, HIV-pos and total by row
    for (i in 1:dim(dta)[1]){ # ugly loop, will optimize later
      # counterfactual (rates)
      cfn <- betaop(dta$inc.nh[i]/m, cfrn, dta$inc.nh.se[i]/m, cfrn.se, op='*', dist=T, nsim=M)
      cfp <- betaop(dta$inc.h[i]/m, cfrp, dta$inc.h.se[i]/m, cfrp.se, op='*', dist=T, nsim=M)
      
      # factual (rates)
      parn <- get.beta(dta$mort.nh[i]/m, dta$mort.nh.se[i]/m)
      fn <- rbeta(M, parn[1], parn[2])
      parp <- get.beta(dta$mort.h[i]/m, dta$mort.h.se[i]/m)
      fp <- rbeta(M, parp[1], parp[2])
      
      # difference counterfactual minus factual (absolute numbers)
      savedn <- pmax((cfn - fn), 0) * dta$e.pop.num[i]
      savedp <- pmax((cfp - fp), 0) * dta$e.pop.num[i]
      saved <- savedn + savedp
      
      dta$savedn[i] <- mean(savedn) 
      dta$savedn.se[i] <- sd(savedn) 
      dta$savedn.lo[i] <- quantile(savedn, probs=0.025)
      dta$savedn.hi[i] <- quantile(savedn, probs=0.975) 
      
      dta$savedp[i] <- mean(savedp)
      dta$savedp.se[i] <- sd(savedp) 
      dta$savedp.lo[i] <- quantile(savedp, probs=0.025) 
      dta$savedp.hi[i] <- quantile(savedp, probs=0.975)
      
      dta$saved[i] <- mean(saved) 
      dta$saved.se[i] <- sd(saved) 
      dta$saved.lo[i] <- quantile(saved, probs=0.025) 
      dta$saved.hi[i] <- quantile(saved, probs=0.975) 
    }
    
    return(dta)
  }
  
  clsaved <- function(dta){
    # cumulative lives saved
    savedn <- sum(dta$savedn)
    savedn.se <- sqrt(sum(dta$savedn.se^2))
    savedn.lo <- savedn - 1.96 * savedn.se
    savedn.hi <- savedn + 1.96 * savedn.se
    
    savedp <- sum(dta$savedp)
    savedp.se <- sqrt(sum(dta$savedp.se^2))
    savedp.lo <- savedp - 1.96 * savedp.se
    savedp.hi <- savedp + 1.96 * savedp.se
    
    saved <- sum(dta$saved)
    saved.se <- sqrt(sum(dta$saved.se^2))
    saved.lo <- saved - 1.96 * saved.se
    saved.hi <- saved + 1.96 * saved.se
    dta <- data.table(savedn, savedn.lo, savedn.hi, savedn.se,
                      savedp, savedp.lo, savedp.hi, savedp.se,
                      saved, saved.lo, saved.hi, saved.se)
    return(dta)
  }
  
  th <- 1000  
  bgd <- lsaved(est[iso3==iso & year>=start])
  bgd.cum <- clsaved(bgd)
  bgd.saved <- bgd[, list(year=as.character(year), savedn, savedn.lo, savedn.hi, savedn.se,
                          savedp, savedp.lo, savedp.hi, savedp.se,
                          saved, saved.lo, saved.hi, saved.se)]
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



