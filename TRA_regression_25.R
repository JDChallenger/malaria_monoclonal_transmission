library(rstan)
library(rethinking) #See README for installation instructions
library(ggplot2)
library(reshape2)

#in 'ViewData.R' I created a single dataframe 'datY' containing all 8 data sets, in the same format
#Load the data
source('ViewData.R')
#View(dat25)

#Data definitions: 'conc' = antibody concentration, 'ID' = experiment or donor ID, 
#'ratio' = ratio of oocyst counts between control group#& intervention group
#'logratio' = log10 (ratio). 'assay' = 1 if assay = DMFA. 'country' = 1 if country = Cameroon
#'sqconc' = square root of concentration

####################################################################################################  
#############################              Pfs25                 ###################################
####################################################################################################

iterz <- 10000

###     M1me: Single slope (with mixed effects)      ###
m1me_25 <- map2stan(
  alist(
    logratio ~ dnorm(mu , sigma),
    mu <-  sqconc*(a + ae[ID]),
    ae[ID] ~ dnorm(0, sigma2),
    a ~ dnorm(0,1), # Prior on slope
    sigma ~ dexp(1), # Prior on sigma
    sigma2 ~ dexp(1) # Prior on sigma
  ),
  data = dat25, iter=iterz , chains=3, cores = 3
)
precis(m1me_25,prob=0.9) # Model summary

###     M2me: Different slopes for the two assays (and mixed effects)    ###
m2me_25 <- map2stan(
  alist(
    logratio ~ dnorm( mu , sigma),
    mu <- sqconc*(a + ae[ID] + b*assay),
    ae[ID] ~ dnorm(0,sigma2), #
    a ~ dnorm(0,1), # Prior on slope
    b ~ dnorm(0,1), # Prior on DMFA slope correction
    sigma ~ dexp(1), # Prior on sigma
    sigma2 ~ dexp(1) # Prior on sigma
  ),
  data = dat25, iter=iterz , chains=3, cores = 3
)
precis(m2me_25,prob=0.9) # Model summary

###     M5me: Different slopes for the two assays, and then also for the two countries (with mixed-effects)   ###
m5me_25 <- map2stan(
  alist(
    logratio ~ dnorm( mu , sigma),
    mu <- sqconc*(a + ae[ID] + b*assay + d*country),
    ae[ID] ~ dnorm(0,sigma2), 
    a ~ dnorm(0,1), # Prior on slope
    b ~ dnorm(0,1), # Prior on DMFA slope correction
    d ~ dnorm(0,1), # Prior on country slope correction
    sigma ~ dexp(1), # Prior on sigma
    sigma2 ~ dexp(1)
  ),
  data = dat25, iter=iterz , chains=3, cores = 3
)
precis(m5me_25,prob=0.9) # Model summary

###     M6: Different slopes for the two countries, not the two assays (with mixed effects)   ###
m6me_25 <- map2stan(
  alist(
    logratio ~ dnorm( mu , sigma),
    mu <- sqconc*(a + ae[ID] + d*country),
    ae[ID] ~ dnorm(0,sigma2), 
    a ~ dnorm(0,1), # Prior on slope
    d ~ dnorm(0,1), # Prior on country slope correction
    sigma ~ dexp(1), # Prior on sigma
    sigma2 ~ dexp(1)
  ),
  data = dat25, iter=iterz , chains=3, cores = 3
)
precis(m6me_25,prob=0.9) # Model summary


#Compare all mixed effects models. (Supplementary Table 5)
compare(m1me_25,m2me_25,m5me_25,m6me_25)


####################################################################################################  
##############     Pfs25: Make Plots and generate interval for $IC_{80}$       #####################
####################################################################################################

#Four models contribute significant WAIC weight.
cf.seq <- seq(0.1,10,0.01)
l1 <- length(cf.seq)
# 10000 samples
ae_zeros <- matrix(0,10000,42)
pred1 <- list( #Fake data for SMFA?
  sqconc = cf.seq,
  logratio = rep(0,l1),#empty outcome
  ratio = rep(0,l1),#empty outcome
  assay = rep(0,l1), #
  country = rep(0,l1), #
  ID = factor(rep("test",l1)), # Is this ok?
  antibody = rep("25",l1),
  conc = cf.seq*cf.seq
)

ensemble25 <- ensemble(m1me_25,m6me_25,m2me_25,m5me_25, data = pred1, n=10000, replace = list(ae = ae_zeros))
mu <- apply(ensemble25$link, 2, mean)
mu.PI <- apply(ensemble25$link, 2, PI, prob=0.90)

mudf1 <- data.frame(sqconc = cf.seq, mu = mu)
PIdf1 <- data.frame(sqconc = cf.seq, PI1 = mu.PI[1,], PI2 = mu.PI[2,])

#SMFA only
#ggplot() + geom_point(data = dat25[which(dat25$assay==0),], aes(x = sqconc, y = logratio)) + themeJDC +
#  geom_line(data = mudf1, aes(x=sqconc, y=mu), colour='blue', alpha=.2) +
#  geom_ribbon(data = PIdf1, aes(x = sqconc, ymin = PI1, ymax = PI2), fill = "blue", alpha = 0.1) +
#  annotate("segment", x = 0.0, xend = 9, y = log10(5), yend = log10(5), linetype = 'dashed')

#Whack it back on original scale
mudf1x <- data.frame(conc = cf.seq*cf.seq, mu = 10**mu)
PIdf1x <- data.frame(conc = cf.seq*cf.seq, PI1 = 10**mu.PI[1,], PI2 = 10**mu.PI[2,])

#Pick out
m1x <- 0
for(i in 1:(length(mu))){
  if(mudf1x$mu[i] > 5){
    print("Mean IC80:")
    print(mudf1x$conc[i])
    m1x <- round(PIdf1x$conc[i],1)
    break
  }
}

#Pick out
m1y <- 0
for(i in 1:(length(mu))){
  if(PIdf1x$PI1[i] > 5){
    print("CI1:")
    print(PIdf1x$conc[i])
    m1y <- round(PIdf1x$conc[i],1)
    break
  }
}

#Pick out
m1z <- 0
for(i in 1:(length(mu))){
  if(PIdf1x$PI2[i] > 5){
    print("CI2:")
    print(PIdf1x$conc[i])
    m1z <- round(PIdf1x$conc[i],1)
    break
  }
}

pl1 <- ggplot() + geom_point(data = dat25[which(dat25$assay==0),], aes(x = conc, y = ratio), color='darkgreen',alpha=0.5) + 
  geom_line(data = dat25[which(dat25$assay==0),], aes(x = conc, y = ratio, group = ID), alpha=.3) +
  scale_x_continuous(trans='sqrt', limits = c(0.5,95), breaks = c(1,4,10,23,50,94)) + 
  scale_y_continuous(trans='log10',limits = c(0.25,240),
  sec.axis = sec_axis(~.*1, name = "TRA (%)", breaks = c(1,5,20,100), labels = c(0,80,95,99))) + themeJDC +
  geom_line(data = mudf1x, aes(x=conc, y=mu), colour='darkgreen', alpha=.75) +
  geom_ribbon(data = PIdf1x, aes(x = conc, ymin = PI1, ymax = PI2), fill = "darkgreen", alpha = 0.3) +
  annotate("segment", x = 0.5, xend = 95, y = 5, yend = 5, linetype = 'dashed') +
  xlab("mAb conc. [\U03BCg/mL]") + ylab("Oocyst ratio(control/test)") + 
  annotate("text", x = 42, y = 0.4, label = deparse(bquote("IC"[80]~"="~.(m1x)~"\U03BCg/mL ["~.(m1z)~","~.(m1y)~"]")),
           colour = "darkgreen", size = 5,parse = T) 
pl1
#ggsave(filename = "whatever1.pdf",height = 5, width = 6,  device = cairo_pdf)

# 10000 samples
ae_zeros <- matrix(0,10000,42)
pred2 <- list( #Fake data for DMFA? BF
  sqconc = cf.seq,
  logratio = rep(0,l1),#empty outcome
  ratio = rep(0,l1),#empty outcome
  assay = rep(1,l1), #
  country = rep(0,l1), #
  ID = factor(rep("test",l1)), # Is this ok?
  antibody = rep("25",l1),
  conc = cf.seq*cf.seq
)
ensemble25_2 <- ensemble(m1me_25,m6me_25,m2me_25,m5me_25, data = pred2, n=10000, replace = list(ae = ae_zeros))
mu_2 <- apply(ensemble25_2$link, 2, mean)
mu.PI_2 <- apply(ensemble25_2$link, 2, PI, prob=0.95)

mudf2x <- data.frame(conc = cf.seq*cf.seq, mu = 10**mu_2)
PIdf2x <- data.frame(conc = cf.seq*cf.seq, PI1 = 10**mu.PI_2[1,], PI2 = 10**mu.PI_2[2,])

#Pick out
m2x <- 0
for(i in 1:(length(mu))){
  if(mudf2x$mu[i] > 5){
    print("Mean IC80:")
    print(mudf2x$conc[i])
    m2x <- round(PIdf2x$conc[i],1)
    break
  }
}

#Pick out
m2y <- 0
for(i in 1:(length(mu))){
  if(PIdf2x$PI1[i] > 5){
    print("CI1:")
    print(PIdf2x$conc[i])
    m2y <- round(PIdf2x$conc[i],1)
    break
  }
}

#Pick out
m2z <- 0
for(i in 1:(length(mu))){
  if(PIdf2x$PI2[i] > 5){
    print("CI2:")
    print(PIdf2x$conc[i])
    m2z <- round(PIdf2x$conc[i],1)
    break
  }
}

# 10000 samples 
ae_zeros <- matrix(0,10000,42)
pred3 <- list( #Synthetic data for DMFA? CAM
  sqconc = cf.seq,
  logratio = rep(0,l1),#empty outcome
  ratio = rep(0,l1),#empty outcome
  assay = rep(1,l1), #
  country = rep(1,l1), #
  ID = factor(rep("test",l1)), # Is this ok?
  antibody = rep("25",l1),
  conc = cf.seq*cf.seq
)
ensemble25_3 <- ensemble(m1me_25,m6me_25,m2me_25,m5me_25, data = pred3, n=10000, replace = list(ae = ae_zeros))

mu_3 <- apply(ensemble25_3$link, 2, mean)
mu.PI_3 <- apply(ensemble25_3$link, 2, PI, prob=0.92)

mudf3x <- data.frame(conc = cf.seq*cf.seq, mu = 10**mu_3)
PIdf3x <- data.frame(conc = cf.seq*cf.seq, PI1 = 10**mu.PI_3[1,], PI2 = 10**mu.PI_3[2,])

#Pick out
m3x <- 0
for(i in 1:(length(mu))){
  if(mudf3x$mu[i] > 5){
    print("Mean IC80:")
    print(mudf3x$conc[i])
    m3x <- round(mudf3x$conc[i],1)
    break
  }
}

#Pick out
m3y <- 0
for(i in 1:(length(mu))){
  if(PIdf3x$PI1[i] > 5){
    print("CI1:")
    print(PIdf3x$conc[i])
    m3y <- round(PIdf3x$conc[i],1)
    break
  }
}

#Pick out
m3z <- 0
for(i in 1:(length(mu))){
  if(PIdf3x$PI2[i] > 5){
    print("CI2:")
    print(PIdf3x$conc[i])
    m3z <- round(PIdf3x$conc[i],1)
    break
  }
}

#DMFA only
pl2 <- ggplot() + geom_point(data = dat25[which(dat25$assay==1),], aes(x = conc, y = ratio, color=factor(country)), alpha=.5) +
  geom_line(data = dat25[which(dat25$assay==1),], aes(x = conc, y = ratio, group = ID), alpha=.15) +
  scale_x_continuous(trans='sqrt', limits = c(0.5,95), breaks = c(1,6,10,31,60,94)) + 
  scale_y_continuous(trans='log10',limits = c(0.25,240), 
  sec.axis = sec_axis(~.*1, name = "TRA (%)", breaks = c(1,5,20,100), labels = c(0,80,95,99))) + themeJDC +
  guides(colour=guide_legend(title = "Country")) + scale_color_discrete(labels = c("BF","CAM")) +
  geom_line(data = mudf2x, aes(x=conc, y=mu), colour='red', alpha=.75) +
  geom_ribbon(data = PIdf2x, aes(x = conc, ymin = PI1, ymax = PI2), fill = "red", alpha = 0.25) +
  geom_line(data = mudf3x, aes(x=conc, y=mu), colour='turquoise', alpha=.75) +
  geom_ribbon(data = PIdf3x, aes(x = conc, ymin = PI1, ymax = PI2), fill = "turquoise", alpha = 0.25) +
  annotate("segment", x = 0.5, xend = 95, y = 5, yend = 5, linetype = 'dashed') +
  xlab(paste0("mAb conc [\U03BCg/mL]")) +
  theme(legend.position = c(0.1,0.85)) + ylab("Oocyst ratio(control/test)") + 
  annotate("text", x = 46, y=0.32, label = deparse(bquote("IC"[80]~"="~.(m2x)~"\U03BCg/mL ["~.(m2z)~","~.(m2y)~"]")),
           colour = "red", size = 5,parse = T) +
  annotate("text", x = 46, y=0.57, label = deparse(bquote("IC"[80]~"="~.(m3x)~"\U03BCg/mL ["~.(m3z)~","~.(m3y)~"]")),
           colour = "turquoise", size = 5,parse = T)
pl2
#ggsave(filename = "whatever2.pdf", height = 5, width = 6, device = cairo_pdf)

