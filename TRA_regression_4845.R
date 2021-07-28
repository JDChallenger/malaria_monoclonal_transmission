library(rstan)
library(rethinking) #See README for installation instructions
library(ggplot2)
library(reshape2)

#in 'ViewData.R' I created a single dataframe 'datY' containing all 8 data sets, in the same format
##Load the data
source("ViewData.R")
#View(dat4845)

#Data definitions: 'conc' = antibody concentration, 'ID' = experiment or donor ID, 
#'ratio' = ratio of oocyst counts between control group#& intervention group
#'logratio' = log10 (ratio). 'assay' = 1 if assay = DMFA. 'country' = 1 if country = Cameroon
#'sqconc' = square root of concentration

####################################################################################################  
######################          Drop certain antibody concentrations           #####################
####################################################################################################

dat4845Z <- dat4845[dat4845$conc >12 | dat4845$conc <4.9,]
dat4845Y <- dat4845Z[dat4845Z$conc < 19,]

####################################################################################################  
#############################              Pfs4845               ###################################
####################################################################################################

iterz <- 8000

###     M1me: Single slope (with mixed effects)      ###
m1me_4845 <- map2stan(
  alist(
    logratio ~ dnorm(mu, sigma),
    mu <-  sqconc*(a + ae[ID]),
    ae[ID] ~ dnorm(0, sigma2),
    a ~ dnorm(0,1), # Prior on slope
    sigma ~ dexp(1), # Prior on sigma
    sigma2 ~ dexp(1) # Prior on sigma2
  ),
  data = dat4845Y, iter=iterz , chains=3, cores = 3
)
precis(m1me_4845,prob=0.9) # Model summary

###     M2: Different slopes for the two assays     ###
m2me_4845 <- map2stan(
  alist(
    logratio ~ dnorm(mu, sigma1),
    mu <- (a + ae[ID] + b*assay)*sqconc ,
    ae[ID] ~ dnorm(0, sigma2),
    a ~ dnorm(0,1), # Prior on slope
    b ~ dnorm(0,1), # Prior on DMFA slope correction
    sigma1 ~ dcauchy(0,1), # Prior on sigma1
    sigma2 ~ dcauchy(0,1) # Prior on sigma2
  ),
  data = dat4845Y, iter=2*iterz , chains=3, cores = 3
)
precis(m2me_4845,prob=0.9) # Model summary

###     M5me: Different slopes for the two assays, and then also for the two countries (with mixed-effects)   ###
m5me_4845 <- map2stan(
  alist(
    logratio ~ dnorm( mu , sigma),
    mu <- sqconc*(a + ae[ID] + b*assay + d*country),
    ae[ID] ~ dnorm(0,sigma2), # Prior on slope
    a ~ dnorm(0,1),
    b ~ dnorm(0,1), # Prior on DMFA slope correction
    d ~ dnorm(0,1), # Prior on country slope correction
    sigma ~ dcauchy(0,1), # Prior on sigma
    sigma2 ~ dcauchy(0,1)
  ),
  data = dat4845Y, iter=3*iterz , chains=4, cores = 4, control = list(adapt_delta = 0.94, max_treedepth = 15)
)
precis(m5me_4845,prob=0.9) # Model summary

###     M6: Different slopes for the two countries, not the two assays (with mixed effects)   ###
m6me_4845 <- map2stan(
  alist(
    logratio ~ dnorm( mu , sigma),
    mu <- sqconc*(a + ae[ID] + d*country),
    ae[ID] ~ dnorm(0,sigma2), # Prior on slope
    a ~ dnorm(0,1),
    d ~ dnorm(0,1), # Prior on country slope correction
    sigma ~ dcauchy(0,1), # Prior on sigma
    sigma2 ~ dcauchy(0,1)
  ),
  data = dat4845Y, iter=3*iterz , chains=6, cores = 6, control = list(adapt_delta = 0.92, max_treedepth = 14)
)
precis(m6me_4845,prob=0.9) # Model summary

#Compare all mixed effects models. (Supplementary Table 5)
compare(m1me_4845,m2me_4845,m5me_4845,m6me_4845)

####################################################################################################  
##############     Pfs4845: Make Plots and generate interval for $IC_{80}$       ###################
####################################################################################################

#Four models contribute significant WAIC weight.
cf.seq <- seq(0.1,4,0.01)
# 5000 samples 
ae_zeros <- matrix(0,5000,62)
l1 <- length(cf.seq)
pred1 <- list( #Fake data for SMFA?
  sqconc = cf.seq,
  logratio = rep(0,l1),#empty outcome
  ratio = rep(0,l1),#empty outcome
  assay = rep(0,l1), #
  country = rep(0,l1), #
  ID = factor(rep("test",l1)), # Is this ok?
  antibody = rep("4845",l1),
  conc = cf.seq*cf.seq
)

ensemble4845 <- ensemble(m5me_4845, m6me_4845, m2me_4845, m1me_4845, data = pred1,
                          replace = list(ae = ae_zeros))
mu <- apply(ensemble4845$link, 2, mean)
mu.PI <- apply(ensemble4845$link, 2, PI,prob=0.95)

mudf1 <- data.frame(sqconc = cf.seq, mu = mu)
PIdf1 <- data.frame(sqconc = cf.seq, PI1 = mu.PI[1,], PI2 = mu.PI[2,])

#Whack it back on original scale
mudf1x <- data.frame(conc = cf.seq*cf.seq, mu = 10**mu)
PIdf1x <- data.frame(conc = cf.seq*cf.seq, PI1 = 10**mu.PI[1,], PI2 = 10**mu.PI[2,])

#Pick out
m1x <- 0
for(i in 1:(length(mu))){
  if(mudf1x$mu[i] > 5){
    print("Mean IC80:")
    print(mudf1x$conc[i])
    m1x <- round(mudf1x$conc[i],1)
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

#Now plotting datY, with additional points as open circles
pl3 <- ggplot() + geom_point(data = dat4845Y[which(dat4845Y$assay==0),], aes(x = conc, y = ratio), colour='darkgreen', alpha =.7) + 
  scale_x_continuous(trans='sqrt', limits = c(0.02,15.1), breaks = c(0.04,0.6,1.25,2.5,5,10,15)) +
  scale_y_continuous(trans='log10',limits = c(0.25,330),
  sec.axis = sec_axis(~.*1, name = "TRA (%)", breaks = c(1,5,20,100), labels = c(0,80,95,99))) +
  geom_point(data = dat4845[which(dat4845$assay==0),], aes(x = conc, y = ratio), colour='darkgreen', shape = 1,alpha=.95) +
  geom_line(data = dat4845[which(dat4845$assay==0),], aes(x = conc, y = ratio, group = ID), alpha=.15) +
  #geom_line(data = dm1SMFA2m, aes(x=conc, y=value, group = variable), colour='blue', alpha=.2) +
  geom_line(data = mudf1x, aes(x=conc, y=mu), colour='darkgreen', alpha=.5) + themeJDC +
  geom_ribbon(data = PIdf1x, aes(x = conc, ymin = PI1, ymax = PI2), fill = "darkgreen", alpha = 0.25) +
  annotate("segment", x = 0.03, xend = 15.1, y = 5, yend = 5, linetype = 'dashed') +
  xlab("mAb conc. [\U03BCg/mL]") + ylab("Oocyst ratio(control/test)") + #ggtitle("Pfs4845 [SMFA]") + 
  annotate("text", x = 7.5, y=0.46, label = deparse(bquote("IC"[80]~"="~.(m1x)~"\U03BCg/mL ["~.(m1z)~","~.(m1y)~"]")),
           colour = "darkgreen", size = 5,parse = T)
pl3
#ggsave(filename = "whatever3.pdf",height = 5, width = 6,  device = cairo_pdf)

# 5000 samples 
ae_zeros <- matrix(0,5000,62)
pred2 <- list( #Fake data for DMFA? BF
  sqconc = cf.seq,
  logratio = rep(0,l1),#empty outcome
  ratio = rep(0,l1),#empty outcome
  assay = rep(1,l1), #
  country = rep(0,l1), #
  ID = factor(rep("test",l1)), # Is this ok?
  antibody = rep("4845",l1),
  conc = cf.seq*cf.seq
)

#ensemble4845_2 <- ensemble(m5me_4845, m6me_4845, m5_4845, m6_4845, m2me_4845, m1me_4845, m1_4845, m2_4845, data = pred2,
#                           replace = list(ae = ae_zeros))
ensemble4845_2 <- ensemble(m5me_4845, m6me_4845, m2me_4845, m1me_4845, data = pred2,
                           replace = list(ae = ae_zeros))
mu_2 <- apply(ensemble4845_2$link, 2, mean)
mu.PI_2 <- apply(ensemble4845_2$link, 2, PI, prob=0.95)

mudf2x <- data.frame(conc = cf.seq*cf.seq, mu = 10**mu_2)
PIdf2x <- data.frame(conc = cf.seq*cf.seq, PI1 = 10**mu.PI_2[1,], PI2 = 10**mu.PI_2[2,])

#Pick out
m2x <- 0
m2x2 <- 0
#sprintf('%.2f',a) # 2 digits after decimal
for(i in 1:(length(mu))){
  if(mudf2x$mu[i] > 5){
    print("Mean IC80:")
    print(mudf2x$conc[i])
    m2x <- round(mudf2x$conc[i],1)
    m2x2 <- sprintf('%.1f',mudf2x$conc[i])
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

# 5000 samples 
ae_zeros <- matrix(0,5000,62)
pred3 <- list( #Fake data for DMFA? CAM
  sqconc = cf.seq,
  logratio = rep(0,l1),#empty outcome
  ratio = rep(0,l1),#empty outcome
  assay = rep(1,l1), #
  country = rep(1,l1), #
  ID = factor(rep("test",l1)), # Is this ok?
  antibody = rep("25",l1),
  conc = cf.seq*cf.seq
)

#ensemble4845_3 <- ensemble(m5me_4845, m6me_4845, m5_4845, m6_4845, m2me_4845, m1me_4845, m1_4845, m2_4845, data = pred3,
#                           replace = list(ae = ae_zeros))
ensemble4845_3 <- ensemble(m5me_4845, m6me_4845, m2me_4845, m1me_4845, data = pred3,
                           replace = list(ae = ae_zeros))
mu_3 <- apply(ensemble4845_3$link, 2, mean)
mu.PI_3 <- apply(ensemble4845_3$link, 2, PI, prob=0.95)

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
pl4 <- ggplot() + geom_point(data = dat4845[which(dat4845$assay==1),], aes(x = conc, y = ratio, color=factor(country)), alpha=.5) +
  geom_line(data = dat4845[which(dat4845$assay==1),], aes(x = conc, y = ratio, group = ID), alpha=.125) +
  scale_x_continuous(trans='sqrt', limits = c(0.07,15.1), breaks = c(0.16,0.7,1.25,2.5,5,10,15)) + 
  themeJDC + scale_y_continuous(trans='log10',limits = c(0.25,330),
  sec.axis = sec_axis(~.*1, name = "TRA (%)", breaks = c(1,5,20,100), labels = c(0,80,95,99))) +
  guides(colour=guide_legend(title = "Country")) + scale_color_discrete(labels = c("BF","CAM")) +
  geom_line(data = mudf2x, aes(x=conc, y=mu), colour='red', alpha=.5) +
  geom_ribbon(data = PIdf2x, aes(x = conc, ymin = PI1, ymax = PI2), fill = "red", alpha = 0.25) +
  geom_line(data = mudf3x, aes(x=conc, y=mu), colour='turquoise', alpha=.5) +
  geom_ribbon(data = PIdf3x, aes(x = conc, ymin = PI1, ymax = PI2), fill = "turquoise", alpha = 0.25) +
  annotate("segment", x = 0.07, xend = 15.1, y = 5, yend = 5, linetype = 'dashed') +
  theme(legend.position = c(0.125,0.85)) +
  xlab("mAb conc. [\U03BCg/mL]") + ylab("Oocyst ratio(control/test)") + #ggtitle("Pfs4845 [DMFA]") +
  annotate("text", x = 8, y=0.3, label = deparse(bquote("IC"[80]~"="~.(m2x2)~"\U03BCg/mL ["~.(m2z)~","~.(m2y)~"]")),
           colour = "red", size = 5,parse = T) +
  annotate("text", x = 8, y=0.65, label = deparse(bquote("IC"[80]~"="~.(m3x)~"\U03BCg/mL ["~.(m3z)~","~.(m3y)~"]")),
           colour = "turquoise", size = 5,parse = T)
pl4
#ggsave(filename = "whatever4.pdf",height = 5, width = 6,  device = cairo_pdf)

