library(rstan)
library(rethinking) #See README for installation instructions
library(ggplot2)
library(reshape2)

#in 'ViewData.R' I created a single dataframe 'datY' containing all 8 data sets, in the same format
#Load the data
source('ViewData.R')
#View(dat230)

#Data definitions: 'conc' = antibody concentration, 'ID' = experiment or donor ID, 
#'ratio' = ratio of oocyst counts between control group#& intervention group
#'logratio' = log10 (ratio). 'assay' = 1 if assay = DMFA. 'country' = 1 if country = Cameroon
#'sqconc' = square root of concentration


#Just use SMFA data for regression modelling
dat230Y <- dat230[dat230$assay==0,]

####################################################################################################  
#############################              Pfs230                 ##################################
####################################################################################################

iterz <- 10000

###     M1me: Single slope (with mixed effects)      ###
m1me_230 <- map2stan(
  alist(
    logratio ~ dnorm(mu , sigma),
    mu <-  sqconc*(a + ae[ID]),
    ae[ID] ~ dnorm(0, sigma2),
    a ~ dnorm(0,1), # Prior on slope
    sigma ~ dexp(1), # Prior on sigma
    sigma2 ~ dexp(1) # Prior on sigma2
  ),
  data = dat230Y, iter=iterz , chains=3, cores = 3
)
precis(m1me_230,prob=0.9,depth=2) # Model summary

# Note: here only one country (BF) is present, so no country-specific effect

####################################################################################################  
##############     Pfs230: Make Plots and generate interval for $IC_{80}$       ####################
####################################################################################################

#Just use model m1me this time- only one mixed effects model

cf.seq <- seq(0.1,5,0.01)
ae_zeros <- matrix(0,5000,9)
l1 <- length(cf.seq)
pred1 <- list( #Fake data for SMFA?
  sqconc = cf.seq,
  logratio = rep(0,l1),#empty outcome
  ratio = rep(0,l1),#empty outcome
  assay = rep(0,l1), #
  country = rep(0,l1), #
  ID = factor(rep("test",l1)), # Is this ok?
  antibody = rep("230",l1),
  conc = cf.seq*cf.seq
)
ensemble230 <- ensemble(m1me_230, data = pred1, replace = list(ae = ae_zeros))# I've taken m1me out of this- b has a big effect. Don't fully understand this

mu <- apply(ensemble230$link, 2, mean)
mu.PI <- apply(ensemble230$link, 2, PI, prob=0.92)

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

pl5 <- ggplot() + geom_point(data = dat230[which(dat230$assay==0),], aes(x = conc, y = ratio),colour = 'darkgreen',alpha=.7) + 
  geom_line(data = dat230[which(dat230$assay==0),], aes(x = conc, y = ratio, group = ID), alpha=.15) +
  scale_x_continuous(trans='sqrt', limits = c(0.1,16.1), breaks = c(0.15,0.6,2.5,5,10,16)) + 
  scale_y_continuous(trans='log10',limits = c(0.21,380),
  sec.axis = sec_axis(~.*1, name = "TRA (%)", breaks = c(1,5,20,100), labels = c(0,80,95,99))) +
  themeJDC + geom_line(data = mudf1x, aes(x=conc, y=mu), colour='darkgreen', alpha=.5) +
  geom_ribbon(data = PIdf1x, aes(x = conc, ymin = PI1, ymax = PI2), fill = "darkgreen", alpha = 0.25) +
  annotate("segment", x = 0.1, xend = 15, y = 5, yend = 5, linetype = 'dashed') +
  xlab("mAb conc. [\U03BCg/mL]") +
  ylab("Oocyst ratio(control/test)") + #ggtitle("Pfs230 [SMFA]") + 
  annotate("text", x = 8.5, y=0.42, label = deparse(bquote("IC"[80]~"="~.(m1x)~"\U03BCg/mL  ["~.(m1z)~","~.(m1y)~"]")),
           colour = "darkgreen", size = 5,parse = T) 
pl5
#ggsave(filename = "whatever5.pdf",height = 5, width = 6,  device = cairo_pdf)

#DMFA only
pl6 <- ggplot() + geom_point(data = dat230[which(dat230$assay==1),], aes(x = conc, y = ratio, colour=factor(country))) + themeJDC +
  geom_line(data = dat230[which(dat230$assay==1),], aes(x = conc, y = ratio, group=ID), alpha = .4) +
  scale_x_continuous(trans='sqrt', limits = c(2.1,15.1), breaks = c(2.5,5,10,15)) + 
  scale_y_continuous(trans='log10',limits = c(0.21,230), 
  sec.axis = sec_axis(~.*1, name = "TRA (%)", breaks = c(1,5,20,100),
                      labels = c(0,80,95,99))) +
  guides(colour=guide_legend(title = "Country")) + scale_color_discrete(labels = c("BF")) +
  annotate("segment", x = 2.1, xend = 15, y = 5, yend = 5, linetype = 'dashed') +
  theme(legend.position = c(0.1,0.85)) +
  xlab("mAb conc. [\U03BCg/mL]") + 
  ylab("Oocyst ratio(control/test)") 
pl6
#ggsave(filename = "whatever6.pdf",height = 5, width = 6,  device = cairo_pdf)