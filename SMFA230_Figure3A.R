library(ggplot2)
library(reshape2)
library(scales)
library(rstan)
library(rethinking)

themeJDC <- theme(panel.background = element_blank(),panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black" ,fill = NA), legend.background = element_blank(),
                  legend.key = element_blank(), legend.title = element_text(), axis.text = element_text(size = 11.5),
                  axis.title = element_text(size=12.2), legend.text = element_text(size = 10.9))

#Check directory
getwd()
datt <- read.csv("SMFA_230_v2.csv") 
str(datt)
table(datt$strain)
colnames(datt)[4] <- "expt"

#recale data as before?
datt$ratio <- 1/(1-datt$value/100)
datt$ratio[is.infinite(datt$ratio)==T] <- 200

#Desired colour scheme

# NF135: 66C0BE
# NF149: E42520
# NF175: 8A1812
# NF176: F5A216
# NF183: 28316F
# NF54: 248B37

# ggplot() + geom_point(data = datt, aes(y=ratio, x=conc_230,color = strain)) + 
#   themeJDC + xlab("mAb conc. [\U03BCg/mL]") + ylab("Oocyst ratio(control/test)") +
#   scale_x_continuous(trans='sqrt', breaks = c(1,5,10,16,32,64,100,128)) + 
#   scale_y_continuous(trans='log10', breaks = c(1,2,5,10,20,50,100,200)) +
#   scale_color_manual(name = "Strain",
#       values=c("#66C0BE","#E42520","#8A1812","#F5A216","#28316F","#248B37")) +
#   annotate("segment", x = 0.03, xend = 130, y = 5, yend = 5, linetype = 'dashed') 

######################################################################

#Rescale data for regression
datt$logratio <- log10(datt$ratio)
datt$sqconc <- sqrt(datt$conc_230)

#Carry out regression independenly for each strain
datt1 <- datt[datt$strain=="NF135",]
dim(datt1)
#datt1$exp <- droplevels(datt1$exp)
datt2 <- datt[datt$strain=="NF149",]
dim(datt2)
datt3 <- datt[datt$strain=="NF175",]
dim(datt3)
datt4 <- datt[datt$strain=="NF176",]
dim(datt4)
datt5 <- datt[datt$strain=="NF183",]
dim(datt5)
datt6 <- datt[datt$strain=="NF54",]
dim(datt6)

iterz <- 10000

###     M1:  NF135     ###
mod1 <- map2stan(
  alist(
    logratio ~ dnorm(mu , sigma),
    mu <-  a*sqconc,
    a ~ dexp(1), 
    sigma ~ dexp(1) # Prior on sigma
  ),
  data = datt1, iter=iterz , chains=3, cores = 3
)
precis(mod1,prob=0.9) # Model summary


###     M1:  NF149     ###
mod2 <- map2stan(
  alist(
    logratio ~ dnorm(mu , sigma),
    mu <-  a*sqconc,
    a ~ dexp(1), # Prior on slope
    sigma ~ dexp(1) # Prior on sigma
  ),
  data = datt2, iter=iterz , chains=3, cores = 3, control = list(adapt_delta = 0.94)
)
precis(mod2,prob=0.9) # Model summary

###     M3:  NF175     ###
mod3 <- map2stan(
  alist(
    logratio ~ dnorm(mu , sigma),
    mu <-  a*sqconc,
    a ~ dexp(1),#dexp(1), #
    sigma ~ dexp(1) # Prior on sigma
  ),
  data = datt3, iter=iterz , chains=3, cores = 3, control = list(adapt_delta = 0.94)
)
precis(mod3,prob=0.9) # Model summary

###     M4:  NF176     ###
mod4 <- map2stan(
  alist(
    logratio ~ dnorm(mu , sigma),
    mu <-  a*sqconc,
    a ~ dexp(1), # 
    sigma ~ dexp(1) # Prior on sigma
  ),
  data = datt4, iter=iterz , chains=3, cores = 3
)
precis(mod4,prob=0.9) # Model summary

###     M5:  NF183     ###
mod5 <- map2stan(
  alist(
    logratio ~ dnorm(mu , sigma),
    mu <-  a*sqconc,
    a ~ dexp(1), # 
    sigma ~ dexp(1) # Prior on sigma
  ),
  data = datt5, iter=iterz , chains=3, cores = 3
)
precis(mod5,prob=0.9) # Model summary

###     M6:  NF54     ###
mod6 <- map2stan(
  alist(
    logratio ~ dnorm(mu , sigma),
    mu <-  sqconc*(a + ae[expt]),
    a ~ dexp(1),
    ae[expt] ~ dnorm(0, sigma2),
    #alpha ~ dnorm(0,1), # Prior on slope
    sigma ~ dexp(1), # Prior on sigma
    sigma2 ~ dexp(1) # Prior on sigma2
  ),
  data = datt6, iter=iterz , chains=3, cores = 3
)
precis(mod6,prob=0.9, depth=2) # Model summary

#
####### 1.
sq.seq1 <- seq(1,5,0.01)
mu1 <- link( mod1, n=10000, data=data.frame(sqconc=sq.seq1) )

mu1.mean <- apply( mu1 , 2 , mean )
mu1.HPDI <- apply( mu1 , 2 , HPDI , prob=0.95 )

#Whack it back on original scale
mu1.meanX <- data.frame(conc = sq.seq1*sq.seq1, mu = 10**mu1.mean)
PI1X <- data.frame(conc = sq.seq1*sq.seq1, PI1 = 10**mu1.HPDI[1,], PI2 = 10**mu1.HPDI[2,])

#Pick out
m1x <- 0
for(i in 1:(length(mu1.mean))){
  if(mu1.meanX$mu[i] > 5){
    print("Mean IC80:")
    print(mu1.meanX$conc[i])
    m1x <- round(mu1.meanX$conc[i],1)
    break
  }
}

#Pick out
m1y <- 0
for(i in 1:(length(mu1.mean))){
  if(PI1X$PI1[i] > 5){
    print("CI1:")
    print(PI1X$conc[i])
    m1y <- round(PI1X$conc[i],1)
    break
  }
}

#Pick out
m1z <- 0
for(i in 1:(length(mu1.mean))){
  if(PI1X$PI2[i] > 5){
    print("CI2:")
    print(PI1X$conc[i])
    m1z <- round(PI1X$conc[i],1)
    break
  }
}

# ggplot() + geom_line(data=mu1.meanX,aes(x=conc, y=mu), color = "#66C0BE") +
#   geom_point(data = datt1, aes(x=conc_230, y=ratio), color = "#66C0BE") +
#   scale_x_continuous(trans='sqrt') + scale_y_continuous(trans='log10') +
#   geom_ribbon(data=PI1X, aes(x=conc, ymin = PI1, ymax = PI2), fill = '#66C0BE',
#               alpha = 0.2) + themeJDC

####### 2.
sq.seq2 <- seq(1.3,11.5,0.2)
#sq.seq2 <- seq(1.3,39.5,0.2)
mu2 <- link( mod2, n=10000, data=data.frame(sqconc=sq.seq2) )

mu2.mean <- apply( mu2 , 2 , mean )
mu2.HPDI <- apply( mu2 , 2 , HPDI , prob=0.95 )

#Whack it back on original scale
mu2.meanX <- data.frame(conc = sq.seq2*sq.seq2, mu = 10**mu2.mean)
PI2X <- data.frame(conc = sq.seq2*sq.seq2, PI1 = 10**mu2.HPDI[1,],
                   PI2 = 10**mu2.HPDI[2,])

#Pick out
m2x <- 0
for(i in 1:(length(mu2.mean))){
  if(mu2.meanX$mu[i] > 5){
    print("Mean IC80:")
    print(mu2.meanX$conc[i])
    m2x <- round(mu2.meanX$conc[i],1)
    break
  }
}

#Pick out
m2y <- 0
for(i in 1:(length(mu2.mean))){
  if(PI2X$PI1[i] > 5){
    print("CI1:")
    print(PI2X$conc[i])
    m2y <- round(PI2X$conc[i],1)
    break
  }
}

#Pick out
m2z <- 0
for(i in 1:(length(mu1.mean))){
  if(PI2X$PI2[i] > 5){
    print("CI2:")
    print(PI2X$conc[i])
    m2z <- round(PI2X$conc[i],1)
    break
  }
}

# ggplot() + geom_line(data=mu2.meanX,aes(x=conc, y=mu), color = "#E42520") +
#   geom_point(data = datt2, aes(x=conc_230, y=ratio), color = "#E42520") +
#   scale_x_continuous(trans='sqrt') + scale_y_continuous(trans='log10') +
#   geom_ribbon(data=PI2X, aes(x=conc, ymin = PI1, ymax = PI2), fill = '#E42520',
#               alpha = 0.2) + themeJDC

####### 3.
sq.seq3 <- seq(1.5,11.5,0.2)
#sq.seq3 <- seq(1.5,911.5,0.2)
mu3 <- link( mod3, n=10000, data=data.frame(sqconc=sq.seq3) )

mu3.mean <- apply( mu3 , 2 , mean )
mu3.HPDI <- apply( mu3 , 2 , HPDI , prob=0.95 )

#Whack it back on original scale
mu3.meanX <- data.frame(conc = sq.seq3*sq.seq3, mu = 10**mu3.mean)
PI3X <- data.frame(conc = sq.seq3*sq.seq3, PI1 = 10**mu3.HPDI[1,],
                   PI2 = 10**mu3.HPDI[2,])

#Pick out
m3x <- 0
for(i in 1:(length(mu3.mean))){
  if(mu3.meanX$mu[i] > 5){
    print("Mean IC80:")
    print(mu3.meanX$conc[i])
    m3x <- round(mu3.meanX$conc[i],1)
    break
  }
}

#Pick out
m3y <- 0
for(i in 1:(length(mu3.mean))){
  if(PI3X$PI1[i] > 5){
    print("CI1:")
    print(PI3X$conc[i])
    m3y <- round(PI3X$conc[i],1)
    break
  }
}

#Pick out
m3z <- 0
for(i in 1:(length(mu3.mean))){
  if(PI3X$PI2[i] > 5){
    print("CI2:")
    print(PI3X$conc[i])
    m3z <- round(PI3X$conc[i],1)
    break
  }
}

# ggplot() + geom_line(data=mu3.meanX,aes(x=conc, y=mu), color = "#8A1812") +
#   geom_point(data = datt3, aes(x=conc_230, y=ratio), color = "#8A1812") +
#   scale_x_continuous(trans='sqrt') + scale_y_continuous(trans='log10') +
#   geom_ribbon(data=PI3X, aes(x=conc, ymin = PI1, ymax = PI2), fill = '#8A1812',
#               alpha = 0.2) + themeJDC

####### 4.
sq.seq4 <- seq(4,11.5,0.1)
#sq.seq4 <- seq(4,911.5,0.1)
mu4 <- link( mod4, n=10000, data=data.frame(sqconc=sq.seq4) )

mu4.mean <- apply( mu4 , 2 , mean )
mu4.HPDI <- apply( mu4 , 2 , HPDI , prob=0.95 )

#Whack it back on original scale
mu4.meanX <- data.frame(conc = sq.seq4*sq.seq4, mu = 10**mu4.mean)
PI4X <- data.frame(conc = sq.seq4*sq.seq4, PI1 = 10**mu4.HPDI[1,],
                   PI2 = 10**mu4.HPDI[2,])


#Pick out
m4x <- 0
for(i in 1:(length(mu4.mean))){
  if(mu4.meanX$mu[i] > 5){
    print("Mean IC80:")
    print(mu4.meanX$conc[i])
    m4x <- round(mu4.meanX$conc[i],1)
    break
  }
}

#Pick out
m4y <- 0
for(i in 1:(length(mu4.mean))){
  if(PI4X$PI1[i] > 5){
    print("CI1:")
    print(PI4X$conc[i])
    m4y <- round(PI4X$conc[i],1)
    break
  }
}

#Pick out
m4z <- 0
for(i in 1:(length(mu4.mean))){
  if(PI4X$PI2[i] > 5){
    print("CI2:")
    print(PI4X$conc[i])
    m4z <- round(PI4X$conc[i],1)
    break
  }
}


# ggplot() + geom_line(data=mu4.meanX,aes(x=conc, y=mu), color = "#F5A216") +
#   geom_point(data = datt4, aes(x=conc_230, y=ratio), color = "#F5A216") +
#   scale_x_continuous(trans='sqrt') + scale_y_continuous(trans='log10') +
#   geom_ribbon(data=PI4X, aes(x=conc, ymin = PI1, ymax = PI2), fill = '#F5A216',
#               alpha = 0.2) + themeJDC

####### 5.
sq.seq5 <- seq(1,4,0.01)
mu5 <- link( mod5, n=10000, data=data.frame(sqconc=sq.seq5) )

mu5.mean <- apply( mu5 , 2 , mean )
mu5.HPDI <- apply( mu5 , 2 , HPDI , prob=0.95 )

#Whack it back on original scale
mu5.meanX <- data.frame(conc = sq.seq5*sq.seq5, mu = 10**mu5.mean)
PI5X <- data.frame(conc = sq.seq5*sq.seq5, PI1 = 10**mu5.HPDI[1,],
                   PI2 = 10**mu5.HPDI[2,])


#Pick out
m5x <- 0
for(i in 1:(length(mu5.mean))){
  if(mu5.meanX$mu[i] > 5){
    print("Mean IC80:")
    print(mu5.meanX$conc[i])
    m5x <- round(mu5.meanX$conc[i],1)
    break
  }
}

#Pick out
m5y <- 0
for(i in 1:(length(mu5.mean))){
  if(PI5X$PI1[i] > 5){
    print("CI1:")
    print(PI5X$conc[i])
    m5y <- round(PI5X$conc[i],1)
    break
  }
}

#Pick out
m5z <- 0
for(i in 1:(length(mu5.mean))){
  if(PI5X$PI2[i] > 5){
    print("CI2:")
    print(PI5X$conc[i])
    m5z <- round(PI5X$conc[i],1)
    break
  }
}

# ggplot() + geom_line(data=mu5.meanX,aes(x=conc, y=mu), color = "#28316F") +
#   geom_point(data = datt5, aes(x=conc_230, y=ratio), color = "#28316F") +
#   scale_x_continuous(trans='sqrt') + scale_y_continuous(trans='log10') +
#   geom_ribbon(data=PI5X, aes(x=conc, ymin = PI1, ymax = PI2), fill = '#28316F',
#               alpha = 0.2) + themeJDC

####### 6.
sq.seq6 <- seq(0.2,4,0.01)
# 1000 samples by 4 actors
ae_zeros <- matrix(0,10000,9)
l6 <- length(sq.seq6)
mu6 <- link( mod6, n=10000, data=data.frame(sqconc=sq.seq6,
          expt = factor(rep("test",l6))), replace=list(ae=ae_zeros) )
#mu6 <- link( mod6 , data=data.frame(sqconc=sq.seq6) )

mu6.mean <- apply( mu6 , 2 , mean )
mu6.HPDI <- apply( mu6 , 2 , HPDI , prob=0.89 )
mu6.HPDI <- apply( mu6 , 2 , PI , prob=0.89 )

#Whack it back on original scale
mu6.meanX <- data.frame(conc = sq.seq6*sq.seq6, mu = 10**mu6.mean)
PI6X <- data.frame(conc = sq.seq6*sq.seq6, PI1 = 10**mu6.HPDI[1,],
                   PI2 = 10**mu6.HPDI[2,])


#Pick out
m6x <- 0
for(i in 1:(length(mu6.mean))){
  if(mu6.meanX$mu[i] > 5){
    print("Mean IC80:")
    print(mu6.meanX$conc[i])
    m6x <- round(mu6.meanX$conc[i],1)
    break
  }
}

#Pick out
m6y <- 0
for(i in 1:(length(mu6.mean))){
  if(PI6X$PI1[i] > 5){
    print("CI1:")
    print(PI6X$conc[i])
    m6y <- round(PI6X$conc[i],1)
    break
  }
}

#Pick out
m6z <- 0
for(i in 1:(length(mu6.mean))){
  if(PI6X$PI2[i] > 5){
    print("CI2:")
    print(PI6X$conc[i])
    m6z <- round(PI6X$conc[i],1)
    break
  }
}


# ggplot() + geom_line(data=mu6.meanX,aes(x=conc, y=mu), color = "#248B37") +
#   geom_point(data = datt6, aes(x=conc_230, y=ratio), color = "#248B37") +
#   scale_x_continuous(trans='sqrt', limits = c(0.1,16.1)) + 
#   scale_y_continuous(trans='log10',limits = c(0.21,320),
#     sec.axis = sec_axis(~.*1, name = "TRA (%)",
#                         breaks = c(1,5,20,100), labels = c(0,80,95,99))) +
#   geom_ribbon(data=PI6X, aes(x=conc, ymin = PI1, ymax = PI2), fill = '#248B37',
#               alpha = 0.2) + themeJDC

### Combine into one plot

#Add legend
# NF135: 66C0BE
# NF149: E42520
# NF175: 8A1812
# NF176: F5A216
# NF183: 28316F
# NF54: 248B37
cols <- c("NF135"="#66C0BE","NF149"="#E42520","NF175"="#8A1812",
          "NF176"="#F5A216","NF183"="#28316F","NF54"="#248B37")
cols2 <- c("NF135"="#66C0BE40","NF149"="#E4252040","NF175"="#8A181240",
          "NF176"="#F5A21640","NF183"="#28316F40","NF54"="#248B3740")

ggplot() +
  geom_line(data=mu1.meanX,aes(x=conc, y=mu, color = "NF135")) +
  geom_point(data = datt1, aes(x=conc_230, y=ratio, color = "NF135"), alpha=0.65) +
  geom_ribbon(data=PI1X, aes(x=conc, ymin = PI1, ymax = PI2, fill = "NF135"),
              alpha = 0.2) +
  geom_line(data=mu2.meanX,aes(x=conc, y=mu, color = "NF149")) +
  geom_point(data = datt2, aes(x=conc_230, y=ratio, color = "NF149"), alpha=0.65) +
  geom_ribbon(data=PI2X, aes(x=conc, ymin = PI1, ymax = PI2, fill = "NF149"),
              alpha = 0.2) +
  geom_line(data=mu3.meanX,aes(x=conc, y=mu, color = "NF175")) +
  geom_point(data = datt3, aes(x=conc_230, y=ratio, color = "NF175"), alpha=0.65) +
  geom_ribbon(data=PI3X, aes(x=conc, ymin = PI1, ymax = PI2, fill = "NF175"),
              alpha = 0.2) +
  geom_line(data=mu4.meanX,aes(x=conc, y=mu, color = "NF176")) +
  geom_point(data = datt4, aes(x=conc_230, y=ratio, color = "NF176"), alpha=0.65) +
  geom_ribbon(data=PI4X, aes(x=conc, ymin = PI1, ymax = PI2, fill = "NF176"),
              alpha = 0.2) +
  geom_line(data=mu5.meanX,aes(x=conc, y=mu, color = "NF183")) +
  geom_point(data = datt5, aes(x=conc_230, y=ratio, color = "NF183"), alpha=0.65) +
  geom_ribbon(data=PI5X, aes(x=conc, ymin = PI1, ymax = PI2, fill = "NF183"),
              alpha = 0.2) +
  xlab("mAb conc. [\U03BCg/mL]") +
  ylab("Oocyst ratio(control/test)") +
  annotate("segment", x = 0.09, xend = 130, y = 5, yend = 5, linetype = 'dashed') +
  geom_line(data=mu6.meanX,aes(x=conc, y=mu, color = "NF54")) +
  geom_point(data = datt6, aes(x=conc_230, y=ratio, color = "NF54"), alpha=0.65) +
  scale_x_continuous(trans='sqrt', breaks = c(1,2.5,5,10,16,32,64,100,128)) +
  scale_y_continuous(trans='log10', breaks = c(1,2,5,10,20,50,100,200),
    sec.axis = sec_axis(~.*1, name = "TRA (%)", breaks = c(1,5,20,100),
                        labels = c(0,80,95,99))) +
  geom_ribbon(data=PI6X, aes(x=conc, ymin = PI1, ymax = PI2, fill = "NF54"),
              alpha = 0.2 ) + themeJDC +
  scale_colour_manual(name="Strain",values=cols) +
  scale_fill_manual(name="Strain",values=cols) + 
  guides(colour = guide_legend(override.aes = list(alpha = .23)))

#Add the IC80s
#Some need to be done manually 
m2x <- 170
m2y <- 980
m2z <- 69
m3x <- 3900
m3z <- 1200
m3y <- "NA"
m4x <- 530
m4z <- 180
m4y <- "NA"
new_plot_v2 <- ggplot() +
  geom_line(data=mu1.meanX,aes(x=conc, y=mu, color = "NF135")) +
  geom_point(data = datt1, aes(x=conc_230, y=ratio, color = "NF135"), alpha=0.65) +
  geom_ribbon(data=PI1X, aes(x=conc, ymin = PI1, ymax = PI2, fill = "NF135"),
              alpha = 0.2) +
  geom_line(data=mu2.meanX,aes(x=conc, y=mu, color = "NF149")) +
  geom_point(data = datt2, aes(x=conc_230, y=ratio, color = "NF149"), alpha=0.65) +
  geom_ribbon(data=PI2X, aes(x=conc, ymin = PI1, ymax = PI2, fill = "NF149"),
              alpha = 0.2) +
  geom_line(data=mu3.meanX,aes(x=conc, y=mu, color = "NF175")) +
  geom_point(data = datt3, aes(x=conc_230, y=ratio, color = "NF175"), alpha=0.65) +
  geom_ribbon(data=PI3X, aes(x=conc, ymin = PI1, ymax = PI2, fill = "NF175"),
              alpha = 0.2) +
  geom_line(data=mu4.meanX,aes(x=conc, y=mu, color = "NF176")) +
  geom_point(data = datt4, aes(x=conc_230, y=ratio, color = "NF176"), alpha=0.65) +
  geom_ribbon(data=PI4X, aes(x=conc, ymin = PI1, ymax = PI2, fill = "NF176"),
              alpha = 0.2) +
  geom_line(data=mu5.meanX,aes(x=conc, y=mu, color = "NF183")) +
  geom_point(data = datt5, aes(x=conc_230, y=ratio, color = "NF183"), alpha=0.65) +
  geom_ribbon(data=PI5X, aes(x=conc, ymin = PI1, ymax = PI2, fill = "NF183"),
              alpha = 0.2) + theme(legend.position = "right") +
  xlab("mAb conc. [\U03BCg/mL]") +
  ylab("Oocyst ratio(control/test)") +
  annotate("segment", x = 0.09, xend = 130, y = 5, yend = 5, linetype = 'dashed') +
  geom_line(data=mu6.meanX,aes(x=conc, y=mu, color = "NF54")) +
  geom_point(data = datt6, aes(x=conc_230, y=ratio, color = "NF54"), alpha=0.65) +
  scale_x_continuous(trans='sqrt', breaks = c(1,2.5,5,10,16,32,64,100,128)) +
  scale_y_continuous(trans='log10', breaks = c(1,2,5,10,20,50,100,200),
                     sec.axis = sec_axis(~.*1, name = "TRA (%)", breaks = c(1,5,20,100),
                                         labels = c(0,80,95,99))) +
  geom_ribbon(data=PI6X, aes(x=conc, ymin = PI1, ymax = PI2, fill = "NF54"),
              alpha = 0.2 ) + themeJDC +
  scale_colour_manual(name="Strain",values=cols) +
  scale_fill_manual(name="Strain",values=cols) +
   annotate("text", x = 70, y=198, 
     label = deparse(bquote("IC"[80]~"="~.(m1x)~"\U03BCg/mL  ["~.(m1z)~","~.(m1y)~"]")),
          colour = "#66C0BE", size = 4.0,parse = T) +
   annotate("text", x = 71, y=145, 
             label = deparse(bquote("IC"[80]~"="~.(m2x)~"\U03BCg/mL ["~.(m2z)~","~.(m2y)~"]")),
             colour = "#E42520", size = 4.0,parse = T) +
    annotate("text", x = 74.5, y=109, 
             label = deparse(bquote("IC"[80]~"="~.(m3x)~"\U03BCg/mL ["~.(m3z)~","~.(m3y)~"]")),
             colour = "#8A1812", size = 4.0,parse = T) +
    annotate("text", x = 71.5, y=79, 
             label = deparse(bquote("IC"[80]~"="~.(m4x)~"\U03BCg/mL ["~.(m4z)~","~.(m4y)~"]")),
             colour = "#F5A216", size = 4.0,parse = T) +
   annotate("text", x = 70, y=56, 
            label = deparse(bquote("IC"[80]~"="~.(m5x)~"\U03BCg/mL ["~.(m5z)~","~.(m5y)~"]")),
            colour = "#28316F", size = 4.0,parse = T) +
   annotate("text", x = 70, y=39, 
            label = deparse(bquote("IC"[80]~"="~.(m6x)~"\U03BCg/mL ["~.(m6z)~","~.(m6y)~"]")),
            colour = "#248B37", size = 4.0,parse = T) 
new_plot_v2
#Note: correctly outputting the greek symbols can be an issue on some computers
ggsave(file="Figure3A.pdf",height = 6, width = 8) 

