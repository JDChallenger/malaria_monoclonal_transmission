

themeJDC <- theme(panel.background = element_blank(),panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black" ,fill = NA), legend.background = element_blank(),
                  legend.key = element_blank(), legend.title = element_text(), axis.text = element_text(size = 11.8),
                  axis.title = element_text(size=12.2), legend.text = element_text(size = 10.9))

#Check directory is OK
getwd()

dat1 <- read.csv("DMFA_230_all.csv") # Updated, without truncating increased transmission
str(dat1)

#Delete column of comments
dat1a <- dat1[,1:3]
dat1a$ratio <- 1/(1-dat1a$thing/100)
dat1a$ratio[is.infinite(dat1a$ratio)==T] <- 200
dat1b <- dat1a[,c(3,1,4)]
dat1b$assay <- 1 #DMFA
dat1b$antibody <- "230"
dat1b$country <- 0 #"BF" #Now question is: is this Cameroon (0/1)

names <- c("conc","ID","ratio","assay","antibody","country") # Note the order!
colnames(dat1b) <- names
dat1b

dat2 <- read.csv("P230_NF54_v2.csv") 
str(dat2)

#Values are expressed as a percentage of the control
j2 <- dim(dat2)[2]
i2 <- dim(dat2)[1]
dat2b <- dat2
for(i in 1:i2){
  for(j in 2:j2){
    print(dat2[i,j])
    if(is.na(dat2[i,j]) == F){

    if (dat2[i,j] < 0.5){
      dat2b[i,j] <- 200
      }  else{
          dat2b[i,j] <- 100/dat2[i,j]
      }
      #print(j)
    }
    #else{
    #  print("NA")
    #}
    }
}

dat2bm <- melt(dat2b, id.vars = "conc_abs") #Standardise names?
#drop NAs, as they interfere with plotting
dat2c <- dat2bm[complete.cases(dat2bm),]
# ggplot(dat2c, aes(x=conc_abs, y=value, colour=variable)) + geom_line() + geom_point() + themeJDC +
#   scale_x_continuous(trans='sqrt') + scale_y_continuous(trans='log10') + ggtitle("Pfs230 (NF54)")

dat2X <-dat2c
dat2X$assay <- 0 #SMFA
dat2X$antibody <- "230"
dat2X$country <- 0#"NA"  #Now question is: is this Cameroon (0/1)

names <- c("conc","ID","ratio","assay","antibody","country") # Note the order!
colnames(dat2X) <- names
dat2X

dat3 <- read.csv("P4845_BF.csv") 
str(dat3)
#Values are expressed as a percentage of the control
j3 <- dim(dat3)[2]
i3 <- dim(dat3)[1]
dat3b <- dat3
for(i in 1:i3){
  for(j in 2:j3){
    print(dat3[i,j])
    if(is.na(dat3[i,j]) == F){
      
      if (dat3[i,j] < 0.5){
        dat3b[i,j] <- 200
      }  else{
        dat3b[i,j] <- 100/dat3[i,j]
      }
      #print(j)
    }
    #else{
    #  print("NA")
    #}
  }
}
#View(dat3b)
dat3bm <- melt(dat3b, id.vars = "conc_4845") #Standardise names?
#drop NAs, as they interfere with plotting
dat3c <- dat3bm[complete.cases(dat3bm),]
# ggplot(dat3c, aes(x=conc_4845, y=value, colour=variable)) + geom_line() + geom_point() + themeJDC +
#   scale_x_continuous(trans='sqrt') + scale_y_continuous(trans='log10') + ggtitle("Pfs45/48 (BF)")

dat3X <-dat3c
dat3X$assay <- 1 #DMFA
dat3X$antibody <- "4845"
dat3X$country <- 0 #"BF" #Now question is: is this Cameroon (0/1)

names <- c("conc","ID","ratio","assay","antibody","country") # Note the order!
colnames(dat3X) <- names
dat3X


dat4 <- read.csv("P4845_CAM.csv") 
str(dat4)
#View(dat4) # This is panel E 
#Values are expressed as a percentage of the control
dat4b <- dat4[1:3,] #What does zero antibodies mean?
j4 <- dim(dat4b)[2]
i4 <- dim(dat4b)[1]
for(i in 1:i4){
  for(j in 2:j4){
    #dat4b[i,j] <- dat4[4,j]/dat4b[i,j] WRONG??
    print(dat4[i,j])
    if(is.na(dat4[i,j]) == F){
      
      if (dat4[i,j] < 0.5){
        dat4b[i,j] <- 200
      }  else{
        dat4b[i,j] <- 100/dat4[i,j]
      }
      #print(j)
    }
    #else{
    #  print("NA")
    #}
  }
}
#dat4b[mapply(is.infinite,dat4b)] <- 200
#View(dat4b)
dat4bm <- melt(dat4b, id.vars = "conc_4845") #Standardise names?
# ggplot(dat4bm, aes(x=conc_4845, y=value, colour=variable)) + geom_line() + geom_point() + themeJDC +
#   scale_x_continuous(trans='sqrt') + scale_y_continuous(trans='log10') + ggtitle("Pfs48/45 (CAM)")
# 
# ggplot() + geom_line(data=dat4bm, aes(x=conc_4845, y=value,group=variable), colour='blue') +  themeJDC +
#   geom_line(data=dat3c, aes(x=conc_4845, y=value, group=variable), colour='red') +
#   scale_x_continuous(trans='sqrt') + scale_y_continuous(trans='log10') + ggtitle("Pfs48/45")

dat4X <-dat4bm
dat4X$assay <- 1 #DMFA
dat4X$antibody <- "4845"
dat4X$country <- 1 #"CAM" #Now question is: is this Cameroon (0/1)

names <- c("conc","ID","ratio","assay","antibody","country") # Note the order!
colnames(dat4X) <- names
dat4X


dat5 <- read.csv("P4845_NF54.csv") 
str(dat5)
#View(dat5) # B? variable 'log10' might be the antibody conc.
#new column
dat5$conc <- 10**(dat5$log10)
#Values are expressed as a percentage of the control
j5 <- dim(dat5)[2]-1 #skipping new variable
i5 <- dim(dat5)[1]
dat5b <- dat5
for(i in 1:i5){
  for(j in 2:j5){
    print(dat5[i,j])
    if(is.na(dat5[i,j]) == F){
      
      if (dat5[i,j] < 0.5){
        dat5b[i,j] <- 200
      }  else{
        dat5b[i,j] <- 100/dat5[i,j]
      }
      #print(j)
    }
    #else{
    #  print("NA")
    #}
  }
}
#drop log10 #Be careful with order of variables
dat5b<- dat5b[ -1]
#View(dat5b)
dat5bm <- melt(dat5b, id.vars = "conc") #Standardise names?
#drop NAs, as they interfere with plotting
dat5c <- dat5bm[complete.cases(dat5bm),]
# ggplot(dat5c, aes(x=conc, y=value, colour=variable)) + geom_line() + geom_point() + themeJDC +
#   scale_x_continuous(trans='sqrt') + scale_y_continuous(trans='log10') + ggtitle("Pfs48/45 (NF54)")

dat5X <-dat5c
dat5X$assay <- 0 #DMFA
dat5X$antibody <- "4845"
dat5X$country <- 0 #"NA" #Now question is: is this Cameroon (0/1)

names <- c("conc","ID","ratio","assay","antibody","country") # Note the order!
colnames(dat5X) <- names
dat5X

dat6 <- read.csv("Pfs25_BF.csv") 
str(dat6)
#View(dat6) # D?
#Values are expressed as a percentage of the control
dat6b <- dat6[1:3,]
j6 <- dim(dat6b)[2]
i6 <- dim(dat6b)[1]
for(i in 1:i6){
  for(j in 2:j6){
   # dat6b[i,j] <- dat6[4,j]/dat6b[i,j] WRONG?
    print(dat6[i,j])
    if(is.na(dat6[i,j]) == F){
      
      if (dat6[i,j] < 0.5){
        dat6b[i,j] <- 200
      }  else{
        dat6b[i,j] <- 100/dat6[i,j]
      }
      #print(j)
    }
    #else{
    #  print("NA")
    #}
  }
}
#dat6b[mapply(is.infinite,dat6b)] <- 200
#View(dat6b)
dat6bm <- melt(dat6b, id.vars = "conc_P25") #Standardise names?
#ggplot(dat6bm, aes(x=conc_P25, y=value, colour=variable)) + geom_line() + geom_point() + themeJDC +
#  scale_x_continuous(trans='sqrt') + scale_y_continuous(trans='log10') + ggtitle("Pfs25 (BF)")

dat6X <-dat6bm
dat6X$assay <- 1 #DMFA
dat6X$antibody <- "25"
dat6X$country <- 0 #"BF" #Now question is: is this Cameroon (0/1)

names <- c("conc","ID","ratio","assay","antibody","country") # Note the order!
colnames(dat6X) <- names
dat6X


dat7 <- read.csv("Pfs25_CAM.csv") 
str(dat7)
#View(dat7) # D?
#These seem to be raw counts. Find ratio with control
dat7b <- dat7[1:3,]
j7 <- dim(dat7b)[2]
i7 <- dim(dat7b)[1]
for(i in 1:i7){
  for(j in 2:j7){
   # dat7b[i,j] <- dat7[4,j]/dat7b[i,j] WRONG?
    print(dat7[i,j])
    if(is.na(dat7[i,j]) == F){
      
      if (dat7[i,j] < 0.5){
        dat7b[i,j] <- 200
      }  else{
        dat7b[i,j] <- 100/dat7[i,j]
      }
      #print(j)
    }
  }
}
#dat7b[mapply(is.infinite,dat7b)] <- 200
#View(dat7b)
dat7bm <- melt(dat7b, id.vars = "conc_P25") #Standardise names?
# ggplot(dat7bm, aes(x=conc_P25, y=value, colour=variable)) + geom_line() + geom_point() + themeJDC +
#   scale_x_continuous(trans='sqrt') + scale_y_continuous(trans='log10') + ggtitle("Pfs25 (CAM)")
# 
# ggplot() + geom_line(data=dat7bm, aes(x=conc_P25, y=value,group=variable), colour='blue') +  themeJDC +
#   geom_line(data=dat6bm, aes(x=conc_P25, y=value, group=variable), colour='red') +
#   scale_x_continuous(trans='sqrt') + scale_y_continuous(trans='log10') + ggtitle("Pfs25")

dat7X <-dat7bm
dat7X$assay <- 1 #DMFA
dat7X$antibody <- "25"
dat7X$country <- 1#"CAM" #Now question is: is this Cameroon (0/1)

names <- c("conc","ID","ratio","assay","antibody","country") # Note the order!
colnames(dat7X) <- names
dat7X


dat8 <- read.csv("Pfs25_NF54.csv") 
str(dat8)
#View(dat8) #A
#Values are expressed as a percentage of the control
dat8b <- dat8
j8 <- dim(dat8b)[2]
i8 <- dim(dat8b)[1]
for(i in 1:i8){
  for(j in 2:j8){
    print(dat8[i,j])
    if(is.na(dat8[i,j]) == F){
      
      if (dat8[i,j] < 0.5){
        dat8b[i,j] <- 200
      }  else{
        dat8b[i,j] <- 100/dat8[i,j]
      }
      #print(j)
    }
  }
}
#View(dat8b)

dat8bm <- melt(dat8b, id.vars = "conc")
# ggplot(dat8bm, aes(x=conc, y=value, colour=variable)) + geom_line() + geom_point() + themeJDC +
#   scale_x_continuous(trans='sqrt') + scale_y_continuous(trans='log10') + ggtitle("Pfs25 (NF54)")

dat8X <-dat8bm
dat8X$assay <- 0 #SMFA
dat8X$antibody <- "25"
dat8X$country <- 0 #"NA"  #Now question is: is this Cameroon (0/1)

names <- c("conc","ID","ratio","assay","antibody","country") # Note the order!
colnames(dat8X) <- names
dat8X

datY <- rbind(dat1b, dat2X, dat3X, dat4X, dat5X, dat6X, dat7X, dat8X)


#transform variables now?
datY$assay <- as.integer(datY$assay)
datY$logratio <- log10(datY$ratio)
datY$sqconc <- sqrt(datY$conc)
namesY <- c("conc","ID","ratio","assay","antibody","country","logratio","sqconc") # Note the order!
names(datY) <- namesY

#Subset by antibody
dat25 <- subset(datY, antibody == "25")

rownames(dat25) <- 1:nrow(dat25)
dat4845 <- subset(datY, antibody == "4845")

rownames(dat4845) <- 1:nrow(dat4845)
dat230 <- subset(datY, antibody == "230")

rownames(dat230) <- 1:nrow(dat230)

#Function to get IC80 conc from slope
ic80 <- function(slope){
  (log10(5)/slope)*(log10(5)/slope)
  #return(result)
}
#check
ic80(0.51)
