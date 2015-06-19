# This function attempts to revise and modularize the MIMIC code to allow for smoother alternative analysis

library(rootSolve)
source('lib/MIMIC2_ODE.R')
source('lib/MIMIC2_Parameters.R')
#---------------------------------------------------------
# (A)       Read in site level data
#---------------------------------------------------------
data <- read.csv("data/LTER_SITE_1.csv") #site level forcing variables
data$ANPP    <- data$ANPP / 2   # convert to gC/m2/y from g/m2/y

siteIndex <- 1

LITtype  <- c('TRAEf', 'PIREf','THPLf','ACSAf','QUPRf','DRGLf')
bagMET   <- c(10.6, 36.2, 37.4, 56.8, 37.1, 49.3) #from Gordon's LitterCharacteristics.txt
bagLIG   <- c(16.2, 19.2, 26.7, 15.9, 23.5, 10.9) # % from Gordon's LitterCharacteristics.txt
bagN     <- c(0.38, 0.59, 0.62, 0.81, 1.03, 1.97) # %N 
bagCN    <- c(133.3,92.7, 83.1, 61.8, 50.5, 24.2)
calcN    <- (1 / bagCN) / 2.5 * 100    
calcMET  <- 0.85 - 0.013 * bagLIG/calcN       #as calculated in DAYCENT
bagMET   <- bagMET / 100
bagMET   <- calcMET
fMET     <- mean(calcMET) 
##################
##I Litter inputs
##################
I        <- array(NA, dim=2)              #Litter inputs to MET/STR
depth    <- 30
EST_LIT  <- data$ANPP[siteIndex] / (365*24) * 1e3 / 1e4  #gC/m2/h (from g/m2/y, Knapp et al. Science 2001)
I[1]     <- (EST_LIT / depth) * fMET      #partitioned to layers, fMET is site specific
I[2]     <- (EST_LIT / depth) * (1-fMET)

par0 <- makeInitalPar(I = I,
                      fMET=fMET, TSOI=data$MAT[siteIndex], 
                      fCLAY=data$CLAY2[siteIndex]/100, #convert from clay fraction to %
                      Tao_MOD1=sqrt(data$ANPP[siteIndex]/100))#basicaily standardize against NWT

y0    <- c( LIT_1 = I[1], LIT_2 = I[2], 
            MIC_1 = I[1], MIC_2 = I[2], 
            SOM_1 = I[1], SOM_2 = I[2], SOM_3 = I[1] )

test  <- stode(y = y0, time = 1e6, fun = MIMIC2_ODE, parms = par0, positive = TRUE)
