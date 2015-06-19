makeInitalPar <- function(I, fMET, TSOI, fCLAY, Tao_MOD1, VMAX=NULL, KM=NULL){

##################
#Calculate Vmax (using parameters from German 2012, as in Wieder et al. 2013 Nature Climate Change)
##################
Vslope   <- array(0.063,dim=6)
Vint     <- 5.47
aV       <- 8e-6
Vmax     <- exp(TSOI * Vslope + Vint) * aV
MOD1     <- c(10, 2, 10, 3, 3, 2) #Modifiers for MIMIC2_b
if(is.null(VMAX)){
  VMAX     <- Vmax * MOD1 
}else{
  if(length(VMAX) != length( Vmax * MOD1 )) stop('VMAX is wrong length')
}

################
#Calculate KM
################
Kslope   <- c(0.017, #META LIT to MIC_1
              0.027, #STRU LIT to MIC_1 
              0.017,  #AVAI SOM to MIC_1 
              0.017, #META LIT to MIC_2
              0.027, #STRU LIT to MIC_2
              0.017) #AVAI SOM to MIC_2
Kint     <- 3.19
aK       <- 10
Km       <- exp(Kslope * TSOI + Kint) * aK
k        <- 2.0    #2.0  		#REDUCED FROM 3 TO 1, REDUCES TEXTURE EFFECTS ON SOMa decay
a        <- 2.0    #2.2			#increased from 4.0 to 4.5
pSCALAR  <- a * exp(-k*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
MOD2     <- c( 8, 2 ,4 * pSCALAR, 2, 4, 6 * pSCALAR)  ##Modifiers for MIMICS2_b
if(is.null(KM)){
  KM       <- Km / MOD2
}else{
  if(length(KM) != length( Km/ MOD1 )) stop('KM is wrong length')
}

########################
## CUE - carbon use effiency
########################
CUE        <- c(0.55, 0.25, 0.75, 0.35)  #for LITm and LITs entering MICr and MICK, respectively

########################
## fPHYS
########################
fPHYS    <- c(0.3 * exp(1.3*fCLAY), 0.2 * exp(0.8*fCLAY))   #fraction to SOMp

########################
## fCHEM
########################
fCHEM    <- c(0.1 * exp(-3*fMET)  , 0.3 * exp(-3*fMET)  ) 	#fraction to SOMc
fAVAI    <- 

########################
## tao
########################
tao      <- c(5.2e-4*exp(0.3*fMET), 2.4e-4*exp(0.1*fMET))  
tao      <- tao * Tao_MOD1

#######################
## desorb
#######################
desorb   <- 9e-4 * exp(-3*(sqrt(fCLAY))) #if modified by MIC!
desorb   <- 3e-4 * exp(-4*(sqrt(fCLAY))) #if stand alone rate constant
desorb   <- 1.5e-5 * exp(-1.5*(fCLAY))      #CHANGED FOR GLOBAL RUN!!!  

return(c( I = I, VMAX = VMAX, KM = KM, CUE = CUE, 
            fPHYS = fPHYS, fCHEM = fCHEM, 
            fAVAI = 1- (fPHYS + fCHEM), FI = c(0.05, 0.05), 
            tao = tao, 
            LITmin = rep(NA, dim=4), SOMmin = rep(NA, dim=2), MICtrn = rep(NA, dim=6), 
            desorb = desorb, 
            DEsorb = rep(NA, dim=1), OXIDAT = rep(NA, dim=1), KO = c(4,4) ))
}