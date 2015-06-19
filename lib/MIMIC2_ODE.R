MIMIC2_ODE <- function(t, y, pars) {
  with (as.list(c(y, pars)),{
    
    #Flows to and from MIC_1
    LITmin[1] = MIC_1 * VMAX1 * LIT_1 / (KM1 + LIT_1)   #MIC_1 decomp of MET lit
    LITmin[2] = MIC_1 * VMAX2 * LIT_2 / (KM2 + LIT_2)   #MIC_1 decomp of STRUC lit
    MICtrn[1] = MIC_1 * tao1  * fPHYS1                  #MIC_1 turnover to PHYSICAL SOM 
    MICtrn[2] = MIC_1 * tao1  * fCHEM1                  #MIC_1 turnover to CHEMICAL SOM  
    MICtrn[3] = MIC_1 * tao1  * fAVAI1                  #MIC_1 turnover to AVAILABLE SOM  
    SOMmin[1] = MIC_1 * VMAX3 * SOM_3 / (KM3 + SOM_3)   #decomp of SOMa by MIC_1
    
    #Flows to and from MIC_2
    LITmin[3] = MIC_2 * VMAX4 * LIT_1 / (KM4 + LIT_1)   #decomp of MET litter
    LITmin[4] = MIC_2 * VMAX5 * LIT_2 / (KM5 + LIT_2)   #decomp of SRUCTURAL litter
    MICtrn[4] = MIC_2 * tao2  * fPHYS2                  #MIC_2 turnover to PHYSICAL  SOM 
    MICtrn[5] = MIC_2 * tao2  * fCHEM2                  #MIC_2 turnover to CHEMICAL  SOM  
    MICtrn[6] = MIC_2 * tao2  * fAVAI2                  #MIC_2 turnover to AVAILABLE SOM  
    SOMmin[2] = MIC_2 * VMAX6 * SOM_3 / (KM6 + SOM_3)   #decomp of SOMa by MIC_2
    
    DEsorb    = SOM_1 * desorb  #* (MIC_1 + MIC_2)		#desorbtion of PHYS to AVAIL (function of fCLAY)
    OXIDAT    = ((MIC_2 * VMAX5 * SOM_2 / (KO2*KM5 + SOM_2)) +
                   (MIC_1 * VMAX2 * SOM_2 / (KO1*KM2 + SOM_2)))  #oxidation of C to A
    #can make fluxes from CHEM a function of microbial biomass size?
    
    dLIT_1 = I1*(1-FI1) - LITmin[1] - LITmin[3]
    dMIC_1 = CUE1*(LITmin[1]+ SOMmin[1]) + CUE2*(LITmin[2]) - sum(MICtrn[1:3])
    dSOM_1 = I1*FI1 + MICtrn[1] + MICtrn[4]- DEsorb 
    
    dLIT_2 = I2 * (1-FI2) - LITmin[2] - LITmin[4]
    dMIC_2 = CUE3*(LITmin[3]+ SOMmin[2]) + CUE4*(LITmin[4]) - sum(MICtrn[4:6])  
    dSOM_2 = I2*FI2 + MICtrn[2] + MICtrn[5] - OXIDAT
    
    dSOM_3  = MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]
    
    list(c(dLIT_1, dLIT_2, dMIC_1, dMIC_2, dSOM_1, dSOM_2, dSOM_3))
  })
}
