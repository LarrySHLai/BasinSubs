library(propagate)

#importing data and parameters with their errors ....

###### IMPORT NORTHERN COASTAL RANGE
#for porosity-depth functions
PoFun = read.csv('Backstripping/Pofun.csv')
#for sea level change
deltaSL = read.csv('Sealevel/SL_Unit_N.csv')
#for stratigraphy data
Strat = read.csv('Thickness/Ages_N_new.csv')
#lithology
Litho <- read.csv('Thickness/Litho_N.csv')
#Paleobathymetry
PB <- read.csv('Thickness/PB_Unit_N.csv')

###### FUNCTIONS
#porosity-depth function ....
porosity_depth <- function(depth, porosity_0, decayC) {
  po_result <- porosity_0 * exp(-depth/decayC)
  return(po_result)
}

#Decompaction ...
decompaction <- function(D0, D0_err, Di, Di_err, T0, T0_err, porosity_0, porosity_0_err, decayC, decayC_err) {
  #D0, D0err are the original "depth" and its s.e. of the top of unit i
  #Di, Dierr are the backstripped "depth" and its s.e. of the top of unit i
  #T0, T0err are the un-decompacted thickness and its s.e. of unit i
  #porosity_0, porosity_0err are the surface porosity and its s.e. for certain lithology
  #decayC, decayCerr are the decay constant and its s.e. for certain lithology
  #calculate decompacted thickness 
  K1 <- (-1)*decayC*porosity_0*exp(-(Di/decayC))
  K2 <- (decayC*porosity_0*exp(-(D0/decayC))*(exp(-(T0/decayC))-1)) + T0 - K1 
  decompacted_T <- T0
  repeat{
    test_T <- K1*exp(-(decompacted_T/decayC)) + K2 
    if (abs(test_T-decompacted_T) < 1e-6 | is.na(test_T)) break
    decompacted_T <- test_T
  }
  oD <- c(D0, D0_err)
  dD <- c(Di, Di_err)
  oT <- c(T0, T0_err)
  dT <- c(decompacted_T, abs(decompacted_T-T0)+2*T0_err) #approximate its error as the absolute error plus T0_err
  P <- c(porosity_0, porosity_0_err)
  C <- c(decayC, decayC_err)
  DF <- cbind(oD, dD, oT, dT, P, C)
  
  if (!is.na(Di) & !is.na(Di_err)){
    Ti <- propagate(expr = expression(((-1)*C*P*exp(-(dD/C)))*exp(-(dT/C)) 
                                      + ((C*P*exp(-(oD/C))*(exp(-(oT/C))-1)) 
                                         + oT - (-1)*C*P*exp(-(dD/C))))
                    , data = DF, type = "stat", 
                    do.sim = F, verbose = TRUE)
    Tiextracted <- c(unname(Ti$prop[1]),unname(Ti$prop[3]))
  }
  else{
    Tiextracted <- c(NA, NA)
  }
  #results <- c(decompacted_T, ERR)

  return(Tiextracted)
  #Ouputs: propagated mean1 and se.1
  }

#Mass for unit i at certain time .... 
Mass <- function(Di, Di_err, Ti, Ti_err, rho_w, rho_w_err, rho_g, rho_g_err, porosity_0, porosity_0_err, decayC, decayC_err) {
  dD <- c(Di, Di_err)
  dT <- c(Ti, Ti_err)
  rhoW <- c(rho_w, rho_w_err)
  rhoG <- c(rho_g, rho_g_err)
  P <- c(porosity_0, porosity_0_err)
  C <- c(decayC, decayC_err)
  DF <- cbind(dD, dT, rhoW, rhoG, P, C)
  
  #find the 1-D mass for unit i
  if (!is.na(Di) & !is.na(Di_err)){
    MassS <- propagate(expr = expression(rhoG*dT + (rhoW - rhoG)*C*P*exp(-(dD/C))*(1 - exp(-(dT/C)))) 
                       , data = DF, type = "stat", 
                       do.sim = F, verbose = TRUE)
    MassSextracted <- c(unname(MassS$prop[1]),unname(MassS$prop[3]))
  }
  else{
    MassSextracted <- c(NA,NA)
  }
  return(MassSextracted)
  #Ouputs: propagated meano and se.1
}

#assuming no error for the density of marine water (rho_w = 1025 +- 0)
#Getting density of each unit by dividing the Mass by Ti (but in calculating total)
	
#calculate Di, depth of t0
Depth <- function(Ti, Ti_err) {
  Di <- c()
  Di_err <- c()
  for (i in 1:length(Ti)) {
    #for unit i-1
    if (i == 1){
      Di[i] <- 0
      Di_err[i] <- 0
    }
    else{
      Di[i] <- Di[i-1] + Ti[i-1]
      Di_err[i] <- sqrt(Ti_err[i-1]^2 + Di_err[i-1]^2)
    }
  }
  Diextracted <- data.frame(Di,Di_err)
  return(Diextracted)
}

#Airy Tectonic subsidence function
TectZ <- function(Si, Si_err, Wi, Wi_err, dSL, dSL_err, rho_w, rho_w_err, rho_a, rho_a_err, rho_s, rho_s_err) {
  S <- c(Si, Si_err)
  W <- c(Wi, Wi_err)
  SL <- c(dSL, dSL_err)
  rhoW <- c(rho_w, rho_w_err)
  rhoA <- c(rho_a, rho_a_err)
  rhoS <- c(rho_s, rho_s_err)
  DF <- cbind(S, W, SL, rhoW, rhoA, rhoS)
  
  #find the 1-D mass for unit i
  if (!is.na(rho_s) & !is.na(rho_s_err)){
    Z <- propagate(expr = expression(S*((rhoA-rhoS)/(rhoA - rhoW))+W-SL*(rhoA/(rhoA-rhoW))) 
                   , data = DF, type = "stat", 
                   do.sim = F, verbose = TRUE)
    Zextracted <- c(unname(Z$prop[1]),unname(Z$prop[3]))
  }
  else{
    Zextracted <- c(NA,NA)
  }
  return(Zextracted)
  #Ouputs: propagated meano and se.1
}


###### INITIAL PARAMETERS
#parameters
T0 = Strat$oT
T0_err = Strat$T0_err
age  = Strat$Age
age_err  = Strat$Age_err
unitno = Strat$Units
rho_w = 1025
rho_w_err = 2
rho_a = 3300
rho_a_err = 0
rho_g = c(t(PoFun$rho_g)) #"Mdt","Sst","AvgMS"
rho_g_err = c(t(PoFun$rho_g_err)) #"Mdt","Sst","AvgMS"
porosity_0 = c(t(PoFun$Porosity_0)) #"Mdt","Sst","AvgMS"
porosity_0_err = c(t(PoFun$Porosity_0_err)) #"Mdt","Sst","AvgMS"
decayC = c(t(PoFun$DecayC)) #"Mdt","Sst","AvgMS"
decayC_err = c(t(PoFun$DecayC_err)) #"Mdt","Sst","AvgMS"
dSL = deltaSL$dSL
dSL_err = deltaSL$dSLerr
W = PB$W
W_err = PB$Werr

D0 <- Depth(T0,T0_err)[,1]
D0_err <- Depth(T0,T0_err)[,2]

Stu <- c()
Stu_err <- c()
Su <- c()
Su_err <- c()

MaxElev = 1334/2
MaxElev_err = 1334/2

#non-decompacted total sediment thickness and subsidencebasement depth
for (i in 1:length(unitno)) {
  Stu[i] <- max(D0)-sum(T0[1:i-1])
  Stu_err[i] <- sqrt(sum(T0_err[i:length(unitno)]^2))
  Su[i] <- (max(D0)-sum(T0[1:i-1]))+W[i]-dSL[i]
  Su_err[i] <- sqrt(sum(T0_err[i:length(unitno)]^2)+W_err[i]^2+dSL_err[i]^2)
}
#Decompaction and backstripping ....
N = length(unitno)
for (i in 1:N) {  #for time i to max_i
    #reset tempt vectors
  Ti <- c()
  Ti_err <- c()
  Titempt <- c()
  Ti_k <- c()
  Ti_k_err <- c()
  Mi <- c()
  Mi_err <- c()
  Mtemp <- c()
  Mi_k <- c()
  Mi_k_err <- c()
  if (i == 1) {
    Di <- 0
    Di_err <- 0
    Ti <- c()
    Ti_err <- c()
    Diall <- c()
    Diall_err <- c()
    Tiall <- c()
    Tiall_err <- c()
    Stiall <- c()
    Stiall_err <-c()
    Siall <- c()
    Siall_err <-c()
    rhoiall <- c()
    rhoiall_err <- c()
  }
    for (j in i:N) { #for the unit j
        for (k in 1:length(Litho[1,2:5])) { #for various lithology
              #fraction of certain lithology
              FL = Litho[1,k+1]
              #calculate decompacted thickness
              Titempt = decompaction(D0[j], D0_err[j], Di[j], Di_err[j], T0[j], T0_err[j], porosity_0[k], porosity_0_err[k], decayC[k], decayC_err[k])
              Ti_k[k] <- Titempt[1]*FL
              Ti_k_err[k] <- Titempt[2]*FL
            }
        Ti[j] <- sum(Ti_k, na.rm = TRUE)
        if (i == 1){
          Ti_err[j] <- T0_err[j]
        }
        else {
          Ti_err[j] <- sum(Ti_k_err, na.rm = TRUE)
        }
        
        if (j <= N-1){
          Di[j+1] <- Di[j] + Ti[j]
          Di_err[j+1] <- sqrt(Di_err[j]^2+Ti_err[j]^2)
        }

        for (k in 1:length(Litho[1,2:5])) { #for various lithology
          #fraction of certain lithology
          FL = Litho[1,k+1]
          #calculate decompacted mass
          Mtempt = Mass(Di[j], Di_err[j], Ti[j], Ti_err[j], rho_w, rho_w_err, rho_g[k], rho_g_err[k], porosity_0[k], porosity_0_err[k], decayC[k], decayC_err[k])
          Mi_k[k] <- FL*Mtempt[1]
          Mi_k_err[k] <- FL*Mtempt[2]
        }
        Mi[j] <- sum(Mi_k)
        Mi_err[j] <- sum(Mi_k_err)
        #compile results
    }
  Tiall <- cbind(Tiall, Ti)
  Diall <- cbind(Diall, Di)
  Tiall_err <- cbind(Tiall_err, Ti_err)
  Diall_err <- cbind(Diall_err, Di_err)
  
  #total decompacted thickness
  Stiall <- append(Stiall, max(Di, na.rm = TRUE), after = length(Stiall))
  Stiall_err <-append(Stiall_err, max(Di_err, na.rm = TRUE), after = length(Stiall_err))
  
  #total subsidence
  Siall <- append(Siall, (max(Di, na.rm=TRUE)+W[i]-dSL[i]), after = length(Siall))
  Siall_err <-append(Siall_err, (sqrt(max(Di_err, na.rm=TRUE)^2+W_err[i]^2+dSL_err[i]^2)), after = length(Siall_err))
  
  #compile results
  if (i == N){
    rhoiall <- append(rhoiall, 0 , after = length(rhoiall))
    rhoiall_err <- append(rhoiall_err, 0, after = length(rhoiall_err))
  }
  else{
    rhoiall <- append(rhoiall, sum(Mi, na.rm = TRUE)/max(Di, na.rm = TRUE) , after = length(rhoiall))
    rhoiall_err <- append(rhoiall_err
                          , (sum(Mi, na.rm = TRUE)/max(Di, na.rm = TRUE)) 
                          * sqrt((sqrt(sum(Mi_err^2, na.rm = TRUE))/sum(Mi, na.rm = TRUE))^2 
                                 + (max(Di_err, na.rm = TRUE)/max(Di, na.rm = TRUE))^2) 
                          , after = length(rhoiall_err))
  }
  #reset
  Di <- c()
  Di_err <- c()
  Di[i+1] <- 0
  Di_err[i+1] <- 0
}


#calculate Tectonic subsidence
Zitempt <- c()
Zi <- c()
Zi_err <- c()
N = length(unitno)
for (i in 1:N) { #for the unit j
    #calculate decompacted thickness
    Zitempt = TectZ(Stiall[i], Stiall_err[i], W[i], W_err[i], dSL[i], dSL_err[i], rho_w, rho_w_err, rho_a, rho_a_err, rhoiall[i], rhoiall_err[i])
    Zi[i] <- Zitempt[1]
    Zi_err[i] <- Zitempt[2]
}


#####
#plotting NORTHERN COASTAL RANGE
pdf("Figures/Subsidence_N.pdf", width=3.5, height=3,pointsize = 9)

#background
plot(Su[2:length(Su)]~age[2:length(age)], type="l", xlab="Age (Ma)", ylab="Elevation (m)", main="Northern Coastal Range",xlim=rev(c(0,5)), ylim=rev(c(-2000,9000)))
grid(nx = -2:8, ny = NULL, col = "grey", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
rug(seq(0, 5, by=0.2), ticksize = -0.02, side = 1)
rug(seq(-2000, 9000, by=500), ticksize = -0.02, side = 2)

#paleobathymetry
lines(W[1:length(W)] ~ age[1:length(W)], col="#69CEEC")
arrows(age-age_err, W, age+age_err, W, length=0, angle=0, code=1, lwd=1, col="#69CEEC")
arrows(age, (W-W_err), age, (W+W_err), length=0, angle=0, code=1, lwd=1, col="#69CEEC")
points(W[2:length(W)]~age[2:length(age)], pch=16, cex=.6, col="#69CEEC")
points(W[1]~age[1], pch=7, cex=.6, col="#69CEEC")

#Tectonic subsidence
lines(Zi[2:length(Zi)]~age[2:length(Zi)], col="red")
arrows(age-age_err, Zi, age+age_err, Zi, length=0, angle=0, code=1, lwd=1, col="red")
arrows(age, (Zi-Zi_err), age, (Zi+Zi_err), length=0, angle=0, code=1, lwd=1, col="red")
points(Zi[2:length(Zi)]~age[2:length(age)], pch=16, cex=.6, col="red")
points(Zi[1]~age[1], pch=7, cex=.6, col="red")
lines(c(age[1:2]),c(Zi[1:2]),lty = "longdash",col="red")

#un-decompacted curve
arrows(age-age_err, Su, age+age_err, Su, length=0, angle=0, code=1, lwd=1)
arrows(age, (Su-Su_err), age, (Su+Su_err), length=0, angle=0, code=1, lwd=1)
points(Su[2:length(Su)]~age[2:length(age)], pch=16, cex=.6)
points(Su[1]~age[1], pch=7, cex=.6)
lines(c(age[1:2]),c(Su[1:2]),lty = "longdash")

#plot total subsidence curve
points(Siall[2:length(Siall)]~age[2:length(age)], type="l", col="blue")
arrows(age[1:(length(Siall)-1)]-age_err[1:(length(Siall)-1)], Siall[1:(length(Siall)-1)], age[1:(length(Siall)-1)]+age_err[1:(length(Siall)-1)], Siall[1:(length(Siall)-1)], length=0, angle=0, code=3, lwd=1, col="blue")
arrows(age[1:(length(Siall)-1)], (Siall[1:(length(Siall)-1)]-Siall_err[1:(length(Siall)-1)]), age[1:(length(Siall)-1)], (Siall[1:(length(Siall)-1)]+Siall_err[1:(length(Siall)-1)]), length=0, angle=0, code=1, lwd=1, col="blue")
points(Siall[2:(length(Siall)-1)]~age[2:(length(age)-1)], pch=16, cex=.6, col="blue")
points(Siall[1]~age[1], pch=7, cex=.6, col="blue")
lines(c(age[1:2]),c(Siall[1:2]),lty = "longdash", col="blue")
lines(c(0,age[1]),c(-MaxElev, Siall[1]),lty = "dashed", col="blue")
lines(c(0,age[2]),c(-MaxElev, Siall[2]),lty = "dashed",col="#b8baff")

#Assumed basement uplift history
arrows(c(mean(c(5.57,4.37)),mean(c(4.31,3.47))), c(64.262-100,0-100), c(mean(c(5.57,4.37)),mean(c(4.31,3.47))), c(64.262+100,0+100), length=0, angle=0, code=3, lwd=1,col="grey")
arrows(c(mean(c(5.57,4.37))-(5.57-4.37)/2,mean(c(4.31,3.47))-(4.31-3.47)/2), c(64.262,0), c(mean(c(5.57,4.37))+(5.57-4.37)/2,mean(c(4.31,3.47))+(4.31-3.47)/2), c(64.262,0), length=0, angle=0, code=3, lwd=1,col="grey")
points(c(64.262,0)~c(mean(c(5.57,4.37)),mean(c(4.31,3.47))), pch=2, cex=.6,col="grey")
lines(c(mean(c(5.57,4.37)),mean(c(4.31,3.47))),c(64.262, 0),lty = "dashed",col="#808080")

#MAX ELEVATION
points( c(-MaxElev) ~ c(0), pch=1, cex=.6)
arrows(c(0), -(MaxElev-MaxElev_err), c(0), -(MaxElev+MaxElev_err), length=0, angle=0, code=3, lwd=1)

#Legend
legend(5, 5500, legend=c("Non-decompacted total subsidence", "Decompacted total subsidence", "Tectonic subsidence","Paleobathymetry"), col=c("black","blue","red","#69CEEC"), lty=1, cex=.6, lwd=1.5)

#calculating uplift rates
MaxUrate = ((Siall[1]+MaxElev)/age[1])/1000
MaxUrate_err = MaxUrate*sqrt((sqrt(Siall_err[1]^2+MaxElev_err^2)/(Siall[1]+MaxElev))^2 + (age_err[1]/age[1])^2)
MinUrate = ((Siall[2]+MaxElev)/age[2])/1000
MinUrate_err = MinUrate*sqrt((sqrt(Siall_err[2]^2+MaxElev_err^2)/(Siall[2]+MaxElev))^2 + (age_err[2]/age[2])^2)
text(0.15, 5000, paste(format(MaxUrate,digits=3) ,"\u00b1",format(MaxUrate_err,digits=2)), col="blue", cex = .6)
text(.15, 5500, paste("mm/yr"), col="blue", cex = .6)
text(.6, 1200, paste(format(MinUrate,digits=2) ,"\u00b1",format(MinUrate_err,digits=1)), col="#b8baff", cex = .6)
text(.6, 1700, paste("mm/yr"), col="#b8baff", cex = .6)

#calculating subsidence rates
Srate_2 = ((Siall[2]-Siall[15])/abs(age[2]-age[15]))/1000
Srate_2_err = Srate_2 * sqrt((sqrt(Siall_err[2]^2+Siall_err[15]^2)/(Siall[2]-Siall[15]))^2 + (sqrt(age_err[2]^2+age_err[15]^2)/(age[2]-age[15]))^2)
text(2.2, 3500, paste(format(Srate_2,digits=2) ,"\u00b1",format(Srate_2_err,digits=1)), col="blue", cex = .6)
text(2.2, 4000, paste("mm/yr"), col="blue", cex = .6)

dev.off()

######
#output results
Results <- data.frame(Unit = unitno, Age = age, Age_err = age_err, T0 = T0, T0_err = T0_err,
                      PB = W, PB_err = W_err, dSL = dSL, dSLerr = dSL_err
                      , undecompactedTotalS = Su, undecompactedTotalS_err = Su_err
                      , TotalS = Siall, TotalS_err = Siall_err
                      , TectS = Zi, TectS_err = Zi_err, RhoS = rhoiall, RhoS_err = rhoiall_err)
print(Results)
write.csv(Results, "Subsidence results_N.csv", row.names = FALSE)

Results2 <- data.frame(SubRate = Srate_2, SubRate_err = Srate_2_err,
                      Elev = MaxElev, Elev_err = MaxElev_err, 
                      MinUpRate = MinUrate, MinUpRate_err = MinUrate_err,
                      EstUpRate = MaxUrate, EstUpRate_err = MaxUrate_err)
print(Results2)
write.csv(Results2, "vertical rates_N.csv", row.names = FALSE)
