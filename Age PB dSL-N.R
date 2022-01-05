#for stratigraphy data
Strat = read.csv('Thickness/Ages_N.csv')

age  = c(Strat$age[3:14])
age_err  = c(Strat$age_err[3:14])
top = c(Strat$mean_T[3:14])
top_err = c(Strat$mean_T_err[3:14])

a = age[1:7]
t = top[1:7]
upperREG <- lm(a ~ t)
plotupperREG <- lm(t ~ a)
a = age[7:8]
t = top[7:8]
upperREG2 <- lm(a ~ t)
plotupperREG2 <- lm(t ~ a)
a = age[8:11]
t = top[8:11]
lowerREG2 <- lm(a ~ t)
plotlowerREG2 <- lm(t ~ a)
a = age[11:12]
t = top[11:12]
lowerREG <- lm(a ~ t)
plotlowerREG <- lm(t ~ a)

#abline(lowerREG)
#abline(upperREG)
summary(upperREG)
summary(upperREG2)
summary(lowerREG2)
summary(lowerREG)

#predict strat top and eroded top
xinput = c(Strat$mean_T[1:2])
pAge_top = predict(upperREG, data.frame(t = xinput), se.fit=T, interval = 'confidence')
Strat$age[1:2] = pAge_top$fit[1:2]
Strat$age_err[1:2] = pAge_top$se.fit

#predict strat top and eroded top
xinput = c(Strat$mean_T[15])
pAge_base = predict(lowerREG, data.frame(t = xinput), se.fit=T, interval = 'confidence')
Strat$age[15] = pAge_base$fit[1]
#Strat$age_err[15] = pAge_base$se.fit
Strat$age_err[15] = Strat$age[15]*.05 #just give 5% error 

#####
#save following plot
pdf("Figures/AgeModel_N.pdf", width=3.5, height=3,pointsize = 9)

#plot data and error bar
plot(top~age, pch = 16,cex=.7, ylab="Stratigraphic Height (m)", xlab="Age (Ma)", main="Northern Coastal Range",ylim=(c(0,8000)), xlim=rev(c(0,5)))
rug(seq(0, 5, by=.5), ticksize = -0.02, side = 1)
rug(seq(0, 8000, by=500), ticksize = -0.02, side = 2)
arrows(age+age_err, top, age-age_err, top, length=0, angle=0, code=3, col ="grey", lwd=1)
arrows(age, top+top_err, age, top-top_err, length=0, angle=0, code=3, col ="grey", lwd=1)
points(top~age,pch=16,cex=.7)
#plot regression lines
abline(plotupperREG)
abline(plotupperREG2)
abline(plotlowerREG2)
abline(plotlowerREG)

#plot predicted data (upper)
xinput = c(Strat$mean_T[c(2,15)])
arrows(Strat$age[c(2,15)]+Strat$age_err[c(2,15)], Strat$mean_T[c(2,15)], Strat$age[c(2,15)]-Strat$age_err[c(2,15)], Strat$mean_T[c(2,15)], length=0, angle=0, code=3, col ="red", lwd=1)
arrows(Strat$age[c(2,15)], Strat$mean_T[c(2,15)]+Strat$mean_T_err[c(2,15)], Strat$age[c(2,15)], Strat$mean_T[c(2,15)]-Strat$mean_T_err[c(2,15)], length=0, angle=0, code=3, col ="red", lwd=1)
points(xinput~Strat$age[c(2,15)],col="red",cex=.7)
#plot predicted data (lower)
xinput = c(Strat$mean_T[1])
arrows(Strat$age[1]+Strat$age_err[1], Strat$mean_T[1], Strat$age[1]-Strat$age_err[1], Strat$mean_T[1], length=0, angle=0, code=3, col ="red", lwd=1)
arrows(Strat$age[1], Strat$mean_T[1]+Strat$mean_T_err[1], Strat$age[1], Strat$mean_T[1]-Strat$mean_T_err[1], length=0, angle=0, code=3, col ="red", lwd=1)
points(xinput~Strat$age[1],pch=7,col="red",cex=.7)

dev.off()

#######
#thickness....
top = Strat$mean_T
top_err = Strat$mean_T_err
T0 <- c(c(Strat$mean_T[1:14] - Strat$mean_T[2:15]),0)
T0err <- c(c(sqrt((Strat$mean_T_err[1:14])^2 + (Strat$mean_T_err[2:15])^2)), Strat$mean_T_err[15])

#Output data
age  = Strat$age[1:15]
age_err  = Strat$age_err[1:15]
unitno = Strat$unit[1:15]

StratUnits <- data.frame(Units = unitno, Age = age, Age_err = age_err, oT = T0, T0_err = T0err)
print(StratUnits)

write.csv(StratUnits, "Thickness/Ages_N_new.csv", row.names = FALSE)

#Find interceptions among regressions
REGintercept1 <- abs((unname(coef(upperREG2)[1])-unname(coef(upperREG)[1]))/(unname(coef(upperREG2)[2])-unname(coef(upperREG)[2])))
REGintercept2 <- abs((unname(coef(lowerREG2)[1])-unname(coef(upperREG2)[1]))/(unname(coef(lowerREG2)[2])-unname(coef(upperREG2)[2])))
REGintercept3 <- abs((unname(coef(lowerREG2)[1])-unname(coef(lowerREG)[1]))/(unname(coef(lowerREG2)[2])-unname(coef(lowerREG)[2])))
plotREGintercept1 <- abs((unname(coef(plotupperREG2)[1])-unname(coef(plotupperREG)[1]))/(unname(coef(plotupperREG2)[2])-unname(coef(plotupperREG)[2])))
plotREGintercept2 <- abs((unname(coef(plotlowerREG2)[1])-unname(coef(plotupperREG2)[1]))/(unname(coef(plotlowerREG2)[2])-unname(coef(plotupperREG2)[2])))
plotREGintercept3 <- abs((unname(coef(plotlowerREG2)[1])-unname(coef(plotlowerREG)[1]))/(unname(coef(plotlowerREG2)[2])-unname(coef(plotlowerREG)[2])))

#output regression functions
Regressions <- data.frame(Model = c("upper ages","upper ages 2","lower ages 2","lower ages")
                         , slope = c(unname(coef(upperREG)[2]),unname(coef(upperREG2)[2]),unname(coef(lowerREG2)[2]),unname(coef(lowerREG)[2]))
                         , slope_err = c(coef(summary(upperREG))[2, "Std. Error"],coef(summary(upperREG2))[2, "Std. Error"],coef(summary(lowerREG2))[2, "Std. Error"],coef(summary(lowerREG))[2, "Std. Error"])
                         , intercept = c(unname(coef(upperREG)[1]),unname(coef(upperREG2)[1]),unname(coef(lowerREG2)[1]),unname(coef(lowerREG)[1]))
                         , intercept_err = c(coef(summary(upperREG))[1, "Std. Error"],coef(summary(upperREG2))[1, "Std. Error"], coef(summary(lowerREG2))[1, "Std. Error"],coef(summary(lowerREG))[1, "Std. Error"])
                         , r_squared = c(summary(upperREG)$r.squared,summary(upperREG2)$r.squared, summary(lowerREG2)$r.squared,summary(lowerREG)$r.squared)
                         , p_slope = c(coef(summary(upperREG))[2, "Pr(>|t|)"],coef(summary(upperREG2))[2, "Pr(>|t|)"],coef(summary(lowerREG2))[2, "Pr(>|t|)"],coef(summary(lowerREG))[2, "Pr(>|t|)"])
                         , lowerupper = c(REGintercept1, REGintercept2, REGintercept3,0)
                         , upperbound = c(max(top), REGintercept1, REGintercept2, REGintercept3))

print(Regressions)
write.csv(Regressions, "Thickness/Age_DepthModel_N.csv", row.names = FALSE)

plotRegressions <- data.frame(Model = c("upper ages","upper ages 2","lower ages 2","lower ages")
                          , slope = c(unname(coef(plotupperREG)[2]),unname(coef(plotupperREG2)[2]),unname(coef(plotlowerREG2)[2]),unname(coef(plotlowerREG)[2]))
                          , slope_err = c(coef(summary(plotupperREG))[2, "Std. Error"],coef(summary(plotupperREG2))[2, "Std. Error"],coef(summary(plotlowerREG2))[2, "Std. Error"],coef(summary(plotlowerREG))[2, "Std. Error"])
                          , intercept = c(unname(coef(plotupperREG)[1]),unname(coef(plotupperREG2)[1]),unname(coef(plotlowerREG2)[1]),unname(coef(plotlowerREG)[1]))
                          , intercept_err = c(coef(summary(plotupperREG))[1, "Std. Error"],coef(summary(plotupperREG2))[1, "Std. Error"], coef(summary(plotlowerREG2))[1, "Std. Error"],coef(summary(plotlowerREG))[1, "Std. Error"])
                          , r_squared = c(summary(plotupperREG)$r.squared,summary(plotupperREG2)$r.squared, summary(plotlowerREG2)$r.squared,summary(plotlowerREG)$r.squared)
                          , p_slope = c(coef(summary(plotupperREG))[2, "Pr(>|t|)"],coef(summary(plotupperREG2))[2, "Pr(>|t|)"],coef(summary(plotlowerREG2))[2, "Pr(>|t|)"],coef(summary(plotlowerREG))[2, "Pr(>|t|)"])
                          , lowerupper = c(plotREGintercept1, plotREGintercept2, plotREGintercept3,max(age))
                          , upperbound = c(0, plotREGintercept1, plotREGintercept2, plotREGintercept3))

print(plotRegressions)
write.csv(plotRegressions, "Thickness/Depth_AgeModel_N.csv", row.names = FALSE)



########################
PBimport = read.csv('Thickness/Bathymetry_N.csv')

allStratH = c(PBimport$StratH)
allPB <- c(PBimport$W)
allPBerr <- c(PBimport$Werr)

PBage <- c()
PBageerr <- c()

for (i in 1:length(allStratH)) {
if (allStratH[i] >= REGintercept1 & !is.na(allStratH[i])) {
  xinput = allStratH[i]
  PredictedAge <- predict(upperREG, data.frame(t = xinput), se.fit=T, interval = 'confidence')
  PBage[i] <- unname(PredictedAge$fit[,1])
  PBageerr[i] <- unname(PredictedAge$se.fit)
}
else if (allStratH[i] < REGintercept1 & allStratH[i] > REGintercept2 & !is.na(allStratH[i])) {
    xinput = allStratH[i]
    PredictedAge <- predict(upperREG2, data.frame(t = xinput), se.fit=T, interval = 'confidence')
    PBage[i] <- unname(PredictedAge$fit[,1])
    #PBageerr[i] <- unname(PredictedAge$se.fit)
    PBageerr[i] <- PBage[i]*0.05 #just give 5% error 
}
else if (allStratH[i] <= REGintercept2 & allStratH[i] >= REGintercept3 & !is.na(allStratH[i])) {
    xinput = allStratH[i]
    PredictedAge <- predict(lowerREG2, data.frame(t = xinput), se.fit=T, interval = 'confidence')
    PBage[i] <- unname(PredictedAge$fit[,1])
    PBageerr[i] <- unname(PredictedAge$se.fit)
  }
else if (allStratH[i] < REGintercept3 & !is.na(allStratH[i])) {
  xinput = allStratH[i]
  PredictedAge <- predict(lowerREG, data.frame(t = xinput), se.fit=T, interval = 'confidence')
  PBage[i] <- unname(PredictedAge$fit[,1])
  #PBageerr[i] <- unname(PredictedAge$se.fit)
  PBageerr[i] <- PBage[i]*0.05 #just give 5% error 
}
else {
    PBage[i] <- NA
    PBageerr[i] <- NA
  }
}

AgePB <- data.frame(Age = PBage, Age_err = PBageerr, W = allPB, Werr = allPBerr, stratH = allStratH, section = PBimport$section, ID = PBimport$id, source = PBimport$Source)
AgePB_r <- AgePB[order(AgePB[,1]),]

print(AgePB_r)
write.csv(AgePB_r, "Thickness/PB_Age_N.csv", row.names = FALSE)

#calculate Paleobathymetry for each unit
ExtractedPB <- c()
ExtractedPBerr <- c()
MeanPB <- c()
sePB <- c()

StratUnits = read.csv('Thickness/Ages_N_new.csv')
age  = c(StratUnits$Age)
age_err  = c(StratUnits$Age_err)
unitno = c(StratUnits$Units)


for (i in 1:length(unitno)) {
    for (j in 1:length(AgePB$Age)) {
      if (abs(AgePB$Age[j] - age[i]) <= 0.2 & !is.na(AgePB$Age[j]) & !is.na(AgePB$W[j]) & AgePB$Werr[j] <= 4500) {
        ExtractedPB <- append(ExtractedPB, AgePB$W[j], after = length(ExtractedPB)) 
        ExtractedPBerr <- append(ExtractedPBerr, AgePB$Werr[j], after = length(ExtractedPBerr))
      }
      #already exclude PB data with error > 4500 (the detective limit of the logistic model by Haywood et al, 2016)
      MeanPB[i] <- mean(ExtractedPB, na.rm = TRUE)
      sePB[i] <- sqrt(sum(ExtractedPBerr^2))/length(ExtractedPBerr)
    }
  ExtractedPB <- c()
  ExtractedPBerr <- c()
}

sePB[sePB == "NaN"] = NA 

#arbitrarily make the last guessed paleobathymetery as the the half of the PB
MeanPB[1] <- MeanPB[3]/2
sePB[1] <- MeanPB[3]/2

print(MeanPB)
print(sePB)

UnitPB <- data.frame(Unit = unitno, W = MeanPB, Werr = sePB)
print(UnitPB)
write.csv(UnitPB, "Thickness/PB_Unit_N.csv", row.names = FALSE)


#####
#plot paleobathymetry-age diagram
pdf("Figures/PB_Age_N.pdf", width=3.5, height=2.25,pointsize = 7.5)

plot(AgePB_r$W~AgePB_r$Age, xlim=rev(c(0,5)), ylim=rev(c(0,3000)), pch=16, cex=.7, col="grey",xlab="Age (Ma)", ylab="Paleobathymetry, W (m)",main="N Coastal Ragne Paleobathymetry")
rug(seq(0, 5, by=0.5), ticksize = -0.02, side = 1)
rug(seq(0, 3000, by=500), ticksize = -0.02, side = 2)
arrows(AgePB_r$Age, (AgePB_r$W-AgePB_r$Werr), AgePB_r$Age, (AgePB_r$W+AgePB_r$Werr), length=0.01, angle=0, code=3, col ="grey")
arrows(AgePB_r$Age-AgePB_r$Age_err, AgePB_r$W, AgePB_r$Age+AgePB_r$Age_err, AgePB_r$W, length=0.01, angle=0, code=3, col ="grey")
#points(AgePB_r$W~AgePB_r$Age,pch=19, cex=.6)

points(MeanPB[2:length(MeanPB)]~age[2:length(age)], pch=16, cex=.8, col="red")
points(MeanPB[1]~age[1], pch=7,cex=.7,col="red")
arrows(age-age_err, MeanPB, age+age_err, MeanPB, length=0, angle=0, code=3, col ="red")
arrows(age, (MeanPB-sePB), age, (MeanPB+sePB), length=0, angle=0, code=3, col ="red")
#points(MeanPB~age,pch=19, cex=.6)

dev.off()


###########################
#deltaSL
Sealevelcurve = read.csv('Sealevel/SeaLevel_Miller2020.csv')
SL <- Sealevelcurve$SeaLevel
SLAge <- Sealevelcurve$Age.Ma.

StratUnits = read.csv('Thickness/Ages_N_new.csv')
age  = c(StratUnits$Age)
age_err  = c(StratUnits$Age_err)
unitno = c(StratUnits$Units)

ExtractedSL <- c()
MeanSL <- c()
seSL <- c()

for (i in 1:length(unitno)) {
  #for unit i-1
  for (j in 1:length(SL)) {
    if (SLAge[j] >= age[i]-age_err[i]-0.01 & SLAge[j] <= age[i]+age_err[i]+0.01) {
      ExtractedSL <- append(ExtractedSL, SL[j], after = length(ExtractedSL)) 
    }
    #already exclude PB data with error > 4500 (the detective limit of the logistic model by Haywood et al, 2016)
    MeanSL[i] <- mean(ExtractedSL)
    seSL[i] <- sd(ExtractedSL)/sqrt(length(ExtractedSL))
  }
  ExtractedSL <- c()
}

UnitSL <- data.frame(Unit = unitno, dSL = MeanSL, dSLerr = seSL)
print(UnitSL)
write.csv(UnitSL, "Sealevel/SL_Unit_N.csv", row.names = FALSE)


#####
#plot deltaSL-age diagram
pdf("Figures/SL_Age_N.pdf", width=3.5, height=2.25, pointsize = 6)

plot(SL~SLAge, xlim=rev(c(0,5)), ylim=c(-150,50),type="l",col="grey", xlab="Age (Ma)", 
     ylab=expression(paste(Delta, "Sea Level (m)")),
     main=expression(paste(Delta, "Sea Level (Miller et al., 2020)")))
rug(seq(0, 5, by=0.5), ticksize = -0.02, side = 1)
rug(seq(-150, 50, by=25), ticksize = -0.02, side = 2)
#lines(deltaSL$name~deltaSL$age)

points(MeanSL[2:length(MeanSL)]~age[2:length(age)], pch=16,cex=1, col="red")
points(MeanSL[1]~age[1], pch=7,cex=1,col="red")
arrows(age-age_err, MeanSL, age+age_err, MeanSL, length=0, angle=0, code=3, col ="red")
arrows(age, MeanSL-seSL, age, MeanSL+seSL, length=0, angle=0, code=3, col ="red")

dev.off()
