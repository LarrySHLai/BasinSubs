#importing others' functions
#for sand
Porosity = read.csv('Backstripping/Porosity_sand.csv')
#for mud
Porosity = read.csv('Backstripping/Porosity_mud.csv')

decayC = Porosity$Decay_c
porosity_0 = Porosity$InitialPorosity
min_z = Porosity$min_z
max_z = Porosity$max_z

#porosity-depth function
porosity_depth <- function(depth, porosity_0, decayC) {
  po_result <- porosity_0 * exp(-depth/decayC)
  return(po_result)
}

#create all functions
for (i in 1:length(porosity_0)) {
  if (i < 2) {
    z = seq(min_z[i],max_z[i],50)
    porosity_z = porosity_depth(z, porosity_0[i], decayC[i])
  } 
  else {
    z_input = seq(min_z[i],max_z[i],50)
    z = append(z, z_input, after = length(z))
    porosity_output = porosity_depth(z_input, porosity_0[i], decayC[i])
    porosity_z = append(porosity_z, porosity_output, after = length(porosity_z))
  }
}
model <- lm(log(porosity_z)~z) 
summary(model)
y.fitted <- model$fitted.values

Meanporosity_0 <- exp(coef(model)[1])
SEporosity_0 <- Meanporosity_0 * coef(summary(model))[1, 2]
MeandecayC <- -1/coef(model)[2]
SEdecayC <- MeandecayC * sqrt((coef(summary(model))[2, 2]/coef(model)[2])^2)

sse <- sum((log(porosity_z) - y.fitted)^2)
mse <- sse / (length(log(porosity_z)) - 2)

t.val <- qt(0.975, length(log(porosity_z)) - 2)

se <- sqrt(mse * sqrt(1 / length(log(porosity_z)) + (z - mean(z))^2 / sum((z - mean(z))^2)))

CIupper <- data.frame(po = exp(suppressWarnings(y.fitted + t.val * se)), z = z)
CIupper <- CIupper[order(CIupper[,2]),]
CIlower <- data.frame(po = exp(suppressWarnings(y.fitted - t.val * se)), z = z)
CIlower <- CIlower[order(CIlower[,2]),]

plot(z/1000~porosity_z, ylim=rev(c(0,6)), xlim=c(0,0.8), xaxt='n',pch=46, col='#6e6e6e', xlab='Porosity', ylab='Depth from Sea Level (km)', main="Carbonate")
axis(1, xlim = c(0,.8))
grid(nx = NULL, ny = 0:0.8, col = "grey", lty = 1,lwd = .5, equilogs = TRUE)
lines(Meanporosity_0*exp(-seq(0,6000)/MeandecayC), seq(0,6000)/1000, col='red',lwd=2)
lines(CIlower$po, CIlower$z/1000, col="blue",lwd=1.5)
lines(CIupper$po, CIupper$z/1000, col="blue",lwd=1.5)
legend(.25, 4.5, legend=c("Sourced functions", "Mean function", "95% confi. interv."), col=c("#6e6e6e","red","blue"), lty=c(3,1,1), cex=.9, lwd=c(1, 2, 1.5))


#for mean marine sediments curve
Meanporosity_0 <- 0.663
SEporosity_0 <- 0.03315
MeandecayC <- 1333
SEdecayC <- 66.65

z = seq(0, 6000, by = 20)
porosity_z = Meanporosity_0 * exp(-z/MeandecayC)

plot(z/1000~porosity_z, ylim=rev(c(0,6)), xlim=c(0,0.8), type='l', col='red', lwd=2, xlab='Porosity', ylab='Depth from Sea Level (km)', main="Avgerage Marine Sediments")
axis(1, xlim = c(0,.8))
grid(nx = NULL, ny = 0:0.8, col = "grey", lty = 1,lwd = .5, equilogs = TRUE)
legend(.3, 5, legend=c("Mean function"), col=c("red"), lty=c(1), cex=.9, lwd=c(2))

