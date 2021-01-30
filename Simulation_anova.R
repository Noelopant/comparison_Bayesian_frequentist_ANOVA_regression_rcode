rm(list = ls())
library(ggplot2)
library(sn)
library(BayesFactor)
library(stats)
library(sjstats)
library(car)
library(gridExtra)


#####


nrepeats <- 100
simulationOptions <- expand.grid(
  n            = c(150),
  distribution = c(1, 30),
  variancesequal = c(1, 2),
  repeats      = 1:nrepeats,
  stringsAsFactors = FALSE
)

simulationOptions

sigma   <- 2                                        # standard deviation
b0      <- 0                                        # intercept
b1      <- 0                                        # level predictor 1
b2      <- 0.5                                      # level predictor 2
b3      <- 1                                        # level predictor 3



simulateData <- function(n, t, var){

  epsilon1 <- rt(n,t)
  epsilon2 <- rt(n,t)*var
  epsilon3 <- rt(n,t)*(var*var)

  A <- c(rep(1, n/3), rep(0, n/3), rep(0, n/3))       # predictor 1
  B <- c(rep(0, n/3), rep(1, n/3), rep(0, n/3))       # predictor 2 
  C <- c(rep(0, n/3), rep(0, n/3), rep(1, n/3))       # predictor 3
 
  y <- b0 + (A * (b1 + epsilon1)) + (B * (b2 + epsilon2)) + (C * (b3 + epsilon3))       # outcome
  
  predictor <- factor(A + 2*B + 3*C)
  
  df = data.frame(
    "predictor" = predictor, 
    "outcome" = y
  )
  return(df)
}

fitFreq      <- function(data){
  freq = aov(data$outcome ~ data$predictor, data=data)
  return(summary(freq)[[1]][1,4])
}

fitBayes     <- function(data){
  bf = anovaBF(outcome ~ predictor, data=data)
  return(extractBF(bf)$bf)
}



nsim <- nrow(simulationOptions)

length(seq_len(nsim))
nrow(data)
nrow(simulationOptions)

resf <- rep(NA, length(seq_len(nsim)))
resb <- rep(NA, length(seq_len(nsim)))

for (i in seq_len(nsim)) {
  
  set.seed(i)
  opts <- simulationOptions[i, ]
  
  data <- simulateData(opts[[1]][1],opts[[2]][1], opts[[3]][1])
  resultFreq  <- fitFreq(data)
  resultBayes <- fitBayes(data)
  
  resf[i] <- resultFreq
  resb[i] <- resultBayes
  
  
  print((100/(length(seq_len(nsim))))*i)
}

results <- data.frame(
  "n" = simulationOptions$n, 
  "eqvar" = simulationOptions$variancesequal,
  "dist" = simulationOptions$distribution,
  "f" = resf, 
  "BF" = resb)


#####
# Make histograms of r2 and bayes factors

selection1 <- results$f[results$n == 150 & results$eqvar == 1 &results$dist == 30]

p1 <- ggplot()+
  aes(selection1)+
  geom_histogram()+
  xlab("All assumptions met")

selection2 <- results$f[results$n == 150 & results$eqvar == 2 &results$dist == 30]

p2 <- ggplot()+
  aes(selection2)+
  geom_histogram()+
  xlab("Equality of variances violated")


selection3 <- results$f[results$n == 150 & results$eqvar == 1 &results$dist == 1]

p3 <- ggplot()+
  aes(selection3)+
  geom_histogram()+
  xlab("Normal distribution violated")

selection4 <- results$f[results$n == 150 & results$eqvar == 2 &results$dist == 1]

p4 <- ggplot()+
  aes(selection4)+
  geom_histogram()+
  xlab("Eq of var and dist violated")


grid.arrange(p1, p2, p3, p4, ncol = 2)

##### 
# BF

selection5 <- results$BF[results$n == 150 & results$eqvar == 1 &results$dist == 30]

p5 <- ggplot()+
  aes(selection5)+
  geom_histogram()+
  xlab("All assumptions met")+
  xlim(0, 1000)


selection6 <- results$BF[results$n == 150 & results$eqvar == 2 &results$dist == 30]

p6 <- ggplot()+
  aes(selection6)+
  geom_histogram()+
  xlab("Equality of variances violated")+
  xlim(0, 30)

selection7 <- results$BF[results$n == 150 & results$eqvar == 1 &results$dist == 1]

p7 <- ggplot()+
  aes(selection7)+
  geom_histogram()+
  xlab("Normal distribution violated")+
  xlim(0, 1)

selection8 <- results$BF[results$n == 150 & results$eqvar == 2 &results$dist == 1]

p8 <- ggplot()+
  aes(selection8)+
  geom_histogram()+
  xlab("Eq of var and dist violated")+
  xlim(0, 1)

grid.arrange(p5, p6, p7, p8, ncol = 2)


#####
# Visualize outliers with boxplots

bp1 <- ggplot()+
  aes(selection1)+
  geom_boxplot()+
  xlab("All assumptions met")

bp2 <- ggplot()+
  aes(selection2)+
  geom_boxplot()+
  xlab("Equality of variances violated")

bp3 <- ggplot()+
  aes(selection3)+
  geom_boxplot()+
  xlab("Normal distribution violated")

bp4 <- ggplot()+
  aes(selection4)+
  geom_boxplot()+
  xlab("Eq of var and dist violated")

grid.arrange(bp1, bp2, bp3, bp4, ncol = 2)

bp5 <- ggplot()+
  aes(selection5)+
  geom_boxplot()+
  xlab("All assumptions met")+
  scale_x_log10()

bp6 <- ggplot()+
  aes(selection6)+
  geom_boxplot()+
  xlab("Equality of variances violated")+
  scale_x_log10()

bp7 <- ggplot()+
  aes(selection7)+
  geom_boxplot()+
  xlab("Normal distribution violated")+
  scale_x_log10()

bp8 <- ggplot()+
  aes(selection8)+
  geom_boxplot()+
  xlab("Eq of var and dist violated")+
  scale_x_log10()

grid.arrange(bp5, bp6, bp7, bp8, ncol = 2)


set.seed(1)
par(mfrow=c(2,2))

test <- simulateData(300, 50, 1)
plot(test$predictor, test$outcome, 
     main = "all assumptions met")
test <- simulateData(300, 1, 1)
plot(test$predictor, test$outcome, 
     main = "normality violated")
test <- simulateData(300, 50, 2)
plot(test$predictor, test$outcome, 
     main = "equality of variances violated")
test <- simulateData(300, 1, 2)
plot(test$predictor, test$outcome, 
     main = "both violated")


test <- simulateData(150, 1, 2)

hist(test$outcome[test$predictor == 1])

ks.test(test$outcome[test$predictor == 1], "pnorm", mean=mean(test$outcome[test$predictor == 1]), sd=sd(test$outcome[test$predictor == 1]))

leveneTest(outcome ~ predictor, data = test)

hist(test$outcome[test$predictor == 3])

ano <- aov(data$outcome ~ data$predictor, data=data)
eta_sq(ano)


mean(results$BF[results$n == 150 & results$eqvar == 1 &results$dist == 1])
sd(results$BF[results$n == 150 & results$eqvar == 1 &results$dist == 1])





