rm(list = ls())
library(ggplot2)
library(sn)
library(BayesFactor)
library(stats)
library(sjstats)
library(car)
library(gridExtra)


visualize <- function(y, x = NULL) {
  
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    layout(1)
  })
  
  par(las = 1, bty = "n", cex = 2)
  if (is.null(x))
    layout(matrix(1:2, ncol = 2))
  else {
    layout(matrix(1:3, ncol = 3))
    plot(x, y)
    lines(x, predict(lm(y ~ x)), col = 2)
  }
  
  boxplot(y)
  densityEstimate <- density(y)
  ex <- range(densityEstimate$x)
  xx <- seq(ex[1], ex[2], length.out = 1e4)
  normalFit <- dnorm(xx, mean(y), sd(y))
  hist(y, probability = TRUE, main = "", ylim = range(c(normalFit, densityEstimate$y)))
  lines(densityEstimate)
  lines(xx, normalFit, col = 2)
  legend("topright", legend = c("density estimate", "normal fit"), col = 1:2, bty = "n", lty = 1)
  return(invisible(NULL))
}

#####


nrepeats <- 100
simulationOptions <- expand.grid(
  n            = c(150, 300, 600),
  linearity = c(1, 0),
  errorvariance = c(1, 0),
  repeats      = 1:nrepeats,
  stringsAsFactors = FALSE
)

simulationOptions

sigma   <- 1                                       # standard deviation
mean    <- 0
x       <- rnorm(n, 0, sigma)
epsilon <- rnorm(n, mean, sigma)
b0      <- 10                                       # intercept
b1      <- 0.2

simulateData <- function(n, lin, evar){
  x       <- rnorm(n, 0, sigma)
  epsilon <- rnorm(n, mean, sigma)
  
  if (lin == 1){
    nonlin <- 1
  }
  else{
    nonlin <- x
  }
  if (evar == 1){
    evariance <- 1
  }
  else{
    evariance <- x
  }
  
  y <- b0 + (b1 * x * nonlin) + (epsilon * evariance)    # outcome
  

  df = data.frame(
    "predictor" = x, 
    "outcome" = y
  )
  return(df)
}




fitFreq      <- function(data){
  freq = lm(data$outcome ~ data$predictor, data=data)
  return(summary(freq)$r.squared)
}


fitBayes     <- function(data){
  bf = regressionBF(outcome ~ predictor, data=data)
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
  "linearity" = simulationOptions$linearity,
  "errorvariances" = simulationOptions$errorvariance,
  "r2" = resf, 
  "BF" = resb)



#####
# Make histograms of r2 and bayes factors

selection1 <- results$r2[results$n == 150 & results$linearity == 1 &results$errorvariances == 1]

p1 <- ggplot()+
  aes(selection1)+
  geom_histogram()+
  xlab("All assumptions met")

selection2 <- results$r2[results$n == 150 & results$linearity == 0 &results$errorvariances == 1]

p2 <- ggplot()+
  aes(selection2)+
  geom_histogram()+
  xlab("Linearity violated")


selection3 <- results$r2[results$n == 150 & results$linearity == 1 &results$errorvariances == 0]

p3 <- ggplot()+
  aes(selection3)+
  geom_histogram()+
  xlab("Errorvariance violated")

selection4 <- results$r2[results$n == 150 & results$linearity == 0 &results$errorvariances == 0]

p4 <- ggplot()+
  aes(selection4)+
  geom_histogram()+
  xlab("Errorvariance & linearity violated")


grid.arrange(p1, p2, p3, p4, ncol = 2)

##### 
# BF

selection5 <- results$BF[results$n == 150 & results$linearity == 1 &results$errorvariances == 1]

p5 <- ggplot()+
  aes(selection5)+
  geom_histogram()+
  xlab("All assumptions met")+
  xlim(0,500)

selection6 <- results$BF[results$n == 150 & results$linearity == 0 &results$errorvariances == 1]

p6 <- ggplot()+
  aes(selection6)+
  geom_histogram()+
  xlab("Linearity violated")+
  xlim(0, 1)

selection7 <- results$BF[results$n == 150 & results$linearity == 1 &results$errorvariances == 0]

p7 <- ggplot()+
  aes(selection7)+
  geom_histogram()+
  xlab("Errorvariance violated")+
  xlim(1,500)

selection8 <- results$BF[results$n == 150 & results$linearity == 0 &results$errorvariances == 0]

p8 <- ggplot()+
  aes(selection8)+
  geom_histogram()+
  xlab("Errorvariance & linearity violated")+
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
  xlab("Linearity violated")

bp3 <- ggplot()+
  aes(selection3)+
  geom_boxplot()+
  xlab("Errorvariance violated")

bp4 <- ggplot()+
  aes(selection4)+
  geom_boxplot()+
  xlab("Errorvariance & linearity violated")

grid.arrange(bp1, bp2, bp3, bp4, ncol = 2)

bp5 <- ggplot()+
  aes(selection5)+
  geom_boxplot()+
  xlab("All assumptions met")+
  scale_x_log10()

bp6 <- ggplot()+
  aes(selection6)+
  geom_boxplot()+
  xlab("Linearity violated")+
  scale_x_log10()

bp7 <- ggplot()+
  aes(selection7)+
  geom_boxplot()+
  xlab("Errorvariance violated")+
  scale_x_log10()

bp8 <- ggplot()+
  aes(selection8)+
  geom_boxplot()+
  xlab("Errorvariance & linearity violated")+
  scale_x_log10()

grid.arrange(bp5, bp6, bp7, bp8, ncol = 2)



b1 <- 0.8
set.seed(1)

par(mfrow=c(2,2))

test <- simulateData(400, 1, 1)
plot(test$predictor, test$outcome, 
     main = "all assumptions met", 
     ylab = "", 
     xlab = "")
test <- simulateData(400, 0, 1)
plot(test$predictor, test$outcome, 
     main = "linearity violated", 
     ylab = "", 
     xlab = "")
test <- simulateData(400, 1, 0)
plot(test$predictor, test$outcome, 
     main = "homoscedacity violated", 
     ylab = "", 
     xlab = "")
test <- simulateData(400, 0, 0)
plot(test$predictor, test$outcome, 
     main = "both violated", 
     ylab = "", 
     xlab = "")


par(mfrow=c(1,1))

test <- simulateData(150, 2, 2)
hist(test$outcome)

ks.test(test$outcome, "pnorm", mean=mean(test$outcome), sd=sd(test$outcome))

visualize(test$outcome, test$predictor)

plot(test$predictor, test$outcome)



citation("stats") 


mean(results$BF[results$n == 150 & results$linearity == 0 &results$errorvariances == 0])
sd(results$BF[results$n == 150 & results$linearity == 0 &results$errorvariances == 0])




