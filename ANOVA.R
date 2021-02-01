rm(list = ls())
library(BayesFactor)
library(car)


nrepeats <- 10000
simulationOptions <- expand.grid(
  n            = c(30, 60, 150),                 # dividable by three and two
  distribution = c("normal", "bimodal"),
  homogeneity  = c(T, F),
  hypothesis   = c("0", "a"),
  repeats      = 1:nrepeats,
  stringsAsFactors = FALSE
)

nrow(simulationOptions)


sd      <- 1    # standard deviation of normal distributions making
                # the bimodal distribution --> also determines sd of 
                # normal distribution as its sd is the same as the 
                # one of the total bimodal distribution.
peaksbmd<- 8    # distance of peaks bimodal distribution
vardiff <- 2    # increase of error in case of homogeneity not met
b0      <- 0    # intercept
mean1   <- 0    # level predictor 1 --> also mean of all if h0 == T
mean2   <- 1.5    # level predictor 2
mean3   <- 3    # level predictor 3



simulateData <- function(n, distr, hom, hypothesis){
  
  sddet <- sd(c(rt(n/2, 10, 0 - peaksbmd/2), rt(n/2, 15, peaksbmd/2)))
  
  if (distr == "bimodal"){
    error <- sample(c(rt(n/2, 8, 0 - peaksbmd/2), rt(n/2, 8, peaksbmd/2)))
  }
  else{
    error <- rnorm(n, 0, sddet)
  }
  
  if (hom == T){
    epsilon1 <- error
    epsilon2 <- error
    epsilon3 <- error
  }
  
  else{
    epsilon1 <- error * 1/vardiff
    epsilon2 <- error 
    epsilon3 <- error * vardiff 
  }
    
  if (hypothesis == "0"){
    b1 <- mean1
    b2 <- mean1
    b3 <- mean1
  }
  else {
    b1 <- mean1
    b2 <- mean2
    b3 <- mean3
  }
    
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

plot(simulateData(60, "bimodal", T, "a"))

fitFreq      <- function(data){
  freq = aov(outcome ~ predictor, data=data)
  return(summary(freq)[[1]][1,5])
}

fitBayes     <- function(data){
  bf = anovaBF(outcome ~ predictor, data=data)
  return(extractBF(bf)$bf)
}


nsim <- nrow(simulationOptions)

resp <- rep(NA, length(seq_len(nsim)))
resb <- rep(NA, length(seq_len(nsim)))
etasq <- rep(NA, length(seq_len(nsim)))


for (i in seq_len(nsim)) {
  
  set.seed(i)
  opts <- simulationOptions[i, ]
  
  data <- simulateData(opts[[1]][1],opts[[2]][1], opts[[3]][1], opts[[4]][1])
  resultFreq  <- fitFreq(data)
  resultBayes <- fitBayes(data)
  etasquared <- eta_sq(aov(outcome ~ predictor, data=data))[[2]]
  
  resp[i] <- resultFreq
  resb[i] <- resultBayes
  etasq[i] <- etasquared
  
  
  print((100/(length(seq_len(nsim))))*i)
}


resultsAn <- data.frame(
  "n" = simulationOptions$n, 
  "distribution" = simulationOptions$distribution,
  "homogeneity" = simulationOptions$homogeneity,
  "hypothesis" = simulationOptions$hypothesis,
  "p" = resp, 
  "bf" = resb, 
  "etasq" = etasq)


# Plot p-values

hyp <- c("a", "0")
n <- c(30, 60, 150)
distr <- c("normal", "bimodal")
hom <- c(T, F)

m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)
layout(m,heights = c(0.3,0.3,0.3))

for (i in hyp){
  if (i == "a"){
    lab <- "Alternative hypothesis;"
    ysc <- c(0, 5)
  }
  else{
    lab <- "Null hypothesis;"
    ysc <- c(0, 1.3)
  }
  for (j in n){
      selectionAnP1 <- resultsAn$p[resultsAn$hypothesis== i & resultsAn$n == j & resultsAn$distribution == "normal" & resultsAn$homogeneity == T]
      
      selectionAnP2 <- resultsAn$p[resultsAn$hypothesis== i & resultsAn$n == j & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == T]
      
      selectionAnP3 <- resultsAn$p[resultsAn$hypothesis== i & resultsAn$n == j & resultsAn$distribution == "normal" & resultsAn$homogeneity == F]
      
      selectionAnP4 <- resultsAn$p[resultsAn$hypothesis== i & resultsAn$n == j & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == F]
      
      if (j == 150 & i == "a"){
        adj <- 10
      }
      else{
        adj <- 1
      }
      
      plot(density(selectionAnP1, adjust = adj), 
           lwd = 1, 
           col = 'black',
           xlab = "p value", 
           xlim = c(0, 1), 
           ylim = ysc,
           frame = F, 
           main = paste(lab, "n = ", toString(j)), 
           zero.line = F, 
           cex.lab = 1.4
      )
      lines(density(selectionAnP2, adjust = adj), lwd = 1, col = 'black', lty = 2)
      lines(density(selectionAnP3, adjust = adj), lwd = 1, col = 'black', lty = 3)
      lines(density(selectionAnP4, adjust = adj), lwd = 1, col = 'black', lty = 4)
      
      
      
    }
}
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",          
       legend = c("All assumptions met", 
                  "Normal error violated", 
                  "Homogeneity of variances violated", 
                  "All assumptions violated"),  
       lty = c(1, 2, 3, 4),          
       lwd = 1, 
       cex = 1.4)

# plot Bayes factors:

minor.ticks.axis <- function(ax,n,t.ratio=0.5,mn,mx,...){
  
  lims <- par("usr")
  if(ax %in%c(1,3)) lims <- lims[1:2] else lims[3:4]
  
  major.ticks <- pretty(lims,n=5)
  if(missing(mn)) mn <- min(major.ticks)
  if(missing(mx)) mx <- max(major.ticks)
  
  major.ticks <- major.ticks[major.ticks >= mn & major.ticks <= mx]
  
  labels <- sapply(major.ticks,function(i)
    as.expression(bquote(10^ .(i)))
  )
  axis(ax,at=major.ticks,labels=labels,...)
  
  n <- n+2
  minors <- log10(pretty(10^major.ticks[1:2],n))-major.ticks[1]
  minors <- minors[-c(1,n)]
  
  minor.ticks = c(outer(minors,major.ticks,`+`))
  minor.ticks <- minor.ticks[minor.ticks > mn & minor.ticks < mx]
  
  
  axis(ax,at=minor.ticks,tcl=par("tcl")*t.ratio,labels=FALSE)
}

hyp <- c("a", "0")
n <- c(30, 60, 150)
distr <- c("normal", "bimodal")
hom <- c(T, F)

m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)
layout(m,heights = c(0.4,0.4,0.3))

for (i in hyp){
  if (i == "a"){
    lab <- "Alternative hypothesis;"
    ysc <- c(0, 0.6)
    xsc <- c(-2, 8)
  }
  else{
    lab <- "Null hypothesis;"
    ysc <- c(0, 1.6)
    xsc <- c(-3, 1)
  }
  for (j in n){
    selectionAnB1 <- resultsAn$bf[resultsAn$hypothesis== i & resultsAn$n == j & resultsAn$distribution == "normal" & resultsAn$homogeneity == T & resultsAn$bf < 1000]
    
    selectionAnB2 <- resultsAn$bf[resultsAn$hypothesis== i & resultsAn$n == j & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == T & resultsAn$bf < 1000]
    
    selectionAnB3 <- resultsAn$bf[resultsAn$hypothesis== i & resultsAn$n == j & resultsAn$distribution == "normal" & resultsAn$homogeneity == F & resultsAn$bf < 1000]
    
    selectionAnB4 <- resultsAn$bf[resultsAn$hypothesis== i & resultsAn$n == j & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == F & resultsAn$bf < 1000]
    
    plot(density(log(selectionAnB1), adjust = 1), 
         lwd = 1, 
         col = 'black',
         xlab = "Bayes factor", 
         xlim = xsc, 
         ylim = ysc,
         frame = F, 
         main = paste(lab, "n = ", toString(j)), 
         zero.line = F, 
         cex.lab=1.4, 
         xaxt="n"
    )
    lines(density(log(selectionAnB2), adjust = 1), lwd = 1, col = 'black', lty = 2)
    lines(density(log(selectionAnB3), adjust = 1), lwd = 1, col = 'black', lty = 3)
    lines(density(log(selectionAnB4), adjust = 1), lwd = 1, col = 'black', lty = 4)
    minor.ticks.axis(1,9,mn=0,mx=8)  
  }
}
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",          
       legend = c("All assumptions met", 
                  "Normal error violated", 
                  "Homogeneity of variances violated", 
                  "All assumptions violated"),  
       lty = c(1, 2, 3, 4),          
       lwd = 1, 
       cex = 1.4
       )



#####

n <- 1000
mean <- rep(NA, n)
sdev <- rep(NA, n)
etasq <- rep(NA, n)
lev <- rep(NA, n)
shap <- rep(NA, n)


for (i in c(1:n)){
  set.seed(i)
  data <- simulateData(60, "bimodal", T, "a")
  sdev[i] <- sd(data$outcome)
  mean[i] <- mean(data$outcome)
  ano <- aov(outcome ~ predictor, data=data)
  etasq[i] <- eta_sq(ano)[[2]]
  lev[i] <- leveneTest(outcome ~ predictor, data = data)[[1,3]]
  shap[i] <- shapiro.test(data$outcome[data$predictor == 1])[[2]]
  
}

mean(mean)
mean(sdev)
mean(etasq)
mean(lev)
mean(shap)


# Tables

type1 <- rep(NA, 24)
etasq <- rep(NA, 24)
samplesize <- c(30, 60, 150)
c <- 0
hyp <- "0"


for (n in samplesize){
  
  # n p-values
  
  npAn1 <- c(length(resultsAn$p[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "normal" & resultsAn$homogeneity == T & resultsAn$p <= 0.05]), mean(resultsAn$etasq[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "normal" & resultsAn$homogeneity == T]))
  
  npAn2 <- c(length(resultsAn$p[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == T & resultsAn$p <= 0.05]), mean(resultsAn$etasq[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == T ]))
  
  npAn3 <- c(length(resultsAn$p[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "normal" & resultsAn$homogeneity == F & resultsAn$p <= 0.05]), mean(resultsAn$etasq[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "normal" & resultsAn$homogeneity == F]))
  
  npAn4 <- c(length(resultsAn$p[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == F & resultsAn$p <= 0.05]), mean(resultsAn$etasq[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == F]))
  
  # n bf ha
  
  nbAn1 <- c(length(resultsAn$bf[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "normal" & resultsAn$homogeneity == T & resultsAn$bf >= 3]), mean(resultsAn$etasq[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "normal" & resultsAn$homogeneity == T]))
  
  nbAn2 <- c(length(resultsAn$bf[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == T & resultsAn$bf >= 3]), mean(resultsAn$etasq[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == T]))
  
  nbAn3 <- c(length(resultsAn$bf[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "normal" & resultsAn$homogeneity == F & resultsAn$bf >= 3]), mean(resultsAn$etasq[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "normal" & resultsAn$homogeneity == F]))
  
  nbAn4 <- c(length(resultsAn$bf[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == F & resultsAn$bf >= 3]), mean(resultsAn$etasq[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == F]))
  
  selections <- list(nbAn1, npAn1, nbAn2, npAn2, nbAn3, npAn3, nbAn4, npAn4)
  
  for (i in selections){
    c <- c+1
    type1[c] <- (i[1])/nrepeats
    etasq[c] <- i[2]
    
  }
}


RtableAnT <- matrix(type1, ncol = 2, byrow = T)
colnames(RtableAnT) <- c("Type1_BF","Type1_p")
rownames(RtableAnT) <- c(rep("n = 30", 4), rep("n = 60", 4), rep("n = 150", 4))
RtableAnT


RtableAnetas <- matrix(etasq, ncol = 2, byrow = T)
colnames(RtableAnetas) <- c("mean_etasq_bf", "mean_etasq_p")
RtableAnetas

RtableAnTEs <- cbind(RtableAnT, RtableAnetas)
RtableAnTEs

# Order of violations per sample size category 
# "Normality",
# "Homogeneity of variances",
# "All violated",

par(mfrow = c(1, 1))
colors <- c("darkred", "goldenrod", "chartreuse4", "deepskyblue4")
size = 1.2

plot(data.frame(RtableAnTEs)$mean_etasq_bf[1:4], data.frame(RtableAnTEs)$Type1_BF[1:4], 
     ylim = c(0, 0.1), 
     xlim = c(0, 0.15), 
     pch = 21, 
     xlab = "Eta-squared", 
     ylab = "Type I error rate", 
     frame = F,
     col = colors, 
     cex = size, 
     cex.lab = 1
     )
points(data.frame(RtableAnTEs)$mean_etasq_bf[5:8], data.frame(RtableAnTEs)$Type1_BF[5:8], 
       pch = 24, 
       col = colors, 
       cex = size
       )
points(data.frame(RtableAnTEs)$mean_etasq_bf[9:12], data.frame(RtableAnTEs)$Type1_BF[9:12], 
       pch = 23, 
       col = colors, 
       cex = size
       )
points(data.frame(RtableAnTEs)$mean_etasq_p[1:4], data.frame(RtableAnTEs)$Type1_p[1:4], 
       pch = 21, 
       col = colors,
       bg = colors, 
       )
points(data.frame(RtableAnTEs)$mean_etasq_p[5:8], data.frame(RtableAnTEs)$Type1_p[5:8], 
       pch = 24,
       col = colors,
       bg = colors
       )
points(data.frame(RtableAnTEs)$mean_etasq_p[9:12], data.frame(RtableAnTEs)$Type1_p[9:12], 
       pch = 23, 
       col = colors,
       bg = colors
       )
legend(x = "topright",          
       legend = c("n = 30", 
                  "n = 60", 
                  "n = 150", 
                  "No violations",
                  "Normality",
                  "Homogeneity of variances",
                  "All violated",
                  "Frequentist", 
                  "Bayesian"),  
       pch = c(1, 2, 5, 15, 15, 15, 15, 16, 1),
       col = c("black", "black", "black", c(colors), "black", "black"),
       cex = 1
)

#####
power <- rep(NA, 24)
etasq <- rep(NA, 24)
samplesize <- c(30, 60, 150)
c <- 0
hyp <- "a"
alpha <- 0.03
bf <- 2

for (n in samplesize){
  
  # n p-values

    npAn1 <- c(length(resultsAn$p[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "normal" & resultsAn$homogeneity == T & resultsAn$p > alpha]), mean(resultsAn$etasq[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "normal" & resultsAn$homogeneity == T & resultsAn$p > alpha]))

  npAn2 <- c(length(resultsAn$p[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == T & resultsAn$p > alpha]), mean(resultsAn$etasq[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == T & resultsAn$p > alpha]))

  npAn3 <- c(length(resultsAn$p[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "normal" & resultsAn$homogeneity == F & resultsAn$p > alpha]), mean(resultsAn$etasq[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "normal" & resultsAn$homogeneity == F & resultsAn$p > alpha]))

  npAn4 <- c(length(resultsAn$p[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == F & resultsAn$p > alpha]), mean(resultsAn$etasq[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == F & resultsAn$p > alpha]))

  # # n bf ha

  nbAn1 <- c(length(resultsAn$bf[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "normal" & resultsAn$homogeneity == T & resultsAn$bf < bf]), mean(resultsAn$etasq[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "normal" & resultsAn$homogeneity == T & resultsAn$bf < bf]))

  nbAn2 <- c(length(resultsAn$bf[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == T & resultsAn$bf < bf]), mean(resultsAn$etasq[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == T & resultsAn$bf < bf]))

  nbAn3 <- c(length(resultsAn$bf[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "normal" & resultsAn$homogeneity == F & resultsAn$bf < bf]), mean(resultsAn$etasq[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "normal" & resultsAn$homogeneity == F & resultsAn$bf < bf]))

  nbAn4 <- c(length(resultsAn$bf[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == F & resultsAn$bf < bf]), mean(resultsAn$etasq[resultsAn$hypothesis== hyp & resultsAn$n == n & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == F & resultsAn$bf < bf]))

  
  selections <- list(nbAn1, npAn1, nbAn2, npAn2, nbAn3, npAn3, nbAn4, npAn4)
  
  for (i in selections){
    c <- c+1
    power[c] <- 1-((i[1])/nrepeats)
    etasq[c] <- i[2]
    
  }
}


RtableAnP <- matrix(power, ncol = 2, byrow = T)
colnames(RtableAnP) <- c("Power_BF","Power_p")
rownames(RtableAnP) <- c(rep("n = 30", 4), rep("n = 60", 4), rep("n = 150", 4))
RtableAnP


RtableAnEs<- matrix(etasq, ncol = 2, byrow = T)
colnames(RtableAnEs) <- c("mean_etasq_bf", "mean_etasq_p")

RtableAnPEs <- cbind(RtableAnP, RtableAnEs)
RtableAnPEs

# Order of violations per sample size category 
# "Normality",
# "Homogeneity of variances",
# "All violated",

# loss of power: 
#b
# n30
1- (0.0743  /0.1329  )

# n30
1- (0.1146    /0.2778    )

# n30
1- (0.3045    /0.6774    )

#f
# n30
1- (0.1201      /0.2222      )

# n30
1- (0.2028      /0.4292      )

# n30
1- (0.5046      /0.8434      )



par(mfrow = c(1,1))
colors <- c("darkred", "goldenrod", "chartreuse4", "deepskyblue4")
size <- 1.4

plot(data.frame(RtableAnPEs)$mean_etasq_bf[1:4], data.frame(RtableAnPEs)$Power_BF[1:4], 
     ylim = c(0, 1), 
     xlim = c(0, 0.15), 
     pch = 21, 
     xlab = "Eta-squared", 
     ylab = "Power", 
     frame = F, 
     col = colors,
     cex = size, 
     cex.lab = 1
)
points(data.frame(RtableAnPEs)$mean_etasq_bf[5:8], data.frame(RtableAnPEs)$Power_BF[5:8], 
       pch = 24, col = colors, cex = size)
points(data.frame(RtableAnPEs)$mean_etasq_bf[9:12], data.frame(RtableAnPEs)$Power_BF[9:12], 
       pch = 23, col = colors, cex = size)
points(data.frame(RtableAnPEs)$mean_etasq_p[1:4], data.frame(RtableAnPEs)$Power_p[1:4], 
       pch = 21, col = colors, bg = colors)
points(data.frame(RtableAnPEs)$mean_etasq_p[5:8], data.frame(RtableAnPEs)$Power_p[5:8], 
       pch = 24, col = colors, bg = colors)
points(data.frame(RtableAnPEs)$mean_etasq_p[9:12], data.frame(RtableAnPEs)$Power_p[9:12], 
       pch = 23, col = colors, bg = colors)
abline(lm(Power_BF[1:4] ~ mean_etasq_bf[1:4], data = data.frame(RtableAnPEs)), 
       lty = 3, col = rgb(0, 0, 0, 0.5))
abline(lm(Power_p[1:4] ~ mean_etasq_p[1:4], data = data.frame(RtableAnPEs)), 
       lty = 2, col = rgb(0, 0, 0, 0.5))
abline(lm(Power_BF[5:8] ~ mean_etasq_bf[5:8], data = data.frame(RtableAnPEs)), 
       lty = 3, col = rgb(0, 0, 0, 0.5))
abline(lm(Power_p[5:8] ~ mean_etasq_p[5:8], data = data.frame(RtableAnPEs)), 
       lty = 2, col = rgb(0, 0, 0, 0.5))
abline(lm(Power_BF[9:12] ~ mean_etasq_bf[9:12], data = data.frame(RtableAnPEs)), 
       lty = 3, col = rgb(0, 0, 0, 0.5))
abline(lm(Power_p[9:12] ~ mean_etasq_p[9:12], data = data.frame(RtableAnPEs)), 
       lty = 2, col = rgb(0, 0, 0, 0.5))
legend(x = "topright",          
       legend = c("n = 30", 
                  "n = 60", 
                  "n = 150", 
                  "No violations",
                  "Normality",
                  "Homogeneity of variances",
                  "All violated",
                  "Filled:  Frequentist", 
                  "Empty: Bayesian"),  
       pch = c(1, 2, 5, 15, 15, 15, 15, NA, NA),
       col = c("black", "black", "black", c(colors), "black", "black"),
       cex = 1
)

#####



# plot data 

par(mfrow = c(2, 2))

data <- simulateData(150, "normal", T, "a")
yl <- c(-20, 25)
yla <- "Outcome"

plot(simulateData(150, "normal", T, "a"), 
     frame = F, 
     main = "All assumptions met", 
     xlab = "Group",
     ylim = yl, 
     ylab = yla)
plot(simulateData(150, "bimodal", T, "a"), 
     frame = F, 
     main = "Normality violated", 
     xlab = "Group",
     ylim = yl, 
     ylab = yla)
plot(simulateData(150, "normal", F, "a"), 
     frame = F, 
     main = "Equality of variances violated", 
     xlab = "Group",
     ylim = yl, 
     ylab = yla)
plot(simulateData(150, "bimodal", F, "a"), 
     frame = F, 
     main = "All assumptions violated", 
     xlab = "Group",
     ylim = yl, 
     ylab = yla)



par(mfrow = c(1,1))

es1 <- (resultsAn$etasq[resultsAn$hypothesis== "a" & resultsAn$n == 150 & resultsAn$distribution == "normal" & resultsAn$homogeneity == T])
es2 <- (resultsAn$etasq[resultsAn$hypothesis== "a" & resultsAn$n == 150 & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == T])
es3 <- (resultsAn$etasq[resultsAn$hypothesis== "a" & resultsAn$n == 150 & resultsAn$distribution == "normal" & resultsAn$homogeneity == F])
es4 <- (resultsAn$etasq[resultsAn$hypothesis== "a" & resultsAn$n == 150 & resultsAn$distribution == "bimodal" & resultsAn$homogeneity == F])


conditions <- c("All met", "Normality", "Homogeneity ", "All violated")
barx <- barplot((c(mean(es1), mean(es2), mean(es3), mean(es4))), names.arg = conditions, ylab = "Eta-squared", ylim = c(0,0.1), border=NA)




