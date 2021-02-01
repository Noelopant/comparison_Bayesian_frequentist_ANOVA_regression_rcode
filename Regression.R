rm(list = ls())
library(BayesFactor)


nrepeats <- 10000
simulationOptions <- expand.grid(
  n            = c(15, 30, 60),
  linearity    = c(T, F),                    # switches between quadratic and linear relationship
  hsty         = c(T, F),                    # homoscedasticity
  normErr      = c(T, F),                    # normality of error variances
  independ     = c(T, F),                    # independence of errors --> 
  hypothesis   = c("0", "a"),
  repeats      = 1:nrepeats,
  stringsAsFactors = FALSE
)

nrow(simulationOptions)

b0 <- 10             # y-intercept
slope <- 0.2           # slope
dft  <- 2.5            # to violate normality of residuals a degree of freedom for rt -> outliers


simulateData <- function(n, lin, hsty, normerr, ind, hypoth){
  x <- runif(n, 0, 10)   # distribution of x values -> should be positive
  if (lin == T){
    l <- x
  }
  else{
    l <- 20*sin(0.3*x-0.2)
  }
  if (normerr == T)
    if (hsty == T){
      h <- rnorm(n)
    }
    else{
      h <- rnorm(n) * 0.6*sqrt(x)
    }
  else{
    if (hsty == T){
      h <- rt(n, dft)
    }
    else{
      h <- rt(n, dft) * sqrt(x)
    }
  }
  if (ind == T){
    i <- 0
  }
  else{
    i <- slope * 4 * sin(x*1.5)
  }
  if (hypoth == "0")
    b1 <- 0
  else{
    b1 <- slope
  }
  
  
  y <- b0 + (b1 * l) + (h + i)    # outcome
  
  
  df = data.frame(
    "predictor" = x, 
    "outcome" = y
  
  )
  return(df)
}


fitFreq <- function(data){
  freq = lm(data$outcome ~ data$predictor, data=data)
  return(summary(freq)$coefficients[2,4])
}


fitBayes     <- function(data){
  bf = regressionBF(outcome ~ predictor, data=data)
  return(extractBF(bf)$bf)
}


nsim <- nrow(simulationOptions)

resp <- rep(NA, length(seq_len(nsim)))
resb <- rep(NA, length(seq_len(nsim)))
rsq <- rep(NA, length(seq_len(nsim)))
size <- 1.2

for (i in seq_len(nsim)) {
  
  set.seed(i)
  opts <- simulationOptions[i, ]
  
  data <- simulateData(opts[[1]][1],opts[[2]][1], opts[[3]][1], opts[[4]][1], opts[[5]][1], opts[[6]][1])
  resultFreq  <- fitFreq(data)
  resultBayes <- fitBayes(data)
  rsquared <- summary(lm(outcome ~ predictor, data = data))$r.squared
  
  resp[i] <- resultFreq
  resb[i] <- resultBayes
  rsq[i] <- rsquared
  
  
  print((100/(length(seq_len(nsim))))*i)
}

# putting results into a dataframe for access

resultsReg <- data.frame(
  "n" = simulationOptions$n, 
  "linearity" = simulationOptions$linearity,
  "homoscedasticity" = simulationOptions$hsty,
  "normalityErrors" = simulationOptions$normErr,
  "independenceErr" = simulationOptions$independ,
  "hypothesis" = simulationOptions$hypothesis,
  "p" = resp, 
  "bf" = resb, 
  "rsq" = rsq)

#####
# Plots

hyp <- c("a", "0")
samplesize <- c(15, 30, 60)

m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)
layout(m,heights = c(0.3,0.3,0.3))

for (i in hyp){
  if (i == "a"){
    lab <- "Alternative hypothesis;"
    ysc <- c(0, 4)
  }
  else{
    lab <- "Null hypothesis;"
    ysc <- c(0, 1.1)
  }
  for (n in samplesize){
    selectionReP1 <- resultsReg$p[resultsReg$hypothesis== i & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T]
    
    selectionReP2 <- resultsReg$p[resultsReg$hypothesis== i & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T]
    
    selectionReP3 <- resultsReg$p[resultsReg$hypothesis== i & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == F & resultsReg$linearity == T & resultsReg$independenceErr == T]
    
    selectionReP4 <- resultsReg$p[resultsReg$hypothesis== i & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == F & resultsReg$linearity == F & resultsReg$independenceErr == F]
    
    if (n == 60 & i == "a"){
      adj <- 10
    }
    else{
      adj <- 1
    }
    
    plot(density(selectionReP1, adjust = adj), 
         lwd = 1, 
         col = 'black',
         xlim = c(0, 1), 
         ylim = ysc,
         frame = F, 
         main = paste(lab, "n = ", toString(n)), 
         xlab = "p value",
         zero.line = F, 
         cex.lab = 1.4
    )
    lines(density(selectionReP2, adjust = 2), lwd = 1, col = 'black', lty = 2)
    lines(density(selectionReP3, adjust = 2), lwd = 1, col = 'black', lty = 3)
    lines(density(selectionReP4, adjust = 2), lwd = 1, col = 'black', lty = 4)
    
    
    
  }
}
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",          
       legend = c("All assumptions met", 
                  "Normal error violated", 
                  "Homoscedasticity and normality violated", 
                  "All assumptions violated"),  
       lty = c(1, 2, 3, 4),          
       lwd = 1, 
       cex = 1.4)

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
samplesize <- c(15, 30, 60)

m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)
layout(m,heights = c(0.3,0.3,0.3))

for (i in hyp){
  if (i == "a"){
    lab <- "Alternative hypothesis;"
    ysc <- c(0, 0.6)
    xsc <- c(-2, 8)
  }
  else{
    lab <- "Null hypothesis;"
    ysc <- c(0, 2)
    xsc <- c(-2, 2)
  }
  for (n in samplesize){
   
     selectionReB1 <- resultsReg$bf[resultsReg$hypothesis== i & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$bf < 1000]
    
    selectionReB2 <- resultsReg$bf[resultsReg$hypothesis== i & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$bf < 1000]
    
    selectionReB3 <- resultsReg$bf[resultsReg$hypothesis== i & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == F & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$bf < 1000]
    
    selectionReB4 <- resultsReg$bf[resultsReg$hypothesis== i & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == F & resultsReg$linearity == F & resultsReg$independenceErr == F & resultsReg$bf < 1000]

        plot(density(log(selectionReB1), adjust = 1), 
         lwd = 1, 
         col = 'black',
         xlab = "Bayes factor", 
         xlim = xsc, 
         ylim = ysc,
         frame = F, 
         main = paste(lab, "n = ", toString(n)), 
         zero.line = F, 
         cex.lab = 1.4, 
         xaxt="n"
         )
    lines(density(log(selectionReB2), adjust = 1), lwd = 1, col = 'black', lty = 2)
    lines(density(log(selectionReB3), adjust = 1), lwd = 1, col = 'black', lty = 3)
    lines(density(log(selectionReB4), adjust = 1), lwd = 1, col = 'black', lty = 4)
    minor.ticks.axis(1,9,mn=0,mx=8)  }
}
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",          
       legend = c("All assumptions met", 
                  "Normal error violated", 
                  "Normal error and homoscedasticity violated", 
                  "All assumptions violated"),  
       lty = c(1, 2, 3, 4),          
       lwd = 1, 
       cex = 1.4)


# Tables

type1Re <- rep(NA, 36)
rsq <- rep(NA, 36)
samplesize <- c(15, 30, 60)
c <- 0
hyp <- "0"

for (n in samplesize){
  nbRe1 <- c(length(resultsReg$bf[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$bf >= 3]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T]))
  
  nbRe2 <- c(length(resultsReg$bf[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == F & resultsReg$bf >= 3]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == F ]))
  
  nbRe3 <- c(length(resultsReg$bf[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == F & resultsReg$independenceErr == T & resultsReg$bf >= 3]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == F & resultsReg$independenceErr == T]))
  
  nbRe4 <- c(length(resultsReg$bf[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == F & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$bf >= 3]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == F & resultsReg$linearity == T & resultsReg$independenceErr == T]))
  
  nbRe5 <- c(length(resultsReg$bf[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$bf >= 3]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T]))
  
  nbRe6 <- c(length(resultsReg$bf[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == F & resultsReg$linearity == F & resultsReg$independenceErr == F & resultsReg$bf >= 3]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == F & resultsReg$linearity == F & resultsReg$independenceErr == F]))
  
  npRe1 <- c(length(resultsReg$p[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$p <= 0.05]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T]))
  
  npRe2 <- c(length(resultsReg$p[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == F & resultsReg$p <= 0.05]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == F]))
  
  npRe3 <- c(length(resultsReg$p[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == F & resultsReg$independenceErr == T & resultsReg$p <= 0.05]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity ==T & resultsReg$linearity == F & resultsReg$independenceErr == T]))
  
  npRe4 <- c(length(resultsReg$p[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == F & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$p <= 0.05]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == F & resultsReg$linearity == T & resultsReg$independenceErr == T]))
  
  npRe5 <- c(length(resultsReg$p[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$p <= 0.05]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T]))
  
  npRe6 <- c(length(resultsReg$p[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == F & resultsReg$linearity == F & resultsReg$independenceErr == F & resultsReg$p <= 0.05]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == F & resultsReg$linearity == F & resultsReg$independenceErr == F]))
  
  selections <- list(nbRe1, npRe1, nbRe2, npRe2, nbRe3, npRe3, nbRe4, npRe4, nbRe5, npRe5, nbRe6, npRe6)
  
  for (i in selections){
    c <- c+1
    type1Re[c] <- i[1]/nrepeats
    rsq[c] <- i[2]
  }
}

RtableReT1 <- matrix(type1Re, ncol = 2, byrow = T)
colnames(RtableReT1) <- c("Type1_BF","Type1_p")
rownames(RtableReT1) <- c(rep("n = 15", 6), rep("n = 30", 6), rep("n = 60", 6))
RtableReT1

RtableReRs <- matrix(rsq, ncol = 2, byrow = T)
colnames(RtableReRs) <- c("mean_rsquared_bf", "mean_rsquared_p")

RtableAnReTRs <- cbind(RtableReT1, RtableReRs)
RtableAnReTRs

# Order of violations per sample size category
# "No violations",
# "Independence",
# "Linearity",
# "Homoscedasticity", 
# "Normality",
# "All violated",

par(mfrow = c(1,1))
colors <- c("darkred", "dimgrey", "goldenrod", "blueviolet", "chartreuse4", "deepskyblue4")
size = 1.2

plot(data.frame(RtableAnReTRs)$mean_rsquared_bf[1:6], data.frame(RtableAnReTRs)$Type1_BF[1:6], 
     ylim = c(0, 0.08), 
     xlim = c(0, 0.15), 
     xlab = "R-squared", 
     ylab = "Type 1", 
     frame = F, 
     pch = 24, 
     col = colors, 
     cex = size
     )
points(data.frame(RtableAnReTRs)$mean_rsquared_bf[7:12], data.frame(RtableAnReTRs)$Type1_BF[7:12], 
       pch = 24, 
       col = colors, 
       cex = size)
points(data.frame(RtableAnReTRs)$mean_rsquared_bf[13:18], data.frame(RtableAnReTRs)$Type1_BF[13:18], 
       pch = 23, 
       col = colors, 
       cex = size)
points(data.frame(RtableAnReTRs)$mean_rsquared_p[1:6], data.frame(RtableAnReTRs)$Type1_p[1:6], 
       pch = 21, 
       col = colors,
       bg = colors)
points(data.frame(RtableAnReTRs)$mean_rsquared_p[7:12], data.frame(RtableAnReTRs)$Type1_p[7:12], 
       pch = 24, 
       col = colors,
       bg = colors)
points(data.frame(RtableAnReTRs)$mean_rsquared_p[13:18], data.frame(RtableAnReTRs)$Type1_p[13:18], 
       pch = 23, 
       col = colors,
       bg = colors)
legend(x = "topright",          
       legend = c("n = 30", 
                  "n = 60", 
                  "n = 150", 
                  "No violations",
                  "Independence",
                  "Linearity",
                  "Homoscedasticity", 
                  "Normality",
                  "All violated",
                  "Frequentist", 
                  "Bayesian"),  
       pch = c(1, 2, 5, 15, 15, 15, 15, 15, 15, 16, 1),
       col = c("black", "black", "black", c(colors), "black", "black"),
       cex = 1
)

#####
powerRe <- rep(NA, 36)
rsq <- rep(NA, 36)
samplesize <- c(15, 30, 60)
c <- 0
hyp <- "a"
alpha <- 0.03
bf <- 2.25

for (n in samplesize){
  nbRe1 <- c(length(resultsReg$bf[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$bf < bf]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$bf < bf]))
  
  nbRe2 <- c(length(resultsReg$bf[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == F & resultsReg$bf < bf]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == F & resultsReg$bf < bf]))
  
  nbRe3 <- c(length(resultsReg$bf[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == F & resultsReg$independenceErr == T & resultsReg$bf < bf]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == F & resultsReg$independenceErr == T & resultsReg$bf < bf]))
  
  nbRe4 <- c(length(resultsReg$bf[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == F & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$bf < bf]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == F & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$bf < bf]))
  
  nbRe5 <- c(length(resultsReg$bf[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$bf < bf]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$bf < bf]))
  
  nbRe6 <- c(length(resultsReg$bf[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == F & resultsReg$linearity == F & resultsReg$independenceErr == F & resultsReg$bf < bf]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == F & resultsReg$linearity == F & resultsReg$independenceErr == F & resultsReg$bf < bf]))
  
  npRe1 <- c(length(resultsReg$p[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$p > alpha]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$p > alpha]))
  
  npRe2 <- c(length(resultsReg$p[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == F & resultsReg$p > alpha]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == F & resultsReg$p > alpha]))
  
  npRe3 <- c(length(resultsReg$p[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == F & resultsReg$independenceErr == T & resultsReg$p > alpha]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity ==T & resultsReg$linearity == F & resultsReg$independenceErr == T & resultsReg$p > alpha]))
  
  npRe4 <- c(length(resultsReg$p[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == F & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$p > alpha]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == F & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$p > alpha]))
  
  npRe5 <- c(length(resultsReg$p[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$p > alpha]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T & resultsReg$p > alpha]))
  
  npRe6 <- c(length(resultsReg$p[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == F & resultsReg$linearity == F & resultsReg$independenceErr == F & resultsReg$p > alpha]), mean(resultsReg$rsq[resultsReg$hypothesis== hyp & resultsReg$n == n & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == F & resultsReg$linearity == F & resultsReg$independenceErr == F & resultsReg$p > alpha]))
  
  selections <- list(nbRe1, npRe1, nbRe2, npRe2, nbRe3, npRe3, nbRe4, npRe4, nbRe5, npRe5, nbRe6, npRe6)
  
  for (i in selections){
    c <- c+1
    powerRe[c] <- 1-(i[1]/nrepeats)
    rsq[c] <- i[2]
  }
}

RtableReP <- matrix(powerRe, ncol = 2, byrow = T)
colnames(RtableReP) <- c("Power_BF","Power_p")
rownames(RtableReP) <- c(rep("n = 15", 6), rep("n = 30", 6), rep("n = 60", 6))
RtableReP

RtableReRs <- matrix(rsq, ncol = 2, byrow = T)
colnames(RtableReRs) <- c("mean_rsquared_bf", "mean_rsquared_p")

RtableRePRs <- cbind(RtableReP, RtableReRs)
RtableRePRs

# Order of violations per sample size category
# "No violations",
# "Independence",
# "Linearity",
# "Homoscedasticity", 
# "Normality",
# "All violated",

# loss of power: 
# bayesian
# 15
1- (0.0728 / 0.3631)

# 30
1- (0.1023/ 0.7345)

# 60
1- (0.1684 / 0.9727)

# frequentist
# 15
1- (0.1248 / 0.5146)

# 30
1- (0.1788/ 0.8412)

# 60
1- (0.2864/ 0.9877)



par(mfrow = c(1,1))
colors <- c("darkred", "dimgrey", "goldenrod", "blueviolet", "chartreuse4", "deepskyblue4")
size = 1.2

plot(data.frame(RtableRePRs)$mean_rsquared_bf[1:6], data.frame(RtableRePRs)$Power_BF[1:6], 
     ylim = c(0, 1), 
     xlim = c(0, 0.25), 
     xlab = "R-squared", 
     ylab = "Power", 
     frame = F, 
     pch = 21, 
     col = colors, 
     cex = size)
points(data.frame(RtableRePRs)$mean_rsquared_bf[7:12], data.frame(RtableRePRs)$Power_BF[7:12], 
       pch = 24, 
       col = colors, 
       cex = size)
points(data.frame(RtableRePRs)$mean_rsquared_bf[13:18], data.frame(RtableRePRs)$Power_BF[13:18], 
       pch = 23, 
       col = colors, 
       cex = size)
points(data.frame(RtableRePRs)$mean_rsquared_p[1:6], data.frame(RtableRePRs)$Power_p[1:6], 
       pch = 21, 
       col = colors,
       bg = colors)
points(data.frame(RtableRePRs)$mean_rsquared_p[7:12], data.frame(RtableRePRs)$Power_p[7:12], 
       pch = 24, 
       col = colors,
       bg = colors)
points(data.frame(RtableRePRs)$mean_rsquared_p[13:18], data.frame(RtableRePRs)$Power_p[13:18], 
       pch = 23, 
       col = colors,
       bg = colors)
abline(lm(Power_BF[1:6] ~ mean_rsquared_bf[1:6], data = data.frame(RtableRePRs)), 
       lty = 3, col = rgb(0, 0, 0, 0.5))
abline(lm(Power_BF[7:12] ~ mean_rsquared_bf[7:12], data = data.frame(RtableRePRs)), 
       lty = 2, col = rgb(0, 0, 0, 0.5))
abline(lm(Power_BF[13:18] ~ mean_rsquared_bf[13:18], data = data.frame(RtableRePRs)), 
       lty = 3, col = rgb(0, 0, 0, 0.5))
abline(lm(Power_p[1:6] ~ mean_rsquared_p[1:6], data = data.frame(RtableRePRs)), 
       lty = 2, col = rgb(0, 0, 0, 0.5))
abline(lm(Power_p[7:12] ~ mean_rsquared_p[7:12], data = data.frame(RtableRePRs)), 
       lty = 3, col = rgb(0, 0, 0, 0.5))
abline(lm(Power_p[13:18] ~ mean_rsquared_p[13:18], data = data.frame(RtableRePRs)), 
       lty = 2, col = rgb(0, 0, 0, 0.5))
legend(x = "topright",          
       legend = c("n = 15", 
                  "n = 30", 
                  "n = 60", 
                  "No violations",
                  "Independence",
                  "Linearity",
                  "Homoscedasticity", 
                  "Normality",
                  "All violated",
                  "Filled: Frequentist", 
                  "Empty: Bayesian"),  
       pch = c(1, 2, 5, 15, 15, 15, 15, 15, 15, NA, NA),
       col = c("black", "black", "black", c(colors), "black", "black"),
       cex = 1
)





#####
  
head(opts)

n <- 1000
r2 <- rep(NA, n)
mean <- rep(NA, n)
sdev <- rep(NA, n)
shap <- rep(NA, n)

for (i in c(1:n)){
  set.seed(i)
  data <- simulateData(60, T, T, F, T, "a")
  r2[i] <- summary(lm(outcome ~ predictor, data = data))$r.squared
  mean <- mean(data$outcome)
  sdev <- sd(data$outcome)
  shap[i] <- shapiro.test(data$outcome)[[2]]
}

mean(mean)
mean(sdev)
mean(r2)
mean(shap)



par(mfrow = c(2, 3))

data <-simulateData(60, T, T, T, T, "a")
size <- 1.4

plot(data, ylim = c(6, 20), 
     main = "All assumptions met",
     xlab = "predictor",
     frame = F, 
     cex.lab = size, 
     ylab = "Outcome"
)
abline(lm(outcome ~ predictor, data = data), col = "red")
plot(simulateData(60, F, T, T, T, "a"), ylim = c(6, 20), main = "Linearity violated", 
     frame = F, cex.lab = size, ylab = "Outcome",xlab = "predictor"
)
curve((4*sin(0.3*x-0.2)+10), from = -1, to = 11, add = T, col = "red")
plot(simulateData(60, T, F, T, T, "a"), ylim = c(6, 20), main = "Homoscedasticity violated", 
     frame = F, cex.lab = size, xlab = "predictor", ylab = "Outcome"
)
plot(simulateData(60, T, T, F, T, "a"), ylim = c(6, 20), main = "Normality violated", 
     frame = F, cex.lab = size, xlab = "predictor", ylab = "Outcome"
)
plot(simulateData(60, T, T, T, F, "a"), ylim = c(6, 20), main = "Independence violatd", 
     frame = F, cex.lab = size, xlab = "predictor", ylab = "Outcome"
)
curve((0.2 * x + 0.2 * 4* sin(x*1.5)+10), from = -1, to = 11, add = T, col = "red")
plot(simulateData(100, F, F, F, F, "a"), ylim = c(6, 20), main = "All assumptions violated",
     frame = F, cex.lab = size, xlab = "predictor", ylab = "Outcome"
)


par(mfrow = c(1,1))

es1 <- mean(resultsReg$rsq[resultsReg$hypothesis== "a" & resultsReg$n == 60 & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T])
es2 <- mean(resultsReg$rsq[resultsReg$hypothesis== "a" & resultsReg$n == 60 & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == F & resultsReg$independenceErr == T])
es3 <- mean(resultsReg$rsq[resultsReg$hypothesis== "a" & resultsReg$n == 60 & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == F])
es4 <- mean(resultsReg$rsq[resultsReg$hypothesis== "a" & resultsReg$n == 60 & resultsReg$normalityErrors == T & resultsReg$homoscedasticity == F & resultsReg$linearity == T & resultsReg$independenceErr == T])
es5 <- mean(resultsReg$rsq[resultsReg$hypothesis== "a" & resultsReg$n == 60 & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == T & resultsReg$linearity == T & resultsReg$independenceErr == T])
es6 <- mean(resultsReg$rsq[resultsReg$hypothesis== "a" & resultsReg$n == 60 & resultsReg$normalityErrors == F & resultsReg$homoscedasticity == F & resultsReg$linearity == F & resultsReg$independenceErr == F])

conditions <- c("All met", "Linearity ", "Independence ", "Homoscedasticity", "Normality", "All violated")
barplot(c(es1, es2, es3, es4, es5, es6), names.arg = conditions, ylab = "R-squared", ylim = c(0,0.3), cex.names = 0.9, cex.axis = 0.9, cex.lab = 0.9)


write.csv(resultsReg, "thesis\\resultsreg.csv")
