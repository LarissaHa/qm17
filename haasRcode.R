###########################################################
#
# title: QM - Data Essay
# author: Larissa Haas, 1417669
# date: December 7, 2017
#
###########################################################
# SETTING UP THE WORKING SPACE

rm(list = ls())
setwd("C:/...") # be sure to set the right working directory
library(MASS)
library(foreign)
library(ggplot2)
library(hydroGOF)
library(stargazer)

###########################################################
# LOADING THE DATA AND GETTING A FEELING FOR IT

dta <- load("Dataessay.Rdata")
#summary(GERdata)
#summary(UKdata)

# RE-CODING SOME OF THE VARIABLES
GERdata$FDP <- ifelse(GERdata$party_affiliation == "FDP", 1, 0)
GERdata$CDU_CSU <- ifelse(GERdata$party_affiliation == "CDU/CSU", 1, 0)
GERdata$SPD <- ifelse(GERdata$party_affiliation == "SPD", 1, 0)
GERdata$Linke <- ifelse(GERdata$party_affiliation == "DIE LINKE.", 1, 0)

GERdata$ingov <- GERdata$SPD + GERdata$CDU_CSU

GERdata$list_candidate_new <- ifelse(GERdata$list_candidate == "District", 0, 1)

#stargazer(GERdata)
#stargazer(UKdata)

hist(GERdata$number_speeches, 
     main = "Frequency of Number of Speeches (Count Data) \n 
     in Germany", 
     xlab = "Number of Speeches (counted per MP)")
hist(UKdata$number_speeches, 
     main = "Frequency of Number of Speeches (Count Data) \n 
     in the United Kingdom",
     xlab = "Number of Speeches (counted per MP)")
#hist(GERdata$ideological_distance)
#hist(UKdata$ideological_distance)

###########################################################
# BASIC MODELS AND LIKELIHOOD-RATIO-TEST

# BASIC POISSON MODEL
m0a <- glm(number_speeches ~ ideological_distance, 
           data = GERdata, family = "poisson")
m0b <- glm(number_speeches ~ ideological_distance, 
           data = UKdata, family = "poisson")

# BASIC NEGATIVE BINOMIAL MODEL
m1a <- glm.nb(number_speeches ~ ideological_distance, 
              data=GERdata, control = glm.control(maxit=100))
m1b <- glm.nb(number_speeches ~ ideological_distance, 
              data=UKdata, control = glm.control(maxit=100))

# FUNCTION: LIKELIHOOD-RATIO-TEST
compute.lrt <- function(model.nb, model.po, sign.level){
  L1 <- logLik(model.po) 
  L2 <- logLik(model.nb) 
  LRT <- -2*L1 + 2*L2
  return(LRT > qchisq(sign.level, df = 1))
}

compute.lrt(m1a, m0a, 0.95) # TRUE --> negative binomial model is better!

# EXPANDED MODELS
m2a <- glm.nb(number_speeches ~ party_leader + 
              ideological_distance,
              data=GERdata, 
              control = glm.control(maxit=100))
m2b <- glm.nb(number_speeches ~ party_leader + 
              ideological_distance, 
              data=UKdata, 
              control = glm.control(maxit=100))

#stargazer(m1a, m1b, m2a, m2b)

###########################################################
# BASIC MODEL TESTING

anova(m1a, m2a)
anova(m1b, m2b)

newGER <- GERdata
newGER$m1a <- predict(m1a, newGER, type = "response")
newGER$m2a <- predict(m2a, newGER, type = "response")

rmse(newGER$m1a, GERdata$number_speeches)
rmse(newGER$m2a, GERdata$number_speeches)

newUK <- UKdata
newUK$m1b <- predict(m1b, newUK, type = "response")
newUK$m2b <- predict(m2b, newUK, type = "response")

rmse(newUK$m1b, UKdata$number_speeches)
rmse(newUK$m2b, UKdata$number_speeches)

###########################################################
# SIMULATING THE DIFFERENCE BETWEEN THE SYSTEMS

gamma.hat.uk1 <- coef(m2b)
V.hat.uk1 <- vcov(m2b)
S.uk1 <- mvrnorm(1000, gamma.hat.uk1, V.hat.uk1)

gamma.hat.ger1 <- coef(m2a)
V.hat.ger1 <- vcov(m2a)
S.ger1 <- mvrnorm(1000, gamma.hat.ger1, V.hat.ger1)

# SCENARIOS: HIGHEST AND LOWEST VALUE OF IDEOLOGICAL DISTANCE,
# SIMULATE THE VALUES FOR EACH COUNTRY SEPARATELY
germin.ideo.scen <- cbind(1, 0, min(GERdata$ideological_distance, na.rm=TRUE))
germax.ideo.scen <- cbind(1, 0, max(GERdata$ideological_distance, na.rm=TRUE))
ukmin.ideo.scen <- cbind(1, 0, min(UKdata$ideological_distance, na.rm=TRUE))
ukmax.ideo.scen <- cbind(1, 0, max(UKdata$ideological_distance, na.rm=TRUE))

Xbeta.germin <- S.ger1 %*% t(germin.ideo.scen)
Xbeta.germax <- S.ger1 %*% t(germax.ideo.scen)
Xbeta.ukmin <- S.uk1 %*% t(ukmin.ideo.scen)
Xbeta.ukmax <- S.uk1 %*% t(ukmax.ideo.scen)

lambda.germin <- exp(Xbeta.germin)
lambda.germax <- exp(Xbeta.germax)
lambda.ukmin <- exp(Xbeta.ukmin)
lambda.ukmax <- exp(Xbeta.ukmax)

theta.ger <- m2a$theta
theta.uk <- m2b$theta

exp.germin <- sapply(lambda.germin, 
                    function(x) mean(rnbinom(100, size = theta.ger, mu = x)))
exp.germax <- sapply(lambda.germax, 
                    function(x) mean(rnbinom(100, size = theta.ger, mu = x)))
exp.ukmin <- sapply(lambda.ukmin, 
                    function(x) mean(rnbinom(100, size = theta.uk, mu = x)))
exp.ukmax <- sapply(lambda.ukmax, 
                    function(x) mean(rnbinom(100, size = theta.uk, mu = x)))

exp.ger <- c(exp.germin, exp.germax)
exp.uk <- c(exp.ukmin, exp.ukmax)
df.ger <- data.frame(exp.ger)
df.uk <- data.frame(exp.uk)

df.ger$id <- c(rep("min", 1000), rep("max", 1000))
df.uk$id <- c(rep("min", 1000), rep("max", 1000))

ggplot(df.ger, aes(x = exp.ger, fill = id)) + 
  geom_density(alpha = 0.4) + 
  guides(fill = guide_legend(title = "")) +
  xlab("Expected Speeches") +
  ylab("Density") + 
  ggtitle("Simulating German Number of Speeches") + 
  theme_bw()

ggplot(df.uk, aes(x = exp.uk, fill = id)) + 
  geom_density(alpha = 0.4) + 
  guides(fill = guide_legend(title = "")) +
  xlab("Expected Speeches") +
  ylab("Density") + 
  ggtitle("Simulating British Number of Speeches") + 
  theme_bw()

###########################################################
# MODELLING THE UK CASE

# MODEL 3: FULL MODEL
m3b <- glm.nb(number_speeches ~ party_leader + 
              ideological_distance + 
              conservative_MP, 
              data=UKdata, 
              control = glm.control(maxit=100))

# MODEL 4: REDUCED (BEST) MODEL
m4b <- glm.nb(number_speeches ~ ideological_distance + 
              conservative_MP, 
              data=UKdata, 
              control = glm.control(maxit=100))

anova(m4b, m3b)

newUK$m3b <- predict(m3b, newUK, type = "response")
rmse(newUK$m3b, UKdata$number_speeches)

newUK$m4b <- predict(m4b, newUK, type = "response")
rmse(newUK$m4b, UKdata$number_speeches)

#stargazer(m3b, m4b)

###########################################################
# MODELLING THE GERMAN CASE

# MODEL 3: EXPANDED MODEL
m3a <- glm.nb(number_speeches ~ party_leader + 
            ideological_distance + 
            list_candidate_new + 
            committee + 
            caolMPoutside + 
            ingov, 
            data = GERdata, 
            control = glm.control(maxit = 100))

anova(m2a, m3a)

newGER$m3a <- predict(m3a, newGER, type = "response")
rmse(newGER$m3a, GERdata$number_speeches)

# MODEL 4: FULL MODEL
m4a <- glm.nb(number_speeches ~ party_leader + 
              ideological_distance + 
              list_candidate_new + 
              committee + 
              caolMPoutside + 
              ingov + 
              FDP + 
              SPD + 
              Linke, 
              data = GERdata, 
              control = glm.control(maxit = 100))

anova(m3a, m4a)

newGER$m4a <- predict(m4a, newGER, type = "response")
rmse(newGER$m4a, GERdata$number_speeches)

# MODEL 5: BEST MODEL
m5a <- glm.nb(number_speeches ~ ideological_distance + 
              committee + 
              ingov, 
              data = GERdata, 
              control = glm.control(maxit = 100))

anova(m5a, m4a)

newGER$m5a <- predict(m5a, newGER, type = "response")
rmse(newGER$m5a, GERdata$number_speeches)

# MODEL 6: PARTY MODEL (FOR SIMULATIONS)
m6a <- glm.nb(number_speeches ~ ideological_distance + 
              committee + 
              FDP + 
              CDU_CSU + 
              Linke + 
              SPD, 
              data = GERdata,
              control = glm.control(maxit = 100))

newGER$m6a <- predict(m6a, newGER, type = "response")
rmse(newGER$m6a, GERdata$number_speeches)

#stargazer(m3a, m4a, m5a, m6a)

###########################################################
# SIMULATION FOR FIRST DIFFERENCES

quants.mean.fun <- function(x){
  c(quants = quantile(x, probs=c(0.025,0.5,0.975), mean = mean(x)))
}

gamma.hat.ger2 <- coef(m3a)
V.hat.ger2 <- vcov(m3a)
S.ger2 <- mvrnorm(1000, gamma.hat.ger2, V.hat.ger2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# scenario: (y/n) party leader, median distance, median list, 
# media committees, median outside, median ingov
scenario.lead <- cbind(1, 1, 0.1253, 1, 3, 0, 1) 
scenario.notl <- cbind(1, 0, 0.1253, 1, 3, 0, 1) 

X.lead <- as.matrix(rbind(scenario.lead, scenario.notl))
lead.combined <- S.ger2 %*% t(X.lead)
fd.lead <- lead.combined[,1] - lead.combined[,2]
quants.fd.lead <- apply(as.matrix(fd.lead), 2, quants.mean.fun)
quants.fd.lead

hist(fd.lead, main = "First Differences between Partyleaders 
     and Non-Partyleaders" )
abline(v = quants.fd.lead, lty = 2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# scenario: no party leader, median distance, (y/n) list candidate, 
# media committees, median outside, median ingov
scenario.list <- cbind(1, 0, 0.1253, 1, 3, 0, 1) 
scenario.dire <- cbind(1, 0, 0.1253, 0, 3, 0, 1) 

X.list <- as.matrix(rbind(scenario.list, scenario.dire))
list.combined <- S.ger2 %*% t(X.list)
fd.list <- list.combined[,1] - list.combined[,2]
quants.fd.list <- apply(as.matrix(fd.list), 2, quants.mean.fun)
quants.fd.list

hist(fd.list, main = "First Differences between List Candidates and 
     Direct Mandates" )
abline(v = quants.fd.list, lty = 2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# scenario: no party leader, median distance, median list, 
# media committees, (y/n) outside, median inreg
scenario.inco <- cbind(1, 0, 0.1253, 1, 3, 0, 1) 
scenario.outc <- cbind(1, 0, 0.1253, 1, 3, 1, 1) 

X.inco <- as.matrix(rbind(scenario.inco, scenario.outc))
inco.combined <- S.ger2 %*% t(X.inco)
fd.inco <- inco.combined[,1] - inco.combined[,2]
quants.fd.inco <- apply(as.matrix(fd.inco), 2, quants.mean.fun)
quants.fd.inco

hist(fd.inco, main = "First Differences between Speekers \n within 
     the Coalition Interval and outside" )
abline(v = quants.fd.inco, lty = 2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# scenario: no party leader, median distance, median list, 
# media committees, median outside, (y/n) inreg
scenario.ingo <- cbind(1, 0, 0.1253, 1, 3, 0, 1) 
scenario.outg <- cbind(1, 0, 0.1253, 1, 3, 0, 0) 

X.ingo <- as.matrix(rbind(scenario.ingo, scenario.outg))
ingo.combined <- S.ger2 %*% t(X.ingo)
fd.ingo <- ingo.combined[,1] - ingo.combined[,2]
quants.fd.ingo <- apply(as.matrix(fd.ingo), 2, quants.mean.fun)
quants.fd.ingo

hist(fd.ingo, main = "First Differences between Speekers \n from 
     the Governing Parties and the Opposition" )
abline(v = quants.fd.ingo, lty = 2)

###########################################################
# SIMULATING OTHER QUANTITIES OF INTEREST

gamma.hat.ger3 <- coef(m6a)
V.hat.ger3 <- vcov(m6a)
S.ger3 <- mvrnorm(1000, gamma.hat.ger3, V.hat.ger3)

# SCENARIOS: TAKE MEDIAN EXCEPT FOR PARTY VALUES
scenario.FDP <- cbind(1, 0.1253, 3, 1, 0, 0, 0) 
scenario.CDU <- cbind(1, 0.1253, 3, 0, 1, 0, 0) 
scenario.LIN <- cbind(1, 0.1253, 3, 0, 0, 1, 0) 
scenario.SPD <- cbind(1, 0.1253, 3, 0, 0, 0, 1)
scenario.GRU <- cbind(1, 0.1253, 3, 0, 0, 0, 0) 

XbetaFDP <- S.ger3 %*% t(scenario.FDP)
XbetaCDU <- S.ger3 %*% t(scenario.CDU)
XbetaLIN <- S.ger3 %*% t(scenario.LIN)
XbetaSPD <- S.ger3 %*% t(scenario.SPD)
XbetaGRU <- S.ger3 %*% t(scenario.GRU)

lambdaFDP <- exp(XbetaFDP)
lambdaCDU <- exp(XbetaCDU)
lambdaLIN <- exp(XbetaLIN)
lambdaSPD <- exp(XbetaSPD)
lambdaGRU <- exp(XbetaGRU)

thetam6a <- m6a$theta

exp.FDP <- sapply(lambdaFDP, 
                  function(x) mean(rnbinom(1000, size = thetam6a, mu = x)))
exp.CDU <- sapply(lambdaCDU, 
                  function(x) mean(rnbinom(1000, size = thetam6a, mu = x)))
exp.LIN <- sapply(lambdaLIN, 
                  function(x) mean(rnbinom(1000, size = thetam6a, mu = x)))
exp.SPD <- sapply(lambdaSPD, 
                  function(x) mean(rnbinom(1000, size = thetam6a, mu = x)))
exp.GRU <- sapply(lambdaGRU, 
                  function(x) mean(rnbinom(1000, size = thetam6a, mu = x)))

exp.values <- c(exp.FDP, exp.CDU, exp.LIN, exp.SPD, exp.GRU)
df.parties <- data.frame(exp.values)
df.parties$id <- c(rep("FDP", 1000), 
                   rep("CDU/CSU", 1000), 
                   rep("LINKE", 1000), 
                   rep("SPD", 1000),
                   rep("GRUENE", 1000))

ggplot(df.parties, aes(x = exp.values, fill = id)) + 
  geom_density(alpha = 0.4) + 
  guides(fill = guide_legend(title = "")) +
  xlab("Expected Speeches") +
  ylab("Density") +
  ggtitle("Simulated expected Number of Speeches, \n separated 
          by Party Affiliation") +
  theme_bw()

###########################################################
