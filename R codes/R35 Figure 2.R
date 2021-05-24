### R35 Figure2

### Blox plot for
### PRO only, Phase I data, and 2 Phase data
### AUC, AUPRC, Sensitivity (10% threshold), FPR (10% threshold)

library(MASS)
library(gam)
library(glmnet)
library(nleqslv)
library(MESS)
library(dplyr)


load('R35.fst.enc.phaseI.risk.Rdata')

## define performance metrics functions

perf.func <- function(yhat) { 
# AUC 
c <- c(seq(0, 1, by=0.01)*max(yhat),1)
denom_tpr <- sum(yhat)
denom_fpr <- sum(1-yhat)
num_tpr <- numeric(length(c))
num_fpr <- numeric(length(c))
for (j in 1:length(c)) {
  ind <- (yhat>=c[j])
  num_tpr[j] <- sum(yhat[ind==1])
  num_fpr[j] <- sum((1-yhat)[ind==1])
}
tpr <- num_tpr/denom_tpr
fpr <- num_fpr/denom_fpr
auc <- MESS::auc(fpr,tpr)

# AUPRC
c <- c(seq(0, 1, by=0.01)*max(yhat))
denom_tpr <- sum(yhat)
num_tpr_ppv <- numeric(length(c))
denom_ppv <- numeric(length(c))
for (j in 1:length(c)) {
  ind <- (yhat>=c[j])
  num_tpr_ppv[j] <- sum(yhat[ind==1])
  denom_ppv[j] <- sum(ind)
}
tpr <- num_tpr_ppv/denom_tpr
ppv <- num_tpr_ppv/denom_ppv
auprc <- MESS::auc(tpr,ppv)


# calculate sensitivty and fpr at cutoff 0.1
ind <- (yhat>=0.1)
num_sensitivity <- sum(yhat[ind==1])
denom_sensitivity <- sum(yhat)
sensitivity_01 <- num_sensitivity/denom_sensitivity
num_fpr <- sum((1-yhat)[ind==1])
denom_fpr <- sum(1-yhat)
fpr_01 <- num_fpr/denom_fpr


return(list(auc=auc,auprc = auprc,
            sensitivity_01=sensitivity_01,fpr_01=fpr_01))

}

## 1. PRO only 

PRO <- subset(data_PL_fit, R == 1)
PRO.names <- dput(names(PRO)[194:205])
PRO[PRO.names] <- lapply(PRO[PRO.names],as.numeric)
PRO.sub <- dput(names(PRO)[194:203])
PRO[PRO.sub] <- 1+PRO[PRO.sub]
summary(PRO[PRO.names])

logit <- glm(Y ~. - Y, data = PRO[194:206], family = "binomial")
PRO$yhat <- logit %>% predict(PRO[194:206], type = "response")

## bootstrapping
samplen = 1000
PERF <- data.frame(matrix(0, nrow = samplen, ncol = 4))
colnames(PERF) <- c("AUC","AUPRC","sensitivity","FPR")

set.seed(1223)

for (j in 1:samplen) {
  temp <- PRO
  bootstrap.rows <- sample(x = 1:nrow(temp), size = nrow(temp), replace = T)
  bootstrap.dat <- temp[bootstrap.rows,]
  
  
  ## metrics 
  
  perf.m <- perf.func(bootstrap.dat$yhat)
  
  perf.tmp <- data.frame(
    perf.m$auc, perf.m$auprc, perf.m$sensitivity_01, perf.m$fpr_01)
  
  PERF[j,] <- perf.tmp
  
} 


## calculate lower bounds, upper bounds, Q1, Q3
perf.se <- apply(PERF, 2, sd, na.rm = TRUE)

perf.df <- as.data.frame(unlist(perf.func(PRO$yhat)))

names(perf.df)[1] <- 'Perf'
perf.df$perf.se <- perf.se
perf.df$lc <- with(perf.df, Perf+qnorm(0.025)*perf.se)
perf.df$hc <- with(perf.df, Perf+qnorm(0.975)*perf.se)
perf.df$q1 <- with(perf.df, Perf+qnorm(0.25)*perf.se)
perf.df$q3 <- with(perf.df, Perf+qnorm(0.75)*perf.se)
perf.df$dat <- 'PRO'

PRO.perf <- perf.df

PRO.perf$perf.names <- c('AUC','AUPRC','Sensitivity', 'FPR')


## 2. Phase I data

## bootstrapping
samplen = 1000
PERF <- data.frame(matrix(0, nrow = samplen, ncol = 4))
colnames(PERF) <- c("AUC","AUPRC","sensitivity","FPR")

set.seed(1223)

for (j in 1:samplen) {
  temp <- data_PL_fit
  bootstrap.rows <- sample(x = 1:nrow(temp), size = nrow(temp), replace = T)
  bootstrap.dat <- temp[bootstrap.rows,]
  
  
  ## metrics 
  
  perf.m <- perf.func(bootstrap.dat$yhat_P1)
  
  perf.tmp <- data.frame(
    perf.m$auc, perf.m$auprc, perf.m$sensitivity_01, perf.m$fpr_01)
  
  PERF[j,] <- perf.tmp
  
} 


## calculate lower bounds, upper bounds, Q1, Q3
perf.se <- apply(PERF, 2, sd, na.rm = TRUE)

perf.df <- as.data.frame(unlist(perf.func(data_PL_fit$yhat_P1)))

names(perf.df)[1] <- 'Perf'
perf.df$perf.se <- perf.se
perf.df$lc <- with(perf.df, Perf+qnorm(0.025)*perf.se)
perf.df$hc <- with(perf.df, Perf+qnorm(0.975)*perf.se)
perf.df$q1 <- with(perf.df, Perf+qnorm(0.25)*perf.se)
perf.df$q3 <- with(perf.df, Perf+qnorm(0.75)*perf.se)
perf.df$dat <- 'EHR only'

PhaseI.perf <- perf.df

PhaseI.perf$perf.names <- c('AUC','AUPRC','Sensitivity', 'FPR')

## 3. 2 phase data
perf.se <- apply(PERF, 2, sd, na.rm = TRUE)

perf.df <- as.data.frame(unlist(perf.func(data_PL_fit$yhat_P1)))

Perf <- c(0.8551153, 0.395996, 0.6372404, 0.1098287)
perf.se <- c(0.002803448, 0.01332925, 0.0103728, 0.003541061)
perf.df <- data.frame(Perf, perf.se)
perf.df$lc <- with(perf.df, Perf+qnorm(0.025)*perf.se)
perf.df$hc <- with(perf.df, Perf+qnorm(0.975)*perf.se)
perf.df$q1 <- with(perf.df, Perf+qnorm(0.25)*perf.se)
perf.df$q3 <- with(perf.df, Perf+qnorm(0.75)*perf.se)
perf.df$dat <- 'EHR+PRO'

perf.df$perf.names <- c('AUC','AUPRC','Sensitivity', 'FPR')

PhaseII.perf <- perf.df

perf <- rbind(PRO.perf, PhaseI.perf, PhaseII.perf)

perf$dat <- factor(perf$dat, levels = c("PRO", "EHR only", "EHR+PRO"))
perf$perf.names <- factor(perf$perf.names, levels = c('AUC','AUPRC','Sensitivity', 'FPR'))
## draw plots
library(ggplot2)
ggplot(perf, aes(x = perf.names, fill= dat, order = dat)) + 
  #geom_errorbar(aes(ymin=lc, ymax=hc, colour = dat), width = 0.1, alpha = 0.5,stat = "identity") +
  geom_boxplot(aes(lower = q1,upper = q3, middle = Perf, ymin = lc, ymax = hc), 
               width = 0.5, alpha = 0.5, stat = "identity")+
  scale_y_continuous(limits = c(min(perf$lc), max(perf$hc)), n.breaks = 15)+
  labs(x = " ")+
  guides(fill=guide_legend(title=" "))+
  theme_bw()

scaleFUN <- function(x) sprintf("%.3f", x)

auc.plot <- subset(perf, perf.names == 'AUC')
ggplot(auc.plot, aes(x = dat, order = dat)) + 
  geom_errorbar(aes(ymin=lc, ymax=hc), width = 0.2, alpha = 1,stat = "identity") +
  geom_boxplot(aes(lower = q1,upper = q3, middle = Perf, ymin = lc, ymax = hc), 
               width = 0.5, alpha = 0.5, stat = "identity")+
  #stat_boxplot(geom = "errorbar", aes(ymin = lc, ymax = hc), width = 0.5) + 
  scale_y_continuous(limits = c(min(auc.plot$lc), max(auc.plot$hc)), n.breaks = 10, labels= scaleFUN)+
  labs(x = " ", title = "Figure2A. AUC")+
  guides(fill=guide_legend(title=" "))+
  theme_bw()

auprc.plot <- subset(perf, perf.names == 'AUPRC')
ggplot(auprc.plot, aes(x = dat, order = dat)) + 
  geom_errorbar(aes(ymin=lc, ymax=hc), width = 0.2, alpha = 1,stat = "identity") +
  geom_boxplot(aes(lower = q1,upper = q3, middle = Perf, ymin = lc, ymax = hc), 
               width = 0.5, alpha = 0.5, stat = "identity")+
  #stat_boxplot(geom = "errorbar", aes(ymin = lc, ymax = hc), width = 0.5) + 
  scale_y_continuous(limits = c(min(auprc.plot$lc), max(auprc.plot$hc)), n.breaks = 10,labels= scaleFUN)+
  labs(x = " ", title = "Figure2B. AUPRC")+
  #guides(fill=guide_legend(title=" "))+
  theme_bw()

sensitivity.plot <- subset(perf, perf.names == 'Sensitivity')
ggplot(sensitivity.plot, aes(x = dat, order = dat)) + 
  geom_errorbar(aes(ymin=lc, ymax=hc), width = 0.2, alpha = 1,stat = "identity") +
  geom_boxplot(aes(lower = q1,upper = q3, middle = Perf, ymin = lc, ymax = hc), 
               width = 0.5, alpha = 0.5, stat = "identity")+
  #stat_boxplot(geom = "errorbar", aes(ymin = lc, ymax = hc), width = 0.5) + 
  scale_y_continuous(limits = c(min(sensitivity.plot$lc), max(sensitivity.plot$hc)), n.breaks = 10,labels= scaleFUN)+
  labs(x = " ", title = "Figure2C. Sensitivity at 10% threshold")+
  #guides(fill=guide_legend(title=" "))+
  theme_bw()

FPR.plot <- subset(perf, perf.names == 'FPR')
ggplot(FPR.plot, aes(x = dat, order = dat)) + 
  geom_errorbar(aes(ymin=lc, ymax=hc), width = 0.2, alpha = 1,stat = "identity") +
  geom_boxplot(aes(lower = q1,upper = q3, middle = Perf, ymin = lc, ymax = hc), 
               width = 0.5, alpha = 0.5, stat = "identity")+
  #stat_boxplot(geom = "errorbar", aes(ymin = lc, ymax = hc), width = 0.5) + 
  scale_y_continuous(limits = c(min(FPR.plot$lc), max(FPR.plot$hc)), n.breaks = 10,labels= scaleFUN)+
  labs(x = " ", title = "Figure2D. FPR at 10% threshold")+
  #guides(fill=guide_legend(title=" "))+
  theme_bw()

