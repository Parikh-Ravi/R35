### R35 manuscript

### subgroup analyses for EHR only algorithm model: cancer type, race, stage, age

library(dplyr)
library(tidyselect)

load('R35.fst.enc.phaseI.risk.Rdata')
load('table1.Rdata')

names(data_PL_fit)

dat <- merge(table1, data_PL_fit[c('yhat_P1','PAT_ID')], by = 'PAT_ID')

dat$subtype[dat$subtype == 'Missing'] <- NA


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
  
  
  # calculate sensitivty,ppv,fpr at cutoff 0.1
  ind <- (yhat>=0.1)
  num_sensitivity <- sum(yhat[ind==1])
  denom_sensitivity <- sum(yhat)
  sensitivity_01 <- num_sensitivity/denom_sensitivity
  num_ppv <- num_sensitivity
  denom_ppv <- sum(ind)
  ppv_01 <- num_ppv/denom_ppv
  num_fpr <- sum((1-yhat)[ind==1])
  denom_fpr <- sum(1-yhat)
  fpr_01 <- num_fpr/denom_fpr
  
  
  return(list(auc=auc,auprc = auprc,
              sensitivity_01=sensitivity_01,ppv_01 = ppv_01, fpr_01=fpr_01))
  
}

## overall

samplen = 1000
PERF <- data.frame(matrix(0, nrow = samplen, ncol = 1))
colnames(PERF) <- c("AUC")

set.seed(1223)

for (j in 1:samplen) {
  temp <- dat
  bootstrap.rows <- sample(x = 1:nrow(temp), size = nrow(temp), replace = T)
  bootstrap.dat <- temp[bootstrap.rows,]
  
  
  ## metrics 
  
  perf.m <- perf.func(bootstrap.dat$yhat_P1)
  
  perf.tmp <- data.frame(
    perf.m$auc)
  
  PERF[j,] <- perf.tmp
  
} 



perf <- unlist(perf.func(dat$yhat_P1))

perf.df <- data.frame(t(perf))


perf.df$perf.se <- apply(PERF, 2, sd, na.rm = TRUE)
perf.df$lc <- with(perf.df, auc+qnorm(0.025)*perf.se)
perf.df$hc <- with(perf.df, auc+qnorm(0.975)*perf.se)

perf.df$CI <- paste0('[',round(perf.df$lc,3), ',', round(perf.df$hc,3), ']')
perf.df$subgroup <- 'overall'

PhaseI.overall <- perf.df[c('subgroup','auc','CI','auprc','sensitivity_01','ppv_01','fpr_01')]



## cancer subtype

cancer <- unique(unlist(dat$subtype))
PERF.CANCER = list()
perf.se = list()
for (i in 1: length(cancer)) {
data <- subset(dat, subtype == cancer[i])

samplen = 1000
PERF <- data.frame(matrix(0, nrow = samplen, ncol = 1))
colnames(PERF) <- c("AUC")

set.seed(1223)

for (j in 1:samplen) {
  temp <- data
  bootstrap.rows <- sample(x = 1:nrow(temp), size = nrow(temp), replace = T)
  bootstrap.dat <- temp[bootstrap.rows,]
  
  
  ## metrics 
  
  perf.m <- perf.func(bootstrap.dat$yhat_P1)
  
  perf.tmp <- data.frame(
    perf.m$auc)
  
  PERF[j,] <- perf.tmp
  
   } 

  perf.se[[i]] <- apply(PERF, 2, sd, na.rm = TRUE)
}

PERF.CANCER = do.call(rbind, perf.se)

perf.cancer <- list()
for (i in 1: length(cancer)) {
  data <- subset(dat, subtype == cancer[i])
  perf.cancer[[i]] <- unlist(perf.func(data$yhat_P1))
}

perf.cancer.d <- do.call(rbind, perf.cancer)
perf.df <- cbind(data.frame(cancer), perf.cancer.d)


perf.df$perf.se <- PERF.CANCER
perf.df$lc <- with(perf.df, auc+qnorm(0.025)*perf.se)
perf.df$hc <- with(perf.df, auc+qnorm(0.975)*perf.se)

perf.df$CI <- paste0('[',round(perf.df$lc,3), ',', round(perf.df$hc,3), ']')

PhaseI.perf.cancer <- perf.df[c('cancer','auc','CI','auprc','sensitivity_01','ppv_01','fpr_01')]



### race subtype
dat$race <- with(dat, ifelse(RACE %in% c('Black'),'Black',
                             ifelse(RACE %in% c('White'), 'White',
                                                ifelse(RACE %in% c('American Indian','Asian','East Indian',
                                                                   'HLB-Hispanic Latino/Black','HLW-Hispanic Latino/White','Other',
                                                                   'Pacific Island'),'Other', NA))))

dat$race.c <- factor(dat$race, levels = c('Black','White','Other'))

race.dat <- na.omit(dat[c('race.c','yhat_P1')])

race <- unique(unlist(race.dat$race.c))
PERF.RACE = list()
perf.se = list()
for (i in 1: length(race)) {
  data <- subset(race.dat, race.c == race[i])
  
  samplen = 1000
  PERF <- data.frame(matrix(0, nrow = samplen, ncol = 1))
  colnames(PERF) <- c("AUC")
  
  set.seed(1223)
  
  for (j in 1:samplen) {
    temp <- data
    bootstrap.rows <- sample(x = 1:nrow(temp), size = nrow(temp), replace = T)
    bootstrap.dat <- temp[bootstrap.rows,]
    
    
    ## metrics 
    
    perf.m <- perf.func(bootstrap.dat$yhat_P1)
    
    perf.tmp <- data.frame(
      perf.m$auc)
    
    PERF[j,] <- perf.tmp
    
  } 
  
  perf.se[[i]] <- apply(PERF, 2, sd, na.rm = TRUE)
}

PERF.RACE = do.call(rbind, perf.se)

perf.race <- list()
for (i in 1: length(race)) {
  data <- subset(race.dat, race.c == race[i])
  perf.race[[i]] <- unlist(perf.func(data$yhat_P1))
}

perf.race.d <- do.call(rbind, perf.race)
perf.df <- cbind(data.frame(race), perf.race.d)


perf.df$perf.se <- PERF.RACE
perf.df$lc <- with(perf.df, auc+qnorm(0.025)*perf.se)
perf.df$hc <- with(perf.df, auc+qnorm(0.975)*perf.se)

perf.df$CI <- paste0('[',round(perf.df$lc,3), ',', round(perf.df$hc,3), ']')

PhaseI.perf.race <- perf.df[c('race','auc','CI','auprc','sensitivity_01','ppv_01','fpr_01')]

### stage subtype
stage1 <- c('Stage IB', 'Stage IA',
            'Stage IA1','Stage IA2',
            'Stage I', 
            'Stage IA3','Stage IC', 'Stage IB2', 
            'Stage IS', 'Stage IEB', 'Stage IB1',
            'Stage IE')

stage2 <- c('Stage IIA', 'Stage II',
            'Stage IIB', 
            'Stage IIC','Stage IIA2')

stage3 <- c('Stage IIIB', 'Stage IIIC', 
            'Stage III', 'Stage IIIA', 'Stage IIIC2', 
            'Stage IIIC1', 
            'Stage IIIA2', 'Stage IIIE',
            'Stage IIID', 'Stage IIIA1')

stage4 <- c('Stage IV', 'Stage IVB', 'Stage IVA',
            'Stage IVC', 'Stage IVA1')
dat$stage <- with(dat, ifelse(STAGE_GROUP %in% stage1, 'I', 
                              ifelse(STAGE_GROUP %in% stage2, 'II',
                                     ifelse(STAGE_GROUP %in% stage3, 'III', 
                                            ifelse(STAGE_GROUP %in% stage4, 'IV', NA)))))

dat$stage.b <- with(dat, ifelse(stage %in% c('I','II','III'), 'I-III',
                                ifelse(stage %in% c('IV'), 'IV', NA)))

dat$stage.b <- factor(dat$stage.b, levels = c('I-III','IV'))

stage.dat <- na.omit(dat[c('stage.b','yhat_P1')])


stage <- unique(unlist(stage.dat$stage.b))
PERF.stage = list()
perf.se = list()
for (i in 1: length(stage)) {
  data <- subset(stage.dat, stage.b == stage[i])
  
  samplen = 1000
  PERF <- data.frame(matrix(0, nrow = samplen, ncol = 1))
  colnames(PERF) <- c("AUC")
  
  set.seed(1223)
  
  for (j in 1:samplen) {
    temp <- data
    bootstrap.rows <- sample(x = 1:nrow(temp), size = nrow(temp), replace = T)
    bootstrap.dat <- temp[bootstrap.rows,]
    
    
    ## metrics 
    
    perf.m <- perf.func(bootstrap.dat$yhat_P1)
    
    perf.tmp <- data.frame(
      perf.m$auc)
    
    PERF[j,] <- perf.tmp
    
  } 
  
  perf.se[[i]] <- apply(PERF, 2, sd, na.rm = TRUE)
}

PERF.stage = do.call(rbind, perf.se)

perf.stage <- list()
for (i in 1: length(stage)) {
  data <- subset(stage.dat, stage.b == stage[i])
  perf.stage[[i]] <- unlist(perf.func(data$yhat_P1))
}

perf.stage.d <- do.call(rbind, perf.stage)
perf.df <- cbind(data.frame(stage), perf.stage.d)


perf.df$perf.se <- PERF.stage
perf.df$lc <- with(perf.df, auc+qnorm(0.025)*perf.se)
perf.df$hc <- with(perf.df, auc+qnorm(0.975)*perf.se)

perf.df$CI <- paste0('[',round(perf.df$lc,3), ',', round(perf.df$hc,3), ']')

PhaseI.perf.stage <- perf.df[c('stage','auc','CI','auprc','sensitivity_01','ppv_01','fpr_01')]

### age subtype
dat$age.b <- with(dat, ifelse(PAT_AGE >70, '>70','<=70'))

dat$age.b <- factor(dat$age.b, levels = c('>70','<=70'))

age <- unique(unlist(dat$age.b))
PERF.age = list()
perf.se = list()
for (i in 1: length(age)) {
  data <- subset(dat, age.b == age[i])
  
  samplen = 1000
  PERF <- data.frame(matrix(0, nrow = samplen, ncol = 1))
  colnames(PERF) <- c("AUC")
  
  set.seed(1223)
  
  for (j in 1:samplen) {
    temp <- data
    bootstrap.rows <- sample(x = 1:nrow(temp), size = nrow(temp), replace = T)
    bootstrap.dat <- temp[bootstrap.rows,]
    
    
    ## metrics 
    
    perf.m <- perf.func(bootstrap.dat$yhat_P1)
    
    perf.tmp <- data.frame(
      perf.m$auc)
    
    PERF[j,] <- perf.tmp
    
  } 
  
  perf.se[[i]] <- apply(PERF, 2, sd, na.rm = TRUE)
}

PERF.age = do.call(rbind, perf.se)

perf.age <- list()
for (i in 1: length(age)) {
  data <- subset(dat, age.b == age[i])
  perf.age[[i]] <- unlist(perf.func(data$yhat_P1))
}

perf.age.d <- do.call(rbind, perf.age)
perf.df <- cbind(data.frame(age), perf.age.d)


perf.df$perf.se <- PERF.age
perf.df$lc <- with(perf.df, auc+qnorm(0.025)*perf.se)
perf.df$hc <- with(perf.df, auc+qnorm(0.975)*perf.se)

perf.df$CI <- paste0('[',round(perf.df$lc,3), ',', round(perf.df$hc,3), ']')

PhaseI.perf.age <- perf.df[c('age','auc','CI','auprc','sensitivity_01','ppv_01','fpr_01')]

names(PhaseI.perf.cancer)[1] <- 'subgroup'
names(PhaseI.perf.race)[1] <- 'subgroup'
names(PhaseI.perf.stage)[1] <- 'subgroup'
names(PhaseI.perf.age)[1] <- 'subgroup'

PhaseI.perf <- rbind(PhaseI.overall,PhaseI.perf.cancer, PhaseI.perf.race, PhaseI.perf.stage, PhaseI.perf.age)

PhaseI.perf <- na.omit(PhaseI.perf)




