## R35 table shell
## Table 2

## mean of each PRO across quartiles of phase I LASSO model risk 


library(MASS)
library(gam)
library(glmnet)
library(nleqslv)
library(MESS)

# read in the data
load("r35.EHR.PRO.Rdata")

# order by patient id/appointment time
r35.EHR.PRO = r35.EHR.PRO[order(r35.EHR.PRO$PAT_ID, r35.EHR.PRO$APPT_TIME),]
r35.EHR.PRO$DEATH_STATUS = as.numeric(r35.EHR.PRO$DEATH_DATE!="")

# subset of phase I variables
phase_I_vars <- c('ALBUMIN..min', 'ALBUMIN..last', 'n_is_Tumor', 'n_is_METS', 'n_is_Tumor_recent', 'CREATININE..count', 
                  'n_is_METS_recent', 'X..LYMPHOCYTES..last', 'ALKALINE.PHOSPHATASE..std', 'ALKALINE.PHOSPHATASE..max', 
                  'PLATELETS..last', 'CHLORIDE..last', 'PAT_AGE', 'LACTATE.DEHYDROGENASE..last', 
                  'ALKALINE.PHOSPHATASE..last', 'POTASSIUM..count', 'PLATELETS..mean', 'CANCER.ANTIGEN.19.9..last', 
                  'n_is_Leukemia.and.Myeloma', 'CARCINOEMBRYONIC.ANT..last', 'X..NEUTROPHILS..last', 'ALBUMIN..mean', 
                  'RED.BLOOD.CELLS..count', 'X..LYMPHOCYTES.MANUAL..count', 'RDW..max', 'RDW..mean', 
                  'RED.BLOOD.CELLS..first', 'HEMOGLOBIN..count', 'ALPHA.1.FETOPROTEINS..first', 'X..EOSINOPHILS..mean', 
                  'X..LYMPHOCYTES..mean', 'CANCER.ANTIGEN.19.9..max', 'HEMATOCRIT..count', 'X..NEUTROPHILS..last', 
                  'n_is_Myeloproliferative.neoplasms', 'X..LYMPHOCYTES..first', 'RED.BLOOD.CELLS..min', 
                  'X..LYMPHOCYTES..min', 'HEMOGLOBIN..last', 'PROTEIN.TOTAL..count', 'THYROID.STIMULATING.HORMONE..std', 
                  'ALT..count', 'MEAN.CELLULAR.HEMOGLOBIN..std', 'UREA.NITROGEN..max', 'R.AXIS..std', 'WBC..count', 
                  'UREA.NITROGEN..mean', 'SODIUM..min', 'CA.125..std', 'SODIUM..mean', 'SEX_C_x', 
                  'X..NEUTROPHILS.MANUAL..count', 'CARCINOEMBRYONIC.ANT..std', 'X..BAND.NEUTROPHILS..std', 
                  'X..BAND.NEUTROPHILS..count', 'RED.BLOOD.CELLS..last', 'CARCINOEMBRYONIC.ANT..min', 
                  'MEAN.CELLULAR.HEMOGLOBIN.CONCENTRATION..count', 'ALPHA.1.FETOPROTEINS..last', 'SODIUM..last', 
                  'RDW..std', 'ANION.GAP..last', 'LACTATE..std', 'X..LYMPHOCYTES..mean', 'X..BAND.NEUTROPHILS..last', 
                  'ALT..last', 'MEAN.CELLULAR.VOLUME..std', 'CANCER.ANTIGEN.19.9..min', 'ALKALINE.PHOSPHATASE..mean', 
                  'CHLORIDE..mean', 'HEMATOCRIT..last', 'FERRITIN..mean', 'AST..last', 'X..MONOCYTES..last', 'n_is_WL', 
                  'WBC..std', 'ALBUMIN..std', 'CANCER.ANTIGEN.19.9..mean', 'RDW..first', 'X..EOSINOPHILS..last', 
                  'PROTEIN.TOTAL..last', 'RDW..count', 'X..BAND.NEUTROPHILS..mean', 'VENTRICULAR.RATE..min', 
                  'LACTATE.DEHYDROGENASE..first', 'PLATELETS..first', 'CARBON.DIOXIDE..max', 'n_is_LD_recent', 
                  'QRS.DURATION..std', 'X..MONOCYTES.MANUAL..max', 'UA.SPECIFIC.GRAVITY..first', 
                  'CARCINOEMBRYONIC.ANT..first', 'BILIRUBIN.INDIRECT..last', 'CARBON.DIOXIDE..min', 
                  'X..BAND.NEUTROPHILS..mean', 'LACTATE.DEHYDROGENASE..max', 'CALCIUM..last', 'RED.BLOOD.CELLS..mean', 
                  'RED.BLOOD.CELLS..max', 'PLATELETS..std', 'X..MONOCYTES..max', 'CARCINOEMBRYONIC.ANT..mean', 
                  'VENTRICULAR.RATE..max', 'AST..max', 'ALKALINE.PHOSPHATASE..first', 'X..NEUTROPHILS.MANUAL..mean', 
                  'THYROID.STIMULATING.HORMONE..last', 'CREATININE..first', 'AST..mean', 'RDW..min', 
                  'CARBOXYHEMOGLOBIN..min', 'X..NEUTROPHILS..mean', 'X..MYELOCYTES..mean', 'X..NEUTROPHILS..min', 
                  'CANCER.ANTIGEN.15.3..first', 'T.AXIS..std', 'Q.T.INTERVAL..std', 'n_is_Coag', 'ALT..max', 
                  'n_is_Lymp_recent', 'X..BAND.NEUTROPHILS..max', 'PTT..first', 'n_is_Weakness', 'CALCIUM..mean', 
                  'FERRITIN..max', 'GLUCOSE.POINT.OF.CARE..max', 'X..BASOPHILS.MANUAL..mean', 'FIBRINOGEN..min', 
                  'CALCIUM..first', 'X..LYMPHOCYTES.MANUAL..last', 'UREA.NITROGEN..std', 'X..LYMPHOCYTES..count', 
                  'Q.T.INTERVAL..max', 'MEAN.CELLULAR.HEMOGLOBIN.CONCENTRATION..max', 'CALCIUM..std', 
                  'X..MONOCYTES.MANUAL..last', 'CALCIUM.IONIZED..first', 'LACTATE..last', 'n_is_Lymp', 
                  'ALPHA.1.FETOPROTEINS..max', 'URIC.ACID..mean', 'URIC.ACID..count', 'UA.SPECIFIC.GRAVITY..count', 
                  'X..EOSINOPHILS..min', 'PLATELETS..min', 'n_is_WL_recent', 'WBC..last', 'X..EOSINOPHILS..mean', 
                  'X..LYMPHOCYTES..max', 'P.R.INTERVAL..std', 'PROTEIN.TOTAL..first', 'X..EOSINOPHILS..std', 
                  'WBC..max', 'CANCER.ANTIGEN.15.3..last', 'CA.125..last', 'Q.T.INTERVAL..first', 
                  'CREATININE..min', 'FERRITIN..first', 'UA.SPECIFIC.GRAVITY..std', 'RDW..last', 'AST..first', 
                  'ALT..min', 'PCO2.ARTERIAL.TEMP.CORRECTED..count', 'MAGNESIUM..min', 'CREATININE..std', 
                  'CHLORIDE..first', 'n_is_Fluid_recent', 'n_is_VD_recent', 'Q.T.INTERVAL..last', 
                  'X..LYMPHOCYTES..first', 'UREA.NITROGEN..first', 'HEMATOCRIT..mean', 'SODIUM..std', 
                  'n_is_Paralysis_recent', 'X..NEUTROPHILS.MANUAL..last', 'X..MONOCYTES..max', 'n_is_VD', 
                  'X..MONOCYTES..mean', 'MEAN.CELLULAR.HEMOGLOBIN.CONCENTRATION..mean', 'PT..min', 
                  'BILIRUBIN.INDIRECT..min', 'BILIRUBIN.DIRECT..count', 'n_is_LD', 'X..BASOPHILS.MANUAL..min', 
                  'P.AXIS..max', 'PTT..min', 'ALT..first', 'ANION.GAP..std', 'PLATELETS..max', 'HEMATOCRIT..first', 
                  'X..MONOCYTES..count', 'X..MONOCYTES..mean')


# remove irrelevant variables
PRO_data = subset(r35.EHR.PRO, select = c(phase_I_vars,'DEATH_DATE','APPT_TIME','PAT_ID','ACTIVITIES & FUNCTION SCORING QUESTION',
                                          'ANXIETY INTERFERENCE SCORE','CONSTIPATION  SCORE','DECREASED APPETITE  SCORE',
                                          'DIARRHEA FREQUENCY SCORE','FATIGUE INTERFERENCE SCORE','NAUSEA FREQUENCY SCORE',
                                          'NUMBNESS & TINGLING  SCORE','RASH YES/NO','SADNESS INTERFERENCE SCORE',
                                          'SHORTNESS OF BREATH INTERFERENCE SCORE','GLOBAL02 SCORING QUESTION')) # count has all 1

# Create outcome variable of death within 180 of appointment
PRO_data$DEATH_DATE <- as.Date(PRO_data$DEATH_DATE, format = "%Y-%m-%d %H:%M:%S")
PRO_data$APPT_TIME <- as.Date(PRO_data$APPT_TIME, format = "%Y-%m-%d %H:%M:%S")
PRO_data$death_days_aft_appt <- PRO_data$DEATH_DATE-PRO_data$APPT_TIME
PRO_data$death_180 <- ifelse(PRO_data$death_days_aft_appt<=180,1,0)
PRO_data$death_180[which(is.na(PRO_data$death_180))] <- 0

# create a dataset that just includes the first/last encounter for each patient so we have one line per patient
last_encounter_pat <- NA
first_encounter_pat <- NA
j <- 1
for (i in unique(PRO_data$PAT_ID)) {
  last_encounter_pat[j] <- max(which(PRO_data$PAT_ID==i))
  first_encounter_pat[j] <- min(which(PRO_data$PAT_ID==i))
  j <- j+1
}
data_first_encounters <- PRO_data[first_encounter_pat,]  

# look for patients with little data before first encounter (n_ variables are 0 and labs are 0)
lab_counts <- c('CREATININE..count','POTASSIUM..count','RED.BLOOD.CELLS..count','X..LYMPHOCYTES.MANUAL..count','HEMOGLOBIN..count',
                'HEMATOCRIT..count', 'PROTEIN.TOTAL..count','ALT..count','WBC..count','X..NEUTROPHILS.MANUAL..count',
                'X..BAND.NEUTROPHILS..count', 'MEAN.CELLULAR.HEMOGLOBIN.CONCENTRATION..count',  'RDW..count', 'X..LYMPHOCYTES..count',
                'URIC.ACID..count', 'UA.SPECIFIC.GRAVITY..count', 'PCO2.ARTERIAL.TEMP.CORRECTED..count', 'BILIRUBIN.DIRECT..count',
                'X..MONOCYTES..count')
lab <- subset(data_first_encounters,select=lab_counts)
length(which(rowSums(lab)==0))
patients_no_labs <- which(rowSums(lab)==0)
data_first_encounters <- data_first_encounters[-patients_no_labs,]

# Create phase I data
data_X <- subset(data_first_encounters,select=phase_I_vars)



# Create phase II data
library(dplyr)
data_Z <- subset(data_first_encounters,select=c(`ACTIVITIES & FUNCTION SCORING QUESTION`,
                                                `ANXIETY INTERFERENCE SCORE`,`CONSTIPATION  SCORE`,`DECREASED APPETITE  SCORE`,
                                                `DIARRHEA FREQUENCY SCORE`,`FATIGUE INTERFERENCE SCORE`,`NAUSEA FREQUENCY SCORE`,
                                                `NUMBNESS & TINGLING  SCORE`,`RASH YES/NO`,`SADNESS INTERFERENCE SCORE`,
                                                `SHORTNESS OF BREATH INTERFERENCE SCORE`,`GLOBAL02 SCORING QUESTION`))
data_Z$rash <- with(data_Z, ifelse(`RASH YES/NO` == 'Y', 1, 0))
data_Z <- subset(data_Z, select = -c(`RASH YES/NO`))

# Create indicator R for selection into the phase II subset
R <- as.numeric(complete.cases(data_Z))
#0.5394186 have PRO data

# Define Y
Y <- data_first_encounters$death_180
PAT_ID <- data_first_encounters$PAT_ID
# Create dataset to use for PL Fit function
data_PL_fit <- cbind(PAT_ID,data_X,data_Z,Y,R)

dim_X <- ncol(data_X)
dim_Z <- ncol(data_Z)


vars <- c(194:205)
data_PL_fit[vars] <- sapply(data_PL_fit[vars], as.character)
data_PL_fit[vars] <- sapply(data_PL_fit[vars], as.numeric)
data_PL_fit[,204] = 6-data_PL_fit[,204]

data_PL_fit[vars] <- sapply(data_PL_fit[vars], as.factor)


## Run Lasso on phase I data
# HighD_X is a matrix of predictors
HighD_X <- as.matrix(data_PL_fit[c(2:193)])
# Here we fit the model using cross validation to choose the optimal penalty term
cv.fit <- cv.glmnet(HighD_X,data_PL_fit$Y,family="binomial")
# Here we obtain the predicted probabilities 
yhat_P1 <- predict(cv.fit,HighD_X,type="response",s="lambda.min")

data_PL_fit$yhat_P1 <- yhat_P1

setwd("/Users/manqingl/Box Sync/R35/Clean data")
save(data_PL_fit, file = 'R35.fst.enc.phaseI.risk.Rdata')

table2 <- data_PL_fit
#table2$yhat_P1 <- yhat_P1
table2$risk.q <- with(table2, ifelse(yhat_P1 <= 0.014, 'Q1',
                                     ifelse(yhat_P1 > 0.014 & yhat_P1 <= 0.026, 'Q2',
                                            ifelse(yhat_P1 > 0.026 & yhat_P1 <= 0.058, 'Q3', 'Q4'))))


PRO.names <- dput(names(table2)[194:205])
table2[PRO.names] <- lapply(table2[PRO.names],as.numeric)
PRO.sub <- dput(names(table2)[194:203])
table2[PRO.sub] <- 1+table2[PRO.sub]
summary(table2[PRO.names])
names(table2) <-  sub("SCORING ", "", names(table2))
names(table2) <-  sub("SCORE", "", names(table2))
names(table2) <-  sub("INTERFERENCE ", "", names(table2))
names(table2) <-  sub("FREQUENCY ", "", names(table2))
names(table2) <- sub("QUESTION", "",names(table2))
names(table2) <- sub("rash", "RASH",names(table2))


names(table2)[194:205] <- stringr::str_to_title(names(table2)[194:205])
names(table2)[204] <- 'Quality of Life'
names(table2)[194] <- 'Performance Status'
PRO.names <- dput(names(table2)[194:205])


library(tableone)

## Construct a table
tab <- CreateTableOne(vars = PRO.names, strata = "risk.q", data = table2, test = F)
## Show table with SMD
tab2 <- print(tab,showAllLevels = TRUE, smd = F)

## trend test (Mann-Kendall Trend Test)
library(data.table)
tab2.t <- as.data.frame(setDT(table2[c(194:205, 209)])[, lapply(.SD, mean, na.rm = T), keyby = risk.q])

tab2.t$risk.n <- with(tab2.t, ifelse(risk.q == 'Q1', 1,
                                     ifelse(risk.q == 'Q2',2,
                                            ifelse(risk.q == 'Q3', 3, 4))))

var <- c("`Performance Status`", "`Anxiety `", "`Constipation  `", "`Decreased Appetite  `", 
         "`Diarrhea `", "`Fatigue `", "`Nausea `", "`Numbness & Tingling  `", 
         "`Sadness `", "`Shortness Of Breath `", "`Quality of Life`", "`Rash`"
)
regresults <-  lapply(var,
                      function(var) {
                        formula    <- as.formula(paste("risk.n ~", var))
                        res.logist <- glm(formula, data = tab2.t, family = gaussian)
                        p.value <- coef(summary(res.logist))[2,4]
                        #res.or <- exp(cbind(OR = coef(res.logist), confint(res.logist)))
                        #res.or[2,]
                      })

allresults <- data.frame(do.call(rbind, regresults), row.names = var)
row.names(allresults) <-  sub("`", "", row.names(allresults))
row.names(allresults) <-  sub("`", "", row.names(allresults))

write.csv(allresults, file = "Table2.trend.csv")
write.csv(tab2, file = "Table2.csv")

res.logist <- glm(risk.n ~ `Nausea `, data = tab2.t, family = gaussian)
summary(res.logist)


