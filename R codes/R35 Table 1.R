
library(dplyr)
library(tidyselect)

load("r35.EHR.PRO.Rdata")

dat <- r35.EHR.PRO %>% group_by(PAT_ID) %>% mutate(count.enc = n())


# order by patient id/appointment time
dat = dat[order(dat$PAT_ID, dat$APPT_TIME),]
dat$DEATH_STATUS = as.numeric(dat$DEATH_DATE!="")



# subset of phase I variables
phase_I_vars <- c('PAT_ID','ALBUMIN..min', 'ALBUMIN..last', 'n_is_Tumor', 'n_is_METS', 'n_is_Tumor_recent', 'CREATININE..count', 
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
                  'X..MONOCYTES..count', 'X..MONOCYTES..mean','count.enc', 
                  'SEX_C_y','RACE','ECOG_strict','subtype','DEPARTMENT_NAME',
                  'PAT_AGE','STAGE_GROUP','fin_class_name','Pred')


# remove irrelevant variables
PRO_data = subset(dat, select = c(phase_I_vars,'DEATH_DATE','APPT_TIME','PAT_ID','ACTIVITIES & FUNCTION SCORING QUESTION',
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
data_Z <- select(data_Z, -`RASH YES/NO`)
# Create indicator R for selection into the phase II subset
R <- as.numeric(complete.cases(data_Z))
#0.5394186 have PRO data

# Define Y
Y <- data_first_encounters$death_180

# Create dataset for table 1
## Covariates
vars <- c('PAT_ID','count.enc', 'SEX_C_y','RACE','ECOG_strict','subtype','DEPARTMENT_NAME',
          'PAT_AGE','STAGE_GROUP','fin_class_name','Pred')

table1 <- subset(data_X, select = vars)
table1$Y <- Y
table1$R <- R

elix <- data_X[,grep("n_is_", names(data_X), value=F)]
table1$elix_count <- apply(X=elix,1,FUN=function(x) length(which(x>=1)))

## clean vars
cont.vars <- c('count.enc', 'PAT_AGE', 'elix_count')
table1[cont.vars] <- lapply(table1[cont.vars],as.numeric)
## age cat
table1$age.c <- with(table1, ifelse(PAT_AGE < 50, '< 50',
                                    ifelse(PAT_AGE >= 50 & PAT_AGE < 60, '50-60',
                                           ifelse(PAT_AGE >= 60 & PAT_AGE < 70, '60-70', '> 70'))))
## gender
table1$gender <- with(table1, ifelse(SEX_C_y == 1, 'Female', 'Male'))
## race
table1$race <- with(table1, ifelse(RACE == 'White', 'White',
                                   ifelse(RACE == 'Black', 'Black',
                                          ifelse(RACE %in% c('HLB-Hispanic Latino/Black',
                                                             'HLW-Hispanic Latino/White'), 'Hispanic', 'Other'))))
## ECOG
table1$ECOG.c <- with(table1, ifelse(ECOG_strict == 0, '0',
                                     ifelse(ECOG_strict == 1, '1', 
                                            ifelse(ECOG_strict >= 2, '>=2',
                                                   ifelse(is.na(ECOG_strict), 'Missing', NA)))))

library(tidyr)
table1$ECOG.c <- replace_na(table1$ECOG.c, 'Missing')

## cancer subtype
table1$subtype <- sub("^$", "Missing", table1$subtype)
table(table1$subtype)

## Tertiary/general -- all are tertiary academic 

# table1$hosp <- with(table1, ifelse(DEPARTMENT_NAME %in% c('GYNECOLOGY ONCOLOGY PERELMAN',
#                                                           'HEM ONC CENTER PERELMAN 2',
#                                                           'HEM ONC CENTER PERELMAN 3',
#                                                           'HEM ONC CENTER PERELMAN 4',
#                                                           'HEM ONC PLEURAL PROGRAM PMUC'), 'Tertiary academic',
#                                    'General oncology'))

## stage
table1$stage4 <- with(table1, ifelse(STAGE_GROUP %in% c('Stage IV', 'Stage IVB', 'Stage IVA',
                                                        'Stage IVC', 'Stage IVA1'), 'Stage IV', 'Not stage IV'))


## insurance
table1$insurance <- with(table1, ifelse(fin_class_name %in% c('Medicare', 'Managed Medicare'), 'Medicare',
                                        ifelse(fin_class_name %in% c('Medicaid (MA)', 'Managed Medicaid'), 'Medicaid',
                                               ifelse(fin_class_name %in% c('Managed Care'), 'Managed Care', 
                                                      ifelse(fin_class_name %in% c('Commercial', 'Blue Shield', 'Blue Cross'),
                                                             'Commercial Insurance',
                                                             ifelse(fin_class_name == 'Self-pay', 'Other', 'Missing'))))))


cat.var <- c('gender', 'race', 'ECOG.c', 'subtype', 'stage4','insurance')
table1[cat.var] <- lapply(table1[cat.var],as.factor)

save(table1, file = 'table1.RData')

library(tableone)

## Construct a table
vars <- c('count.enc', 'PAT_AGE', 'age.c','gender', 'elix_count','ECOG.c', 'race', 'insurance','stage4', 'subtype')
tab <- CreateTableOne(vars = vars, strata = "R", factorVars = cat.var, data = table1, test = T)
## Show table with SMD
tab1 <- print(tab,showAllLevels = TRUE, nonnormal = cont.vars,smd = TRUE)

setwd("/Users/manqingl/Box Sync/R35/Results")
write.csv(tab1, file = "Table1.csv")

### overall table
vars <- c('count.enc', 'PAT_AGE', 'age.c','gender', 'elix_count','ECOG.c', 'race', 'insurance','stage4', 'subtype')
tab <- CreateTableOne(vars = vars, factorVars = cat.var, data = table1, test = T)
## Show table with SMD
tab1 <- print(tab,showAllLevels = TRUE, nonnormal = cont.vars,smd = TRUE)

setwd("/Users/manqingl/Box Sync/R35/Results")
write.csv(tab1, file = "Table1.overall.csv")



