## R35 figures

### figure 1

load('R35.fst.enc.phaseI.risk.Rdata')

figure2 <- data_PL_fit

figure2$logodds <- log(figure2$yhat_P1/(1-figure2$yhat_P1))


names(figure2) <-  sub("SCORING ", "", names(figure2))
names(figure2) <-  sub("SCORE", "", names(figure2))
names(figure2) <-  sub("INTERFERENCE ", "", names(figure2))
names(figure2) <-  sub("FREQUENCY ", "", names(figure2))
names(figure2) <- sub("QUESTION", "",names(figure2))
names(figure2) <- sub("rash", "RASH",names(figure2))


names(figure2)[194:205] <- stringr::str_to_title(names(figure2)[194:205])
names(figure2)[204] <- 'Quality of Life'
names(figure2)[194] <- 'Performance Status'

var <- c("`Performance Status`", "`Anxiety `", "`Constipation  `", "`Decreased Appetite  `", 
         "`Diarrhea `", "`Fatigue `", "`Nausea `", "`Numbness & Tingling  `", 
         "`Sadness `", "`Shortness Of Breath `", "`Quality of Life`", "`Rash`"
)


PRO.names <- dput(names(figure2)[194:205])
figure2[PRO.names] <- lapply(figure2[PRO.names],as.numeric)
PRO.sub <- dput(names(figure2)[194:203])
figure2[PRO.sub] <- 1+figure2[PRO.sub]
summary(figure2[PRO.names])

## unadjusted 

regresults <-  lapply(var,
                      function(var) {
                        formula    <- as.formula(paste("Y ~", var))
                        res.logist <- glm(formula, data = figure2, family = binomial)
                        res.or <- exp(cbind(OR = coef(res.logist), confint(res.logist)))
                        res.or[2,]
                      })

allresults <- data.frame(do.call(rbind, regresults), row.names = var)
row.names(allresults) <-  sub("`", "", row.names(allresults))
row.names(allresults) <-  sub("`", "", row.names(allresults))

allresults <- round(allresults, 2)

unadOR <- allresults
unadOR$CI <- with(unadOR, paste0('[',X2.5..,',', X97.5..,']'))
unadOR <- unadOR[c(1,4)]


write.csv(unadOR, file = "unadOR.csv")


## adjusted 

regresults <-  lapply(var,
                      function(var) {
                        formula    <- as.formula(paste("Y ~logodds+", var))
                        res.logist <- glm(formula, data = figure2, family = binomial)
                        res.or <- exp(cbind(OR = coef(res.logist), confint(res.logist)))
                        res.or[3,]
                      })

allresults <- data.frame(do.call(rbind, regresults), row.names = var)
row.names(allresults) <-  sub("`", "", row.names(allresults))
row.names(allresults) <-  sub("`", "", row.names(allresults))

allresults$PRO <- row.names(allresults)
allresults$PRO <- factor(allresults$PRO,levels=rev(unique(allresults$PRO)))

allresults <- round(allresults, 2)

adOR <- allresults
adOR$CI <- with(adOR, paste0('[',X2.5..,',', X97.5..,']'))
adOR <- adOR[c(1,4)]


write.csv(adOR, file = "adOR.csv")



