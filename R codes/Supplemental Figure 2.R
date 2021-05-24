## R35 figures

### supplemental figure 1
### univariate association between PRO trends and downstream mortality

load('pro.slope.RData')

var <- c("`Performance Status`", "`Anxiety `", "`Constipation  `", "`Decreased Appetite  `", 
         "`Diarrhea `", "`Fatigue `", "`Nausea `", "`Numbness & Tingling  `", 
         "`Sadness `", "`Shortness Of Breath `", "`Quality of Life`", "`Rash `"
)

regresults <-  lapply(var,
                      function(var) {
                        formula    <- as.formula(paste("Y ~", var))
                        res.logist <- glm(formula, data = pro.slope, family = binomial)
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


library(ggplot2)
library(forcats)

allresults$PRO <- row.names(allresults)
allresults$PRO2 <- factor(allresults$PRO,levels=c("Performance Status", "Quality of Life", "Fatigue ",  "Shortness Of Breath ",
                                                  "Anxiety ",  "Sadness ","Constipation  ", "Decreased Appetite  ", 
                                                  "Diarrhea ", "Nausea ", "Numbness & Tingling  ", 
                                                  "Rash "
))

allresults$PRO2 <- factor(allresults$PRO2, levels=rev(levels(allresults$PRO2)))

allresults %>% 
  ggplot(aes(OR,PRO2)) +
  # geom_point(aes(row.names(COR1),correlation)) +
  geom_point(size = 1.5)+
  geom_linerange(aes(xmin = X2.5.., xmax = X97.5..))+
  geom_vline(xintercept = 1, color = "grey")+
  ylab("") +
  xlab("Unadjusted Odds Ratio")+
  #ylim(0.05,0.4)+
  scale_x_continuous(breaks=seq(0.5,1.25, by = 0.05))+
  theme_bw()
