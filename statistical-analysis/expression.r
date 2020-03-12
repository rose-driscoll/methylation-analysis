###
## Analysis of qPCR expression data
###
rm(list=ls())
qPCR <- read.csv("expression_qPCR.csv")

library(car)
library(lme4);library(lmerTest)

qPCR.fry <- subset(qPCR,sex=="fry")
qPCR.adult <- subset(qPCR,sex!="fry")

adult.A <- aov(log(rel_cyp19a1A)~sex*tissue,data=qPCR.adult)
Anova(adult.A)
TukeyHSD(adult.A)

adult.B <- aov(log(rel_cyp19a1B)~sex*tissue,data=qPCR.adult)
Anova(adult.B)
TukeyHSD(adult.B)

fry.A.model<-lmer(data=qPCR.fry,log(rel_cyp19a1A)~tissue+(1|set_specific)+(1|plate))
summary(fry.A.model)

fry.B.model<-lmer(data=qPCR.fry,log(rel_cyp19a1B)~tissue+(1|set_specific)+(1|plate))
summary(fry.B.model)
