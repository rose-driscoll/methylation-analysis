## expression graphs
rm(list=ls())
expression.dat <- read.csv("expression_qPCR.csv")
attach(expression.dat)

### function to calculate SEM
sem <- function(x){sd(x,na.rm=T)/sqrt(length(x[!is.na(x)]))} 

par(mfrow=c(2,1))

par(oma=c(0,0,0,0))
par(mar=c(4,4,3,2))
means <- c(mean(rel_cyp19a1A[tissue=="body" & sex=="F"],na.rm=T),
           mean(rel_cyp19a1A[tissue=="body" & sex=="Y"],na.rm=T),
           mean(rel_cyp19a1A[tissue=="body" & sex=="R"],na.rm=T),
           NA,
           mean(rel_cyp19a1A[tissue=="head" & sex=="F"],na.rm=T),
           mean(rel_cyp19a1A[tissue=="head" & sex=="Y"],na.rm=T),
           mean(rel_cyp19a1A[tissue=="head" & sex=="R"],na.rm=T))

sems <- c(sem(rel_cyp19a1A[tissue=="body" & sex=="F"]),
          sem(rel_cyp19a1A[tissue=="body" & sex=="Y"]),
          sem(rel_cyp19a1A[tissue=="body" & sex=="R"]),
          NA,
          sem(rel_cyp19a1A[tissue=="head" & sex=="F"]),
          sem(rel_cyp19a1A[tissue=="head" & sex=="Y"]),
          sem(rel_cyp19a1A[tissue=="head" & sex=="R"]))

bars <- barplot(means,xlab="",ylab="Relative cyp19a1A expression",ylim=c(0,3),
                names=c("Females","Yellow","Red",NA,"Females","Yellow","Red"),
                col=c("black","goldenrod1","red3",NA,"black","goldenrod1","red3"))
arrows(bars,means+sems,bars,means, length = 0.1, angle = 90, code = 1)
arrows(bars,means,bars,means-sems, length = 0.1, angle = 90, code = 2, col=c("white","black","black",NA,"white","black","black","white"))
text(bars[1]-0.6,3,"A)",pos=1)
text(bars[5]-0.6,3,"B)",pos=1)
lines(c(bars[1],bars[2]),c(2.88,2.88))
text((bars[1]+bars[2])/2,2.94,"***",cex=2)
lines(c(bars[1],bars[3]),c(2.72,2.72))
text(bars[2],2.78,"***",cex=2)

mtext("Ovaries",side=1,at=bars[1],line=2)
mtext("Testes",side=1,at=(bars[2]+bars[3])/2,line=2)
mtext("Brain",side=1,at=bars[6],line=2)

means <- c(mean(rel_cyp19a1B[tissue=="body" & sex=="F"],na.rm=T),
           mean(rel_cyp19a1B[tissue=="body" & sex=="Y"],na.rm=T),
           mean(rel_cyp19a1B[tissue=="body" & sex=="R"],na.rm=T),
           NA,
           mean(rel_cyp19a1B[tissue=="head" & sex=="F"],na.rm=T),
           mean(rel_cyp19a1B[tissue=="head" & sex=="Y"],na.rm=T),
           mean(rel_cyp19a1B[tissue=="head" & sex=="R"],na.rm=T))

sems <- c(sem(rel_cyp19a1B[tissue=="body" & sex=="F"]),
          sem(rel_cyp19a1B[tissue=="body" & sex=="Y"]),
          sem(rel_cyp19a1B[tissue=="body" & sex=="R"]),
          NA,
          sem(rel_cyp19a1B[tissue=="head" & sex=="F"]),
          sem(rel_cyp19a1B[tissue=="head" & sex=="Y"]),
          sem(rel_cyp19a1B[tissue=="head" & sex=="R"]))

bars <- barplot(means,xlab="",ylab="Relative cyp19a1B expression",ylim=c(0,80),
                names=c("Females","Yellow","Red",NA,"Females","Yellow","Red"),
                col=c("black","goldenrod1","red3",NA,"black","goldenrod1","red3"))
arrows(bars,means+sems,bars,means, length = 0.1, angle = 90, code = 1)
arrows(bars,means,bars,means-sems, length = 0.1, angle = 90, code = 2, col=c("white","black","black",NA,"white","black","black","white"))
text(bars[1]-0.6,80,"C)",pos=1)
text(bars[5]-0.6,80,"D)",pos=1)

mtext("Ovaries",side=1,at=bars[1],line=2)
mtext("Testes",side=1,at=(bars[2]+bars[3])/2,line=2)
mtext("Brain",side=1,at=bars[6],line=2)

dev.print(device=pdf,file="Fig3-Adult-Expression-by-sex.pdf",height=8,width=8)
dev.print(device=jpeg,file="Fig3-Adult-Expression-by-sex.jpg",width=800)

###
###
###

par(mfrow=c(1,1))
par(mar=c(4,4,3,2))
means <- c(mean(rel_cyp19a1A[tissue=="body" & sex=="fry"],na.rm=T),
           mean(rel_cyp19a1A[tissue=="head" & sex=="fry"],na.rm=T),
           NA,
           mean(rel_cyp19a1B[tissue=="body" & sex=="fry"],na.rm=T),
           mean(rel_cyp19a1B[tissue=="head" & sex=="fry"],na.rm=T))
sems <- c(sem(rel_cyp19a1A[tissue=="body" & sex=="fry"]),
          sem(rel_cyp19a1A[tissue=="head" & sex=="fry"]),
          NA,
          sem(rel_cyp19a1B[tissue=="body" & sex=="fry"]),
          sem(rel_cyp19a1B[tissue=="head" & sex=="fry"]))
bars <- barplot(means,xlab="",ylim=c(0,0.8),ylab="Relative cyp19a1 expression",
                names=c("Trunk","Head",NA,"Trunk","Head"),
                col=c("white","grey",NA,"white","grey"))
arrows(bars,means-sems,bars,means+sems,length=0.1,angle=90,code=3)

lines(c(bars[4],bars[5]),c(0.73,0.73))
text((bars[4]+bars[5])/2,0.74,"***",cex=2)

lines(c(bars[1],bars[4]),c(0.19,0.19))
text((bars[1]+bars[4])/2,0.20,"***",cex=2)

lines(c(bars[2],bars[5]),c(0.78,0.78))
text((bars[2]+bars[5])/2,0.79,"***",cex=2)

mtext("A copy",side=1,at=(bars[1]+bars[2])/2,line=2)
mtext("B copy",side=1,at=(bars[4]+bars[5])/2,line=2)

dev.print(device=pdf,file="Fig4-Fry-Expression.pdf",width=8)
dev.print(device=jpeg,file="Fig4-Fry-Expression.jpg",width=800)
