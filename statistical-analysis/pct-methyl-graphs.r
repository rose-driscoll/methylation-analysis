### Code for Pct Methyl Figs 
rm(list=ls())
### function to calculate SEM for graph whiskers
sem <- function(x){sd(x,na.rm=T)/sqrt(length(x[!is.na(x)]))} 

methyl.mat <- read.csv("pct-methylation.csv", na.strings=c("NA","Na","na","N/A","N/a","n/a","","-"),stringsAsFactors=FALSE)

attach(methyl.mat)
### calculate mean methylation for each amplicon
m.AI <- apply(cbind(AI130, AI135, AI143, AI150, AI163, AI181, AI187, AI228),1,mean,na.rm=TRUE)
m.AP <- apply(cbind(AP145, AP159, AP188, AP279, AP317, AP326),1,mean,na.rm=TRUE)
m.BI <- apply(cbind(BI217, BI333, BI358, BI370),1,mean,na.rm=TRUE)
m.B2 <- apply(cbind(B2174, B2187, B2213, B2238, B2253, B2302, B2314),1,mean,na.rm=TRUE)
methyl.mat <- cbind(methyl.mat,m.AI,m.AP,m.BI,m.B2)
rm(m.AI,m.AP,m.BI,m.B2)
detach(methyl.mat);attach(methyl.mat)

#par(mar=c(5,4,4,2)); par(oma=c(3,3,3,3));   #this is default
par(mfrow=c(2,1))
par(mar=c(3,4,1,0.5)); par(oma=c(1,1,1,1))   #this is tuned for the jpg

means <- c(
           mean(m.AI[tissue=="body" & sex=="F"],na.rm=T),
           mean(m.AI[tissue=="body" & sex=="Y"],na.rm=T),
           mean(m.AI[tissue=="body" & sex=="R"],na.rm=T),
           mean(m.AI[tissue=="body" & sex=="B"],na.rm=T),
           NA,
           mean(m.AP[tissue=="body" & sex=="F"],na.rm=T),
           mean(m.AP[tissue=="body" & sex=="Y"],na.rm=T),
           mean(m.AP[tissue=="body" & sex=="R"],na.rm=T),
           mean(m.AP[tissue=="body" & sex=="B"],na.rm=T)
           )

sems <- c(
          sem(m.AI[tissue=="body" & sex=="F"]),
          sem(m.AI[tissue=="body" & sex=="Y"]),
          sem(m.AI[tissue=="body" & sex=="R"]),
          sem(m.AI[tissue=="body" & sex=="B"]),
          NA,
          sem(m.AP[tissue=="body" & sex=="F"]),
          sem(m.AP[tissue=="body" & sex=="Y"]),
          sem(m.AP[tissue=="body" & sex=="R"]),
          sem(m.AP[tissue=="body" & sex=="B"])
)

bars <- barplot(means,xlab="",ylab="Pct. methylation (gonads)",ylim=c(0,100),
names=c("Female","Yellow","Red","Green",NA,"Female","Yellow","Red","Green"),
                col=c("black","goldenrod1","red3","springgreen3",NA,"black","goldenrod1","red3","springgreen3"))
arrows(bars,means+sems,bars,means, length = 0.1, angle = 90, code = 1)
arrows(bars,means,bars,means-sems, length = 0.1, angle = 90, code = 2, col=c("white","black","black","black",NA,"white","black","black","black","white"))
mtext("AI amplicon",side=1,at=(bars[2]+bars[3])/2,line=2)
mtext("AP amplicon",side=1,at=(bars[7]+bars[8])/2,line=2)
text(0.1,95,"a)")
lines(c(bars[1]-0.3,bars[1]+0.3),c(100,100))
lines(c(bars[1]-0.3,bars[1]+0.3),c(97,97))
lines(c(bars[1]-0.3,bars[1]+0.3),c(94,94))
lines(c(bars[2]-0.3,bars[4]+0.3),c(100,100))
lines(c(bars[2]-0.3,bars[3]+0.3),c(97,97))
lines(c(bars[2]-0.3,bars[2]+0.3),c(94,94))
text((bars[1]+bars[2])/2,97,"***",cex=2)

lines(c(bars[6]-0.3,bars[6]+0.3),c(100,100))
lines(c(bars[6]-0.3,bars[6]+0.3),c(97,97))
lines(c(bars[6]-0.3,bars[6]+0.3),c(94,94))
lines(c(bars[7]-0.3,bars[9]+0.3),c(100,100))
lines(c(bars[7]-0.3,bars[8]+0.3),c(97,97))
lines(c(bars[7]-0.3,bars[7]+0.3),c(94,94))
text((bars[6]+bars[7])/2,97,"***",cex=2)

means <- c(
           mean(m.AI[tissue=="head" & sex=="F"],na.rm=T),
           mean(m.AI[tissue=="head" & sex=="Y"],na.rm=T),
           mean(m.AI[tissue=="head" & sex=="R"],na.rm=T),
           mean(m.AI[tissue=="head" & sex=="B"],na.rm=T),
           NA,
           mean(m.AP[tissue=="head" & sex=="F"],na.rm=T),
           mean(m.AP[tissue=="head" & sex=="Y"],na.rm=T),
           mean(m.AP[tissue=="head" & sex=="R"],na.rm=T),
           mean(m.AP[tissue=="head" & sex=="B"],na.rm=T))

sems <- c(
          sem(m.AI[tissue=="head" & sex=="F"]),
          sem(m.AI[tissue=="head" & sex=="Y"]),
          sem(m.AI[tissue=="head" & sex=="R"]),
          sem(m.AI[tissue=="head" & sex=="B"]),
          NA,
          sem(m.AP[tissue=="head" & sex=="F"]),
          sem(m.AP[tissue=="head" & sex=="Y"]),
          sem(m.AP[tissue=="head" & sex=="R"]),
          sem(m.AP[tissue=="head" & sex=="B"])
)
bars <- barplot(means,xlab="",ylab="Pct. methylation (brain)",ylim=c(0,100),
names=c("Female","Yellow","Red","Green",NA,"Female","Yellow","Red","Green"),
                col=c("black","goldenrod1","red3","springgreen3",NA,"black","goldenrod1","red3","springgreen3"))
arrows(bars,means+sems,bars,means, length = 0.1, angle = 90, code = 1)
arrows(bars,means,bars,means-sems, length = 0.1, angle = 90, code = 2, col=c("white","black","black","black",NA,"white","black","black","black","white"))
mtext("AI amplicon",side=1,at=(bars[2]+bars[3])/2,line=2)
mtext("AP amplicon",side=1,at=(bars[7]+bars[8])/2,line=2)
text(0.2,95,"b)")

dev.print(device=pdf,file="Fig5-pct-Amethylation-adults.pdf")
dev.print(device=jpeg,file="Fig5-pct-Amethylation-adults.jpg",height=600,width=900)

par(mfrow=c(2,1))
par(mar=c(3,4,1,0.5)); par(oma=c(1,1,1,1))   #this is tuned for the jpg

### Adult B pct methylation by sex/colour
means <- c(
           mean(m.BI[tissue=="body" & sex=="F"],na.rm=T),
           mean(m.BI[tissue=="body" & sex=="Y"],na.rm=T),
           mean(m.BI[tissue=="body" & sex=="R"],na.rm=T),
           mean(m.BI[tissue=="body" & sex=="B"],na.rm=T),
           NA,
           mean(m.B2[tissue=="body" & sex=="F"],na.rm=T),
           mean(m.B2[tissue=="body" & sex=="Y"],na.rm=T),
           mean(m.B2[tissue=="body" & sex=="R"],na.rm=T),
           mean(m.B2[tissue=="body" & sex=="B"],na.rm=T)
           )

sems <- c(
          sem(m.BI[tissue=="body" & sex=="F"]),
          sem(m.BI[tissue=="body" & sex=="Y"]),
          sem(m.BI[tissue=="body" & sex=="R"]),
          sem(m.BI[tissue=="body" & sex=="B"]),
          NA,
          sem(m.B2[tissue=="body" & sex=="F"]),
          sem(m.B2[tissue=="body" & sex=="Y"]),
          sem(m.B2[tissue=="body" & sex=="R"]),
          sem(m.B2[tissue=="body" & sex=="B"])
)

bars <- barplot(means,xlab="",ylab="Pct. methylation (gonads)",ylim=c(0,100),
                names=c("Female","Yellow","Red","Green",NA,"Female","Yellow","Red","Green"),
                col=c("black","goldenrod1","red3","springgreen3",NA,"black","goldenrod1","red3","springgreen3"))
arrows(bars,means+sems,bars,means, length = 0.1, angle = 90, code = 1)
arrows(bars,means,bars,means-sems, length = 0.1, angle = 90, code = 2, col=c("white","black","black","black",NA,"white","black","black","black","white"))
mtext("BI amplicon",side=1,at=(bars[2]+bars[3])/2,line=2)
mtext("BP amplicon",side=1,at=(bars[7]+bars[8])/2,line=2)
text(0.1,95,"b)")

means <- c(
           mean(m.BI[tissue=="head" & sex=="F"],na.rm=T),
           mean(m.BI[tissue=="head" & sex=="Y"],na.rm=T),
           mean(m.BI[tissue=="head" & sex=="R"],na.rm=T),
           mean(m.BI[tissue=="head" & sex=="B"],na.rm=T),
           NA,
           mean(m.B2[tissue=="head" & sex=="F"],na.rm=T),
           mean(m.B2[tissue=="head" & sex=="Y"],na.rm=T),
           mean(m.B2[tissue=="head" & sex=="R"],na.rm=T),
           mean(m.B2[tissue=="head" & sex=="B"],na.rm=T))

sems <- c(
          sem(m.BI[tissue=="head" & sex=="F"]),
          sem(m.BI[tissue=="head" & sex=="Y"]),
          sem(m.BI[tissue=="head" & sex=="R"]),
          sem(m.BI[tissue=="head" & sex=="B"]),
          NA,
          sem(m.B2[tissue=="head" & sex=="F"]),
          sem(m.B2[tissue=="head" & sex=="Y"]),
          sem(m.B2[tissue=="head" & sex=="R"]),
          sem(m.B2[tissue=="head" & sex=="B"])
)
bars <- barplot(means,xlab="",ylab="Pct. methylation (brain)",ylim=c(0,100),
                names=c("Female","Yellow","Red","Green",NA,"Female","Yellow","Red","Green"),
                col=c("black","goldenrod1","red3","springgreen3",NA,"black","goldenrod1","red3","springgreen3"))
arrows(bars,means+sems,bars,means, length = 0.1, angle = 90, code = 1)
arrows(bars,means,bars,means-sems, length = 0.1, angle = 90, code = 2, col=c("white","black","black","black",NA,"white","black","black","black","white"))
mtext("BI amplicon",side=1,at=(bars[2]+bars[3])/2,line=2)
mtext("BP amplicon",side=1,at=(bars[7]+bars[8])/2,line=2)
text(0.1,95,"b)")

dev.print(device=pdf,file="Fig6-pct-Bmethylation-adults.pdf")
dev.print(device=jpeg,file="Fig6-pct-Bmethylation-adults.jpg",height=600,width=900)

###
###
###  FRY
###
###
par(mfrow=c(2,1))
par(mar=c(3,4,1,0.5)); par(oma=c(1,1,1,1))   #this is tuned for the jpg

means <- c(mean(m.AI[tissue=="body" & sex=="Fry"],na.rm=T),
           mean(m.AI[tissue=="head" & sex=="Fry"],na.rm=T),
           NA,
           mean(m.AP[tissue=="body" & sex=="Fry"],na.rm=T),
           mean(m.AP[tissue=="head" & sex=="Fry"],na.rm=T)
           )

sems <- c(sem(m.AI[tissue=="body" & sex=="Fry"]),
          sem(m.AI[tissue=="head" & sex=="Fry"]),
          NA,
          sem(m.AP[tissue=="body" & sex=="Fry"]),
          sem(m.AP[tissue=="head" & sex=="Fry"])
          )
bars <- barplot(means-50,offset=50,xlab="",ylab="Pct. methylation",ylim=c(50,100),
                names=c("Body","Head",NA,"Body","Head"),col=c("white","grey",NA,"white","grey"))
arrows(bars,means-sems,bars,means+sems,length=0.05,angle=90,code=3)
mtext("AI",side=1,at=(bars[1]+bars[2])/2,line=2)
mtext("AP",side=1,at=(bars[4]+bars[5])/2,line=2)

lines(c(bars[1]-0.4,bars[1]+0.4),c(85,85))
lines(c(bars[2]-0.4,bars[2]+0.4),c(85,85))
text((bars[1]+bars[2])/2,85,"***",cex=2)


means <- c(mean(m.BI[tissue=="body" & sex=="Fry" & sex=="Fry"],na.rm=T),
           mean(m.BI[tissue=="head" & sex=="Fry" & sex=="Fry"],na.rm=T),
           NA,
           mean(m.B2[tissue=="body" & sex=="Fry" & sex=="Fry"],na.rm=T),
           mean(m.B2[tissue=="head" & sex=="Fry" & sex=="Fry"],na.rm=T)
           )

sems <- c(sem(m.BI[tissue=="body" & sex=="Fry" & sex=="Fry"]),
          sem(m.BI[tissue=="head" & sex=="Fry" & sex=="Fry"]),
          NA,
          sem(m.B2[tissue=="body" & sex=="Fry" & sex=="Fry"]),
          sem(m.B2[tissue=="head" & sex=="Fry" & sex=="Fry"])
          )

bars <- barplot(means-50,offset=50,xlab="",ylab="Pct. methylation",ylim=c(50,100),
                names=c("Body","Head",NA,"Body","Head"),col=c("white","grey",NA,"white","grey"))
arrows(bars,means-sems,bars,means+sems,length=0.05,angle=90,code=3)
mtext("BI",side=1,at=(bars[1]+bars[2])/2,line=2)
mtext("BP",side=1,at=(bars[4]+bars[5])/2,line=2)

dev.print(device=pdf,file="Fig7-pct-methylation-fry.pdf")
dev.print(device=jpeg,file="Fig7-pct-methylation-fry.jpg",height=700,width=600)
