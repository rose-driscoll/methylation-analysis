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

bars <- barplot(means,xlab="",ylab="Percent methylation (gonads)",ylim=c(0,100),
#col="white",
names=c("Female","Yellow","Red","Green",NA,"Female","Yellow","Red","Green"),
                col=c("black","goldenrod1","red3","springgreen3",NA,"black","goldenrod1","red3","springgreen3"))
arrows(bars,means+sems,bars,means, length = 0.1, angle = 90, code = 1)
arrows(bars,means,bars,means-sems, length = 0.1, angle = 90, code = 2, col=c("white","black","black","black",NA,"white","black","black","black","white"))
mtext("AI amplicon",side=1,at=(bars[2]+bars[3])/2,line=2)
mtext("AP amplicon",side=1,at=(bars[7]+bars[8])/2,line=2)
text(bars[1]-0.6,95,"A)")
text(bars[6]-0.6,95,"B)")
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
bars <- barplot(means,xlab="",ylab="Percent methylation (brain)",ylim=c(0,100),
names=c("Female","Yellow","Red","Green",NA,"Female","Yellow","Red","Green"),
                col=c("black","goldenrod1","red3","springgreen3",NA,"black","goldenrod1","red3","springgreen3"))
arrows(bars,means+sems,bars,means, length = 0.1, angle = 90, code = 1)
arrows(bars,means,bars,means-sems, length = 0.1, angle = 90, code = 2, col=c("white","black","black","black",NA,"white","black","black","black","white"))
mtext("AI amplicon",side=1,at=(bars[2]+bars[3])/2,line=2)
mtext("AP amplicon",side=1,at=(bars[7]+bars[8])/2,line=2)
text(bars[1]-0.6,95,"C)")
text(bars[6]-0.6,95,"D)")

dev.print(device=jpeg,file="Fig5-pct-Amethylation-adults.jpg",height=600,width=900) ## this is pretty
dev.print(device=pdf,file="Fig5-pct-Amethylation-adults.pdf", height=6,width=9) ## this works well

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

bars <- barplot(means,xlab="",ylab="Percent methylation (gonads)",ylim=c(0,100),
                names=c("Female","Yellow","Red","Green",NA,"Female","Yellow","Red","Green"),
                col=c("black","goldenrod1","red3","springgreen3",NA,"black","goldenrod1","red3","springgreen3"))
arrows(bars,means+sems,bars,means, length = 0.1, angle = 90, code = 1)
arrows(bars,means,bars,means-sems, length = 0.1, angle = 90, code = 2, col=c("white","black","black","black",NA,"white","black","black","black","white"))
mtext("BI amplicon",side=1,at=(bars[2]+bars[3])/2,line=2)
mtext("BP amplicon",side=1,at=(bars[7]+bars[8])/2,line=2)
text(bars[1]-0.6,95,"A)")
text(bars[6]-0.6,95,"B)")

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
bars <- barplot(means,xlab="",ylab="Percent methylation (brain)",ylim=c(0,100),
                names=c("Female","Yellow","Red","Green",NA,"Female","Yellow","Red","Green"),
                col=c("black","goldenrod1","red3","springgreen3",NA,"black","goldenrod1","red3","springgreen3"))
arrows(bars,means+sems,bars,means, length = 0.1, angle = 90, code = 1)
arrows(bars,means,bars,means-sems, length = 0.1, angle = 90, code = 2, col=c("white","black","black","black",NA,"white","black","black","black","white"))
mtext("BI amplicon",side=1,at=(bars[2]+bars[3])/2,line=2)
mtext("BP amplicon",side=1,at=(bars[7]+bars[8])/2,line=2)
text(bars[1]-0.6,95,"C)")
text(bars[6]-0.6,95,"D)")

dev.print(device=pdf,file="Fig6-pct-Bmethylation-adults.pdf",height=6,width=9)
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
bars <- barplot(means,xlab="",ylab="Percent methylation",ylim=c(0,100),
                names=c("Body","Head",NA,"Body","Head"),col=c("white","grey",NA,"white","grey"))
arrows(bars,means-sems,bars,means+sems,length=0.05,angle=90,code=3)
mtext("AI",side=1,at=(bars[1]+bars[2])/2,line=2)
mtext("AP",side=1,at=(bars[4]+bars[5])/2,line=2)
text(bars[1]-0.55,100,"A)",pos=1)
text(bars[4]-0.6,100,"B)",pos=1)

lines(c(bars[1],bars[2]),c(84,84))
text((bars[1]+bars[2])/2,86,"***",cex=2)


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

bars <- barplot(means,xlab="",ylab="Percent methylation",ylim=c(0,100),
                names=c("Body","Head",NA,"Body","Head"),col=c("white","grey",NA,"white","grey"))
arrows(bars,means-sems,bars,means+sems,length=0.05,angle=90,code=3)
mtext("BI",side=1,at=(bars[1]+bars[2])/2,line=2)
mtext("BP",side=1,at=(bars[4]+bars[5])/2,line=2)
text(bars[1]-0.55,100,"C)",pos=1)
text(bars[4]-0.6,100,"D)",pos=1)

dev.print(device=pdf,file="Fig7-pct-methylation-fry.pdf",height=7,width=6)
dev.print(device=jpeg,file="Fig7-pct-methylation-fry.jpg",height=700,width=600)



####
#### Figure 8
####
rm(list=ls())
source("pct-methyl.r")

par(mfrow=c(1,2))
plot(sex.diff.dat.A.gonads$d,fry.pca.dat.A.body$PC1,xlab="Sex difference effect size",ylab="PC1 of fry body epialle variation")
text(sex.diff.dat.A.gonads["percent.CCCCCC",]$d,fry.pca.dat.A.body["percent.CCCCCC",]$PC1,"A",pos=1)
text(sex.diff.dat.A.gonads["percent.CCCTTT",]$d,fry.pca.dat.A.body["percent.CCCTTT",]$PC1,"B",pos=3)
text(sex.diff.dat.A.gonads["percent.CTCCCC",]$d,fry.pca.dat.A.body["percent.CTCCCC",]$PC1,"C",pos=1)

text(sex.diff.dat.A.gonads["percent.CCCCCCCC",]$d,fry.pca.dat.A.body["percent.CCCCCCCC",]$PC1,"D",pos=3)
text(sex.diff.dat.A.gonads["percent.TTTTTTTT",]$d,fry.pca.dat.A.body["percent.TTTTTTTT",]$PC1,"E",pos=3)
text(sex.diff.dat.A.gonads["percent.TCTTTTTT",]$d,fry.pca.dat.A.body["percent.TCTTTTTT",]$PC1,"F",pos=3)
text(sex.diff.dat.A.gonads["percent.CCCCCCTC",]$d,fry.pca.dat.A.body["percent.CCCCCCTC",]$PC1,"G",pos=1)
text(sex.diff.dat.A.gonads["percent.TCCCCCCC",]$d,fry.pca.dat.A.body["percent.TCCCCCCC",]$PC1,"H",pos=3)
text(sex.diff.dat.A.gonads["percent.CTTTTTTT",]$d,fry.pca.dat.A.body["percent.CTTTTTTT",]$PC1,"I",pos=3)
text(sex.diff.dat.A.gonads["percent.CCCTCCCC",]$d,fry.pca.dat.A.body["percent.CCCTCCCC",]$PC1,"J",pos=3)
text(sex.diff.dat.A.gonads["percent.CTCCCCCC",]$d,fry.pca.dat.A.body["percent.CTCCCCCC",]$PC1,"K",pos=3)
text(sex.diff.dat.A.gonads["percent.CCCCTCCC",]$d,fry.pca.dat.A.body["percent.CCCCTCCC",]$PC1,"L",pos=3)

text(sex.diff.dat.A.gonads["percent.CCCCTT",]$d,fry.pca.dat.A.body["percent.CCCCTT",]$PC1,"M",pos=3)
text(sex.diff.dat.A.gonads["percent.CCCCCT",]$d,fry.pca.dat.A.body["percent.CCCCCT",]$PC1,"N",pos=1)
text(sex.diff.dat.A.gonads["percent.CCTTTT",]$d,fry.pca.dat.A.body["percent.CCTTTT",]$PC1,"P",pos=1)

plot(sex.diff.dat.A.gonads$d,fry.pca.dat.A.body$PC2,xlab="Sex difference effect size",ylab="PC2 of fry body epialle variation")
text(sex.diff.dat.A.gonads["percent.CCCCCC",]$d,fry.pca.dat.A.body["percent.CCCCCC",]$PC2,"A",pos=3)
text(sex.diff.dat.A.gonads["percent.CCCTTT",]$d,fry.pca.dat.A.body["percent.CCCTTT",]$PC2,"B",pos=1)
text(sex.diff.dat.A.gonads["percent.CTCCCC",]$d,fry.pca.dat.A.body["percent.CTCCCC",]$PC2,"C",pos=3)

text(sex.diff.dat.A.gonads["percent.CCCCCCCC",]$d,fry.pca.dat.A.body["percent.CCCCCCCC",]$PC2,"D",pos=3)
text(sex.diff.dat.A.gonads["percent.TTTTTTTT",]$d,fry.pca.dat.A.body["percent.TTTTTTTT",]$PC2,"E",pos=3)
text(sex.diff.dat.A.gonads["percent.TCTTTTTT",]$d,fry.pca.dat.A.body["percent.TCTTTTTT",]$PC2,"F",pos=1)
text(sex.diff.dat.A.gonads["percent.CCCCCCTC",]$d,fry.pca.dat.A.body["percent.CCCCCCTC",]$PC2,"G",pos=3)
text(sex.diff.dat.A.gonads["percent.TCCCCCCC",]$d,fry.pca.dat.A.body["percent.TCCCCCCC",]$PC2,"H",pos=3)
text(sex.diff.dat.A.gonads["percent.CTTTTTTT",]$d,fry.pca.dat.A.body["percent.CTTTTTTT",]$PC2,"I",pos=3)
text(sex.diff.dat.A.gonads["percent.CCCTCCCC",]$d,fry.pca.dat.A.body["percent.CCCTCCCC",]$PC2,"J",pos=3)
text(sex.diff.dat.A.gonads["percent.CTCCCCCC",]$d,fry.pca.dat.A.body["percent.CTCCCCCC",]$PC2,"K",pos=1)
text(sex.diff.dat.A.gonads["percent.CCCCTCCC",]$d,fry.pca.dat.A.body["percent.CCCCTCCC",]$PC2,"L",pos=1)

text(sex.diff.dat.A.gonads["percent.CCCCTT",]$d,fry.pca.dat.A.body["percent.CCCCTT",]$PC2,"M",pos=3)
text(sex.diff.dat.A.gonads["percent.CCCCCT",]$d,fry.pca.dat.A.body["percent.CCCCCT",]$PC2,"N",pos=3)
text(sex.diff.dat.A.gonads["percent.CCTTTT",]$d,fry.pca.dat.A.body["percent.CCTTTT",]$PC2,"P",pos=3)

dev.print(device=pdf,file="Fig8-pc12-in-fry-vs-sex-diff-gonad-epitype-v2.pdf",height=5,width=8)
dev.print(device=jpeg,file="Fig8-pc12-in-fry-vs-sex-diff-gonad-epitype-v2.jpg",height=500,width=800)

###
#### Figure 9
###

rm(list=ls())

methyl.mat <- read.csv("pct-methylation.csv", na.strings=c("NA","Na","na","N/A","N/a","n/a","","-"),stringsAsFactors=FALSE)
attach(methyl.mat)
m.AI <- apply(cbind(AI130, AI135, AI143, AI150, AI163, AI181, AI187, AI228),1,mean,na.rm=TRUE)
m.AP <- apply(cbind(AP145, AP159, AP188, AP279, AP317, AP326),1,mean,na.rm=TRUE)
m.BI <- apply(cbind(BI217, BI333, BI358, BI370),1,mean,na.rm=TRUE)
m.B2 <- apply(cbind(B2174, B2187, B2213, B2238, B2253, B2302, B2314),1,mean,na.rm=TRUE)
methyl.mat <- cbind(methyl.mat,m.AI,m.AP,m.BI,m.B2)
rm(m.AI,m.AP,m.BI,m.B2)
detach(methyl.mat);attach(methyl.mat)

means.AI.adult <- c(
           mean(m.AI[tissue=="body" & sex=="F"],na.rm=T),
           mean(m.AI[tissue=="body" & sex=="Y"],na.rm=T),
           mean(m.AI[tissue=="body" & sex=="R"],na.rm=T),
           mean(m.AI[tissue=="head" & sex=="F"],na.rm=T),
           mean(m.AI[tissue=="head" & sex=="Y"],na.rm=T),
           mean(m.AI[tissue=="head" & sex=="R"],na.rm=T)
           )
means.AP.adult <- c(
           mean(m.AP[tissue=="body" & sex=="F"],na.rm=T),
           mean(m.AP[tissue=="body" & sex=="Y"],na.rm=T),
           mean(m.AP[tissue=="body" & sex=="R"],na.rm=T),
           mean(m.AP[tissue=="head" & sex=="F"],na.rm=T),
           mean(m.AP[tissue=="head" & sex=="Y"],na.rm=T),
           mean(m.AP[tissue=="head" & sex=="R"],na.rm=T)
           )

means.BI.adult<- c(
           mean(m.BI[tissue=="body" & sex=="F"],na.rm=T),
           mean(m.BI[tissue=="body" & sex=="Y"],na.rm=T),
           mean(m.BI[tissue=="body" & sex=="R"],na.rm=T),
           mean(m.BI[tissue=="head" & sex=="F"],na.rm=T),
           mean(m.BI[tissue=="head" & sex=="Y"],na.rm=T),
           mean(m.BI[tissue=="head" & sex=="R"],na.rm=T)
           )
means.B2.adult<- c(
           mean(m.B2[tissue=="body" & sex=="F"],na.rm=T),
           mean(m.B2[tissue=="body" & sex=="Y"],na.rm=T),
           mean(m.B2[tissue=="body" & sex=="R"],na.rm=T),
           mean(m.B2[tissue=="head" & sex=="F"],na.rm=T),
           mean(m.B2[tissue=="head" & sex=="Y"],na.rm=T),
           mean(m.B2[tissue=="head" & sex=="R"],na.rm=T)
           )

means.AI.fry <- c(
    mean(m.AI[tissue=="body" & sex=="Fry"],na.rm=T),
    mean(m.AI[tissue=="head" & sex=="Fry"],na.rm=T)
)
means.AP.fry<-c(
    mean(m.AP[tissue=="body" & sex=="Fry"],na.rm=T),
    mean(m.AP[tissue=="head" & sex=="Fry"],na.rm=T)
)
means.BI.fry <- c(
    mean(m.BI[tissue=="body" & sex=="Fry" ],na.rm=T),
    mean(m.BI[tissue=="head" & sex=="Fry" ],na.rm=T)
)
means.B2.fry <- c(          
    mean(m.B2[tissue=="body" & sex=="Fry" ],na.rm=T),
    mean(m.B2[tissue=="head" & sex=="Fry" ],na.rm=T)
)



expression.dat <- read.csv("expression_qPCR.csv")
detach(methyl.mat);attach(expression.dat)
meansExp.A.adult<- c(
    mean(rel_cyp19a1A[tissue=="body" & sex=="F"],na.rm=T),
    mean(rel_cyp19a1A[tissue=="body" & sex=="Y"],na.rm=T),
    mean(rel_cyp19a1A[tissue=="body" & sex=="R"],na.rm=T),
    mean(rel_cyp19a1A[tissue=="head" & sex=="F"],na.rm=T),
    mean(rel_cyp19a1A[tissue=="head" & sex=="Y"],na.rm=T),
    mean(rel_cyp19a1A[tissue=="head" & sex=="R"],na.rm=T)
)

meansExp.B.adult <- c(
    mean(rel_cyp19a1B[tissue=="body" & sex=="F"],na.rm=T),
    mean(rel_cyp19a1B[tissue=="body" & sex=="Y"],na.rm=T),
    mean(rel_cyp19a1B[tissue=="body" & sex=="R"],na.rm=T),
    mean(rel_cyp19a1B[tissue=="head" & sex=="F"],na.rm=T),
    mean(rel_cyp19a1B[tissue=="head" & sex=="Y"],na.rm=T),
    mean(rel_cyp19a1B[tissue=="head" & sex=="R"],na.rm=T)
)

meansExp.A.fry<- c(
    mean(rel_cyp19a1A[tissue=="body" & sex=="fry"],na.rm=T),
    mean(rel_cyp19a1A[tissue=="head" & sex=="fry"],na.rm=T)
)
meansExp.B.fry <- c(
    mean(rel_cyp19a1B[tissue=="body" & sex=="fry"],na.rm=T),
    mean(rel_cyp19a1B[tissue=="head" & sex=="fry"],na.rm=T)
)

par(mfrow=c(2,1))
par(oma=c(0,0,0,0))
par(mar=c(4,4,2,2))
plot(meansExp.A.adult, means.AI.adult, type = "n", ylab = "Percent methylation", xlab = "Relative cyp19a1A expression", ylim =c(0,100))
abline (65, -40, col = "grey60", lty = 2)
abline (120, -40, col = "grey60", lty = 2)
xx <- c(-0.2, -0.2, 0.25, 2.4, 2.4, 1.75, -0.2)
yy <- c( 73, 110, 110, 24, -5, -5, 73)
#xx<- c(0,0,0.5,2.2,2.2,1.65,0)  # fix numbers to fill
#yy <- c(65,100,100,31,0,0,65)  # fix numbers to fill
polygon(xx, yy, col="grey95", border = NA)
points(meansExp.A.adult, means.AI.adult , cex = 1.5, pch = c(19,19,19,15,15,15), col = c("black","red3","goldenrod1"), lwd = 2)
points(  meansExp.A.adult, means.AP.adult , cex = 1.5, pch = c(1,1,1,0,0,0),  col = c("black","red3","goldenrod1"), lwd = 2)
points ( meansExp.A.fry*10, means.AI.fry, pch = c(19,15), col = "grey60", lwd = 2) # note fry expression an order of magnitude lower, multiply by 10 to put on same plot
points ( meansExp.A.fry*10, means.AP.fry, pch = c(1,0), col = "grey60", lwd = 2)
#text(1.5,100,"filled = AI")
#text(1.5,90,"open = AP")
#text(1.5,80,"square = brain")
#text(1.5,70, "circle =gonad")
mtext("A)",side=3,line=0,at=0)

## make fonts and axis labels match the other figures
plot( meansExp.B.adult, means.BI.adult , type = "n", ylab = "Percent methylation", xlab = "Relative cyp19a1B expression", ylim =c(0,100))
points(meansExp.B.adult, means.BI.adult , cex = 1.5, pch = c(19,19,19,15,15,15), col = c("black","red3","goldenrod1"), lwd = 2)
points(  meansExp.B.adult, means.B2.adult , cex = 1.5, pch = c(1,1,1,0,0,0),  col = c("black","red3","goldenrod1"), lwd = 2)
points ( meansExp.B.fry*10, means.BI.fry, pch = c(19,15), col = "grey60", lwd = 2) # note fry expression an order of magnitude lower, multiply by 10 to put on same plot
points ( meansExp.B.fry*10, means.B2.fry, pch = c(1,0), col = "grey60", lwd = 2)
#text(1.5,100,"filled = BI")
#text(1.5,90,"open = BP")
#text(1.5,80,"square = brain")
#text(1.5,70, "circle =gonad")
mtext("B)",side=3,line=0,at=0)

dev.print(device=pdf,file="Fig9-ExpByMeth.pdf",height=8, width=6)
dev.print(device=jpeg,file="Fig9-ExpByMeth.jpg",height=600, width=400)

