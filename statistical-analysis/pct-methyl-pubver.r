### Code for Tables 4,5 & 6, and Tables 4 & 5 post-hocs in text
rm(list=ls())
library(car)

methyl.mat <- read.csv("pct-methylation.csv", na.strings=c("NA","Na","na","N/A","N/a","n/a","","-"),stringsAsFactors=FALSE)

attach(methyl.mat)
### calculate mean methylation for each amplicon
m.AI <- apply(cbind(AI130, AI135, AI143, AI150, AI163, AI181, AI187, AI228),1,mean,na.rm=TRUE)
m.AP <- apply(cbind(AP145, AP159, AP188, AP279, AP317, AP326),1,mean,na.rm=TRUE)
m.BI <- apply(cbind(BI217, BI333, BI358, BI370),1,mean,na.rm=TRUE)
m.B2 <- apply(cbind(B2174, B2187, B2213, B2238, B2253, B2302, B2314),1,mean,na.rm=TRUE)
methyl.mat <- cbind(methyl.mat,m.AI,m.AP,m.BI,m.B2)
rm(m.AI,m.AP,m.BI,m.B2)
detach(methyl.mat)
    
methyl.mat.ad <- subset(methyl.mat,sex!="Fry")
methyl.mat.fry <- subset(methyl.mat,sex=="Fry")

AI.mod <- aov(m.AI~tissue*sex,data=methyl.mat.ad)
Anova(AI.mod)
#Anova Table (Type II tests)
#
#Response: m.AI
#           Sum Sq Df F value    Pr(>F)    
#tissue     6643.9  1 212.635 2.006e-13 ***
#sex        4077.7  3  43.501 7.360e-10 ***
#tissue:sex 4596.0  3  49.031 2.182e-10 ***
#Residuals   749.9 24                      
TukeyHSD(AI.mod,which="tissue:sex")
#  Tukey multiple comparisons of means
#    95% family-wise confidence level
#
#Fit: aov(formula = m.AI ~ tissue * sex, data = methyl.mat.ad)
#
#$`tissue:sex`
#                    diff         lwr        upr     p adj
#head:B-body:B -49.348310 -63.4877232 -35.208897 0.0000000
#body:F-body:B -60.221592 -74.3610054 -46.082179 0.0000000
#head:F-body:B -48.218270 -62.3576830 -34.078857 0.0000000
#body:R-body:B  -8.120018 -22.2594311   6.019395 0.5625256
#head:R-body:B -46.573323 -60.7127356 -32.433910 0.0000000
#body:Y-body:B -10.107896 -24.2473088   4.031517 0.3005918
#head:Y-body:B -51.567095 -65.0869596 -38.047230 0.0000000
#body:F-head:B -10.873282 -23.9638350   2.217271 0.1553843
#head:F-head:B   1.130040 -11.9605126  14.220593 0.9999896
#body:R-head:B  41.228292  28.1377393  54.318845 0.0000000
#head:R-head:B   2.774988 -10.3155652  15.865540 0.9961246
#body:Y-head:B  39.240414  26.1498616  52.330967 0.0000000
#head:Y-head:B  -2.218785 -14.6375734  10.200004 0.9986716
#head:F-body:F  12.003322  -1.0872304  25.093875 0.0890070
#body:R-body:F  52.101574  39.0110215  65.192127 0.0000000
#head:R-body:F  13.648270   0.5577171  26.738823 0.0367330
#body:Y-body:F  50.113697  37.0231438  63.204249 0.0000000
#head:Y-body:F   8.654498  -3.7642912  21.073286 0.3295084
#body:R-head:F  40.098252  27.0076991  53.188805 0.0000000
#head:R-head:F   1.644947 -11.4456053  14.735500 0.9998672
#body:Y-head:F  38.110374  25.0198215  51.200927 0.0000000
#head:Y-head:F  -3.348825 -15.7676136   9.069964 0.9839155
#head:R-body:R -38.453304 -51.5438572 -25.362752 0.0000000
#body:Y-body:R  -1.987878 -15.0784304  11.102675 0.9995354
#head:Y-body:R -43.447077 -55.8658655 -31.028288 0.0000000
#body:Y-head:R  36.465427  23.3748740  49.555980 0.0000001
#head:Y-head:R  -4.993772 -17.4125610   7.425017 0.8776564
#head:Y-body:Y -41.459199 -53.8779878 -29.040410 0.0000000

AP.mod <- aov(m.AP~tissue*sex,data=methyl.mat.ad)
Anova(AP.mod)
#Anova Table (Type II tests)
#
#Response: m.AP
#           Sum Sq Df F value    Pr(>F)    
#tissue      57.82  1  2.8418   0.10480    
#sex        427.82  3  7.0087   0.00151 ** 
#tissue:sex 767.46  3 12.5728 3.856e-05 ***
#Residuals  488.33 24                      
TukeyHSD(AP.mod,which="tissue:sex")
#Fit: aov(formula = m.AP ~ tissue * sex, data = methyl.mat.ad)
#
#$`tissue:sex`
#                     diff         lwr        upr     p adj
#head:B-body:B  -6.8655682 -18.2756302  4.5444939 0.5065020
#body:F-body:B -18.9784881 -30.3885502 -7.5684261 0.0002684
#head:F-body:B  -4.9977755 -16.4078376  6.4122866 0.8247509
#body:R-body:B   2.8269906  -8.5830715 14.2370527 0.9901221
#head:R-body:B  -8.2860780 -19.6961400  3.1239841 0.2832361
#body:Y-body:B  -0.7485578 -12.1586198 10.6615043 0.9999984
#head:Y-body:B  -7.6730191 -18.5831252  3.2370870 0.3190170
#body:F-head:B -12.1129199 -22.6765847 -1.5492551 0.0168220
#head:F-head:B   1.8677927  -8.6958721 12.4314575 0.9987570
#body:R-head:B   9.6925588  -0.8711060 20.2562236 0.0886533
#head:R-head:B  -1.4205098 -11.9841746  9.1431550 0.9997915
#body:Y-head:B   6.1170104  -4.4466544 16.6806752 0.5526951
#head:Y-head:B  -0.8074509 -10.8290233  9.2141214 0.9999935
#head:F-body:F  13.9807126   3.4170478 24.5443774 0.0041946
#body:R-body:F  21.8054787  11.2418139 32.3691435 0.0000112
#head:R-body:F  10.6924102   0.1287454 21.2560750 0.0458145
#body:Y-body:F  18.2299303   7.6662655 28.7935951 0.0001622
#head:Y-body:F  11.3054690   1.2838967 21.3270414 0.0193826
#body:R-head:F   7.8247661  -2.7388987 18.3884309 0.2620058
#head:R-head:F  -3.2883025 -13.8519673  7.2753623 0.9646942
#body:Y-head:F   4.2492177  -6.3144471 14.8128825 0.8774775
#head:Y-head:F  -2.6752436 -12.6968159  7.3463288 0.9848110
#head:R-body:R -11.1130686 -21.6767334 -0.5494038 0.0342821
#body:Y-body:R  -3.5755484 -14.1392132  6.9881164 0.9456539
#head:Y-body:R -10.5000097 -20.5215821 -0.4784373 0.0353749
#body:Y-head:R   7.5375202  -3.0261446 18.1011850 0.3026901
#head:Y-head:R   0.6130589  -9.4085135 10.6346312 0.9999990
#head:Y-body:Y  -6.9244613 -16.9460337  3.0971110 0.3394057

BI.mod <- aov(m.BI~tissue*sex,data=methyl.mat.ad)
Anova(BI.mod)
#Anova Table (Type II tests)
#
#Response: m.BI
#            Sum Sq Df F value    Pr(>F)    
#tissue      7257.0  1 17.0837 0.0003761 ***
#sex         2347.5  3  1.8420 0.1665006    
#tissue:sex  1792.1  3  1.4062 0.2652811    
#Residuals  10195.0 24                      
TukeyHSD(BI.mod,which="tissue:sex")
#  Tukey multiple comparisons of means
#    95% family-wise confidence level
#
#Fit: aov(formula = m.BI ~ tissue * sex, data = methyl.mat.ad)
#
#$`tissue:sex`
#                    diff          lwr         upr     p adj
#head:B-body:B -53.074320 -105.2089167  -0.9397228 0.0439257
#body:F-body:B -28.287806  -80.4224026  23.8467913 0.6278377
#head:F-body:B -40.771774  -92.9063707  11.3628231 0.2077579
#body:R-body:B   4.788443  -47.3461539  56.9230399 0.9999841
#head:R-body:B -32.932960  -85.0675570  19.2016368 0.4471131
#body:Y-body:B  -5.372836  -57.5074333  46.7617605 0.9999651
#head:Y-body:B -27.240754  -77.0909632  22.6094557 0.6199044
#body:F-head:B  24.786514  -23.4807436  73.0537718 0.6867819
#head:F-head:B  12.302546  -35.9647118  60.5698037 0.9883511
#body:R-head:B  57.862763    9.5955050 106.1300205 0.0112409
#head:R-head:B  20.141360  -28.1258981  68.4086174 0.8565049
#body:Y-head:B  47.701483   -0.5657744  95.9687411 0.0543518
#head:Y-head:B  25.833566  -19.9567752  71.6239073 0.5833734
#head:F-body:F -12.483968  -60.7512259  35.7832896 0.9873237
#body:R-body:F  33.076249  -15.1910091  81.3435064 0.3490502
#head:R-body:F  -4.645154  -52.9124121  43.6221033 0.9999781
#body:Y-body:F  22.914969  -25.3522884  71.1822270 0.7613106
#head:Y-body:F   1.047052  -44.7432893  46.8373932 1.0000000
#body:R-head:F  45.560217   -2.7070409  93.8274745 0.0740952
#head:R-head:F   7.838814  -40.4284440  56.1060715 0.9992806
#body:Y-head:F  35.398937  -12.8683203  83.6661951 0.2725851
#head:Y-head:F  13.531020  -32.2593212  59.3213613 0.9732793
#head:R-body:R -37.721403  -85.9886608  10.5458547 0.2084046
#body:Y-body:R -10.161279  -58.4285371  38.1059784 0.9962860
#head:Y-body:R -32.029197  -77.8195380  13.7611445 0.3252464
#body:Y-head:R  27.560124  -20.7071340  75.8273814 0.5693085
#head:Y-head:R   5.692206  -40.0981349  51.4825476 0.9998765
#head:Y-body:Y -21.867917  -67.6582586  23.9224239 0.7561258

BP.mod <- aov(m.B2~tissue*sex,data=methyl.mat.ad)
Anova(BP.mod)
#Anova Table (Type II tests)
#
#Response: m.B2
#           Sum Sq Df F value  Pr(>F)  
#tissue      313.6  1  2.1830 0.15255  
#sex         878.1  3  2.0377 0.13538  
#tissue:sex 1163.9  3  2.7007 0.06816 .
#Residuals  3447.6 24                  
TukeyHSD(BP.mod,which="tissue:sex")
#  Tukey multiple comparisons of means
#    95% family-wise confidence level
#
#Fit: aov(formula = m.B2 ~ tissue * sex, data = methyl.mat.ad)
#
#$`tissue:sex`
#                     diff         lwr       upr     p adj
#head:B-body:B -10.8929593 -41.2103641 19.424446 0.9271203
#body:F-body:B -22.1891386 -52.5065435  8.128266 0.2747855
#head:F-body:B  -7.9593078 -38.2767127 22.358097 0.9861896
#body:R-body:B   4.3642151 -25.9531898 34.681620 0.9996729
#head:R-body:B  -7.1677713 -37.4851762 23.149634 0.9925065
#body:Y-body:B   5.2708681 -25.0465368 35.588273 0.9988841
#head:Y-body:B -11.2035155 -40.1924992 17.785468 0.8974993
#body:F-head:B -11.2961794 -39.3646422 16.772283 0.8772072
#head:F-head:B   2.9336514 -25.1348114 31.002114 0.9999616
#body:R-head:B  15.2571743 -12.8112885 43.325637 0.6258289
#head:R-head:B   3.7251879 -24.3432749 31.793651 0.9998089
#body:Y-head:B  16.1638273 -11.9046355 44.232290 0.5592549
#head:Y-head:B  -0.3105562 -26.9386381 26.317526 1.0000000
#head:F-body:F  14.2298308 -13.8386320 42.298294 0.6997936
#body:R-body:F  26.5533537  -1.5151091 54.621817 0.0730300
#head:R-body:F  15.0213673 -13.0470955 43.089830 0.6430422
#body:Y-body:F  27.4600067  -0.6084561 55.528470 0.0583150
#head:Y-body:F  10.9856231 -15.6424587 37.613705 0.8632950
#body:R-head:F  12.3235229 -15.7449399 40.391986 0.8230790
#head:R-head:F   0.7915365 -27.2769263 28.859999 1.0000000
#body:Y-head:F  13.2301759 -14.8382869 41.298639 0.7675154
#head:Y-head:F  -3.2442076 -29.8722895 23.383874 0.9998921
#head:R-body:R -11.5319864 -39.6004492 16.536476 0.8657002
#body:Y-body:R   0.9066530 -27.1618098 28.975116 1.0000000
#head:Y-body:R -15.5677305 -42.1958124 11.060351 0.5412566
#body:Y-head:R  12.4386394 -15.6298234 40.507102 0.8163971
#head:Y-head:R  -4.0357441 -30.6638260 22.592338 0.9995413
#head:Y-body:Y -16.4743836 -43.1024655 10.153698 0.4725682

###
#### Fry pct-methylated
###
###

summary(aov(m.AI~tissue+Error(set),data=methyl.mat.fry))
#Error: set
#          Df Sum Sq Mean Sq F value Pr(>F)
#Residuals  2  221.5   110.7               
#
#Error: Within
#          Df Sum Sq Mean Sq F value   Pr(>F)    
#tissue     1  680.4   680.4   149.8 2.72e-11 ***
#Residuals 22   99.9     4.5                     

summary(aov(m.AP~tissue+Error(set),data=methyl.mat.fry))
#Error: set
#          Df Sum Sq Mean Sq F value Pr(>F)
#Residuals  2  193.9   96.95               
#
#Error: Within
#          Df Sum Sq Mean Sq F value Pr(>F)
#tissue     1   36.4   36.37   2.463  0.131
#Residuals 22  324.8   14.76               

summary(aov(m.BI~tissue+Error(set),data=methyl.mat.fry))
#Error: set
#          Df Sum Sq Mean Sq F value Pr(>F)
#Residuals  2   1460   729.9               
#
#Error: Within
#          Df Sum Sq Mean Sq F value Pr(>F)
#tissue     1    338   338.2   1.225   0.28
#Residuals 22   6073   276.0               

summary(aov(m.B2~tissue+Error(set),data=methyl.mat.fry))
#Error: set
#          Df Sum Sq Mean Sq F value Pr(>F)
#Residuals  2  378.8   189.4               
#
#Error: Within
#          Df Sum Sq Mean Sq F value Pr(>F)
#tissue     1  207.6  207.63    2.71  0.114
#Residuals 22 1685.5   76.61               
