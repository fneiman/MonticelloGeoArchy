
library(rioja)
library(dplyr)
library(mgcv)


pp1 <- read.csv('https://raw.github.com/fneiman/MonticelloGeoArchy/master/data/pawPawUnit1Pollen.csv', 
                stringsAsFactors = F)


# check to see what it looks like
str(pp1)
# what vars are in what rows?
cbind(1:ncol(pp1), colnames(pp1))



# leave out "unknown x", other, indeterminant, and sums
taxonCounts <- pp1[,c(9:116,122)]
# get rid of taxa with all NAs
taxonCounts <- taxonCounts[,!is.na(colSums(taxonCounts))]
# get rid of taxa with all 0s 
taxonCounts <- taxonCounts[, colSums(taxonCounts) > 0] 

# histogram of number of occurrences
hist(colSums(taxonCounts > 0), breaks=1:50, col='grey') 

# get rid of taxa with fewer than 5 occurrences
taxonCounts <- taxonCounts[,colSums(taxonCounts > 0) > 5]


# write a function to subset pollen taxa based on minimum percent frequency 
# and a minimum number of samples in which that frequnecy is attained
moreThan2px3 <- function(taxonCounts, minPercent, minNumSamples ){
  taxonPcts <- prop.table(as.matrix(taxonCounts), 1)*100
  gt2p <- colSums(taxonPcts >=2) >=3
  return(taxonCounts[,gt2p])
}


taxonCounts1 <- moreThan2px3(taxonCounts, 2, 1)



taxonPcts <- prop.table(as.matrix(taxonCounts1), 1)*100
rowSums(taxonPcts)

elevation <- pp1$Elevation
depth <- pp1$Depth 


# add a dendrogram from constrained cluster analysis
diss <- dist(taxonPcts, method="euclidean" )
clust <- chclust(diss, method="coniss")
# broken stick model suggest 3 significant zones
bstick(clust)

colnames(taxonPcts)

title= 'Paw Paw: Unit 1'
xx <- 1:ncol(taxonPcts)
tColors <- ifelse(xx < 7, "darkgreen", "darkred")
x <- strat.plot(taxonPcts, 
        
           yvar = depth, 
           y.rev=TRUE, 
           #y.tks= depth,
           ylabel="Depth (feet)", 
           cex.ylabel = 1.5, 
           cex.yaxis = 1,
           srt.xlabel = 90,
           cex.xlabel= 1.5,
           
           scale.percent=TRUE, 
         
           plot.poly=TRUE, 
           col.poly= tColors, 
           col.poly.line=NA,
           exag=TRUE,
           exag.mult =10,
           col.exag="auto",
           clust=clust)

addClustZone(x, clust, 5, col="black")


# pull out individual taxa a fit quasi binomial GAMs

# asteraceae
Oak.Hickory.Chestnut <- pp1$Quercus + pp1$Carya + pp1$Castanea 

Asters <- pp1$Asteraceae.LS +  pp1$Asteraceae.HS 
asterIndex <- data.frame(Asters, Oak.Hickory.Chestnut, Depth = -depth)
asterModel <- gam( cbind(Asters,Oak.Hickory.Chestnut)~ s(Depth),
                  link = logit,
                  family = quasibinomial,
                  data= asterIndex)
summary(asterModel)
anova(asterModel)
acf(residuals(asterModel),ci.type = "ma")

invLogit <- function(x){ exp(x)/(1+exp(x)) }
pred <- predict(asterModel, type='link', se.fit=T) 
pHat <- invLogit(pred$fit)
upperCL <- invLogit(pred$fit + 2*pred$se.fit)
lowerCL <- invLogit(pred$fit - 2*pred$se.fit)

par(mar=c(6,6,2,2))
with(asterIndex, plot(Depth, Asters/(Asters + Oak.Hickory.Chestnut),
                      pch=21, bg='green', cex=2,
                      cex.lab= 2,
                      cex.axis=1.5,
                      xlim= c(-5,0)
                       ))

lines(asterIndex$Depth, pHat, lwd=3, col= 'darkgreen')
lines(asterIndex$Depth, upperCL, lwd=1, col= 'darkgreen', lty=2)
lines(asterIndex$Depth, lowerCL, lwd=1, col= 'darkgreen', lty=2)
abline( v=c(-.5, -2.5, -3.1), lty=1, lwd=2
        , col='grey')


# Pine

Pine <-  pp1$Pinus
pineIndex <- data.frame(Pine, Oak.Hickory.Chestnut, Depth = -depth)
pineModel <- gam( cbind(Pine,Oak.Hickory.Chestnut)~ s(Depth),
                   link = logit,
                   family = binomial,
                   method='REML',
                   #correlation = corAR1(form =  ~ Depth),
                   data= pineIndex)

summary(pineModel)
anova(pineModel)
plot(pineModel)


acf(residuals(pineModel),ci.type = "ma")

pred <- predict(pineModel, type='link', se.fit=T) 
pHat <- invLogit(pred$fit)
upperCL <- invLogit(pred$fit + 2*pred$se.fit)
lowerCL <- invLogit(pred$fit - 2*pred$se.fit)

par(mar=c(6,6,2,2))
with(pineIndex, plot(Depth, Pine/(Pine + Oak.Hickory.Chestnut),
                      pch=21, bg='green', cex=2,
                      cex.lab= 2,
                      cex.axis=1.5,
                      xlim = c(-5, 0)
))

lines(pineIndex$Depth, pHat, lwd=3, col= 'darkgreen')
lines(pineIndex$Depth, upperCL, lwd=1, col= 'darkgreen', lty=2)
lines(pineIndex$Depth, lowerCL, lwd=1, col= 'darkgreen', lty=2)
abline( v=c(-.5, -2.5, -3.1), lty=1, lwd=2, col='grey')


# Hickory Chestnut
Hickory.Chestnut <-pp1$Carya + pp1$Castanea 

HCIndex <- data.frame(Hickory.Chestnut, Oak.Hickory.Chestnut, Depth = -depth)
HCModel <- gam( cbind(Hickory.Chestnut, Oak.Hickory.Chestnut)~ s(Depth),
                   link = logit,
                   family = quasibinomial,
                   data= HCIndex)
summary(HCModel)
anova(HCModel)
plot(HCModel)
acf(residuals(HCModel),ci.type = "ma")

invLogit <- function(x){ exp(x)/(1+exp(x)) }
pred <- predict(HCModel, type='link', se.fit=T) 
pHat <- invLogit(pred$fit)
upperCL <- invLogit(pred$fit + 2*pred$se.fit)
lowerCL <- invLogit(pred$fit - 2*pred$se.fit)

par(mar=c(6,6,2,2))
with(HCIndex, plot(Depth, Hickory.Chestnut/(Hickory.Chestnut + Oak.Hickory.Chestnut),
                      pch=21, bg='green', cex=2,
                      cex.lab= 2,
                      cex.axis=1.5,
                      xlim= c(-5,0)
))

lines(HCIndex$Depth, pHat, lwd=3, col= 'darkgreen')
lines(HCIndex$Depth, upperCL, lwd=1, col= 'darkgreen', lty=2)
lines(HCIndex$Depth, lowerCL, lwd=1, col= 'darkgreen', lty=2)
abline( v=c(-.5, -2.5, -3.1), lty=1, lwd=2
        , col='grey')



# ChenoAms


ChenoAm <-  pp1$ChenoAm
chenoAmIndex <- data.frame(ChenoAm, Asters, Depth = -depth)
chenoAmModel <- gam( cbind(ChenoAm,Asters)~ s(Depth),
                  link = logit,
                  family = quasibinomial,
                  #correlation = corAR1(form =  Depth),
                  data= chenoAmIndex)

summary(chenoAmModel)
anova(chenoAmModel)
plot(chenoAmModel)

acf(residuals(chenoAmModel),ci.type = "ma")

pred <- predict(chenoAmModel, type='link', se.fit=T) 
pHat <- invLogit(pred$fit)
upperCL <- invLogit(pred$fit + 2*pred$se.fit)
lowerCL <- invLogit(pred$fit - 2*pred$se.fit)

par(mar=c(6,6,2,2))
with(chenoAmIndex, plot(Depth, ChenoAm/(ChenoAm + Asters),
                     pch=21, bg='green', cex=2,
                     cex.lab= 2,
                     ylim= c(0.,.4),
                     cex.axis=1.5,
                     xlim = c(-5, 0)))

lines(chenoAmIndex$Depth, pHat, lwd=3, col= 'darkgreen')
lines(chenoAmIndex$Depth, upperCL, lwd=1, col= 'darkgreen', lty=2)
lines(chenoAmIndex$Depth, lowerCL, lwd=1, col= 'darkgreen', lty=2)
abline( v=c(-.5, -2.5, -3.1), lty=1, lwd=2, col='grey')


par(mfrow=c(1,1))

# trifolum
######################################################
Trifolium <-  pp1$Trifolium + pp1$Trifolium.repens
trifoliumIndex <- data.frame(Trifolium, Asters, Depth = -depth)
trifoliumModel <- gam( cbind(Trifolium,Asters) ~ s(Depth),
                     family = quasibinomial (link=logit),
                     method= 'REML',
                     data= trifoliumIndex)

summary(trifoliumModel)
anova(trifoliumModel)
plot(trifoliumModel)


acf(residuals(trifoliumModel),ci.type = "ma")


invLogit <- function(x){ exp(x)/(1+exp(x)) }

pred <- predict(trifoliumModel, type='link', se.fit=T) 
pHat <- invLogit(pred$fit)
upperCL <- invLogit(pred$fit + 2*pred$se.fit)
lowerCL <- invLogit(pred$fit - 2*pred$se.fit)

par(mar=c(6,6,2,2))
with(trifoliumIndex, plot(Depth, Trifolium/(Trifolium + Asters),
                        pch=21, bg='green', cex=2,
                        cex.lab= 2,
                        cex.axis=1.5,
                        xlim=c(-5,0)))


lines(trifoliumIndex$Depth, pHat, lwd=3, col= 'darkgreen')
lines(trifoliumIndex$Depth, upperCL, lwd=1, col= 'darkgreen', lty=2)
lines(trifoliumIndex$Depth, lowerCL, lwd=1, col= 'darkgreen', lty=2)
abline( v=c(-.5, -2.5, -3.1), lty=1, lwd=2, col='grey')


                                   
