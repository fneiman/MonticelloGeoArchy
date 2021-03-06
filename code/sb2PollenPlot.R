
library(rioja)
library(dplyr)
library(mgcv)


sb2 <- read.csv('data/SlateBranchUnit2Pollen.csv', stringsAsFactors = F)
# get rid of row of NAs
sb2 <- sb2[1:43,]
# fix Depth
sb2$Depth <- -sb2$Depth
# sort on depth






# check to see what it looks like
str(sb2)
# what vars are in what rows?
cbind(1:ncol(sb2), colnames(sb2))



# leave out "unknown x", other, indeterminant, and sums
taxonCounts <- sb2[,c(5:87,90)]
# get rid of taxa with all NAs
taxonCounts <- taxonCounts[,!is.na(colSums(taxonCounts))]
# get rid of taxa with all 0s 
taxonCounts <- taxonCounts[, colSums(taxonCounts) > 0] 

# find the sample with no pollen and get rid of them
rowSums(taxonCounts) 
sb2<- sb2[(rowSums(taxonCounts) != 0),]
taxonCounts<- taxonCounts[(rowSums(taxonCounts) != 0),]



# histgram of number of occurrences
hist(colSums(taxonCounts > 0), breaks=1:50, col='grey') 

# get rid of taxa with fewer than 5 occurrences
taxonCounts <- taxonCounts[,colSums(taxonCounts > 0) > 5]




# create appropriately sized graphics window
#windows(width=16, height=10)

taxonPcts <- prop.table(as.matrix(taxonCounts), 1)*100
rowSums(taxonPcts)

elevation <- sb2$Elevation
depth <- sb2$Depth 


# add a dendrogram from constrained cluster analysis
diss <- dist(taxonPcts, method="euclidean" )
clust <- chclust(diss, method="coniss")
# broken stick model suggest 6 significant zones
bstick(clust)

colnames(taxonPcts)

title= 'Slate Branch: Unit 2'
xx <- 1:ncol(taxonPcts)
tColors <- ifelse(xx < 24, "darkgreen", "darkred")
x <- strat.plot(taxonPcts, 
        
           yvar = depth, 
           y.rev=TRUE, 
           #y.tks= depth,
           ylabel="Depth (feet)", 
           cex.ylabel = 1.5, 
           cex.yaxis = 1,
           srt.xlabel = 90,
           cex.xlabel= 1.25,
           
           scale.percent=TRUE, 
         
           plot.poly=TRUE, 
           col.poly= tColors, 
           col.poly.line=NA,
           exag=TRUE,
           exag.mult =10,
           col.exag="auto",
           clust=clust)

addClustZone(x, clust, 6, col="black")

par(mfrow=c(1,1))
OakHickoryChestnut <- sb2$Quercus + sb2$Carya + sb2$Castanea 

Asters <- sb2$Asteraceae.LS +  sb2$Asteraceae.HS + sb2$Asteraceae.Liguliflorae 
asterIndex <- data.frame(Asters, OakHickoryChestnut, Depth = -depth)
asterModel <- gam( cbind(Asters,OakHickoryChestnut)~ s(Depth),
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
with(asterIndex, plot(Depth, Asters/(Asters + OakHickoryChestnut),
                      pch=21, bg='green', cex=2,
                      cex.lab= 2,
                      cex.axis=1.5,
                      xlim= c(-5,0)
                       ))

lines(asterIndex$Depth, pHat, lwd=3, col= 'darkgreen')
lines(asterIndex$Depth, upperCL, lwd=1, col= 'darkgreen', lty=2)
lines(asterIndex$Depth, lowerCL, lwd=1, col= 'darkgreen', lty=2)
#abline( v=c(-.5, -2.5, -3.1), lty=1, lwd=2
        , col='grey')



Pine <-  sb2$Pinus
pineIndex <- data.frame(Pine, OakHickoryChestnut, Depth = -depth)
pineModel <- gam( cbind(Pine,OakHickoryChestnut)~ s(Depth),
                   link = logit,
                   family = quasibinomial,
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
with(pineIndex, plot(Depth, Pine/(Pine + OakHickoryChestnut),
                      pch=21, bg='green', cex=2,
                      cex.lab= 2,
                      cex.axis=1.5,
                      xlim = c(-5, 0)
))

lines(pineIndex$Depth, pHat, lwd=3, col= 'darkgreen')
lines(pineIndex$Depth, upperCL, lwd=1, col= 'darkgreen', lty=2)
lines(pineIndex$Depth, lowerCL, lwd=1, col= 'darkgreen', lty=2)
#abline( v=c(-.5, -2.5, -3.1), lty=1, lwd=2, col='grey')




#####################################################

ChenoAm <-  sb2$Cheno.Am
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
#abline( v=c(-.5, -2.5, -3.1), lty=1, lwd=2, col='grey')


par(mfrow=c(1,1))

######################################################
Trifolium <-  sb2$Trifolium + sb2$Trifolium.repens
trifoliumIndex <- data.frame(Trifolium, Asters, Depth = -depth)
trifoliumModel <- gam( cbind(Trifolium,Asters) ~ s(Depth),
                     link = logit,
                     family = binomial,
                    # method= 'REML',
                     #correlation = corAR1(form = Depth),
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
#abline( v=c(-.5, -2.5, -3.1), lty=1, lwd=2, col='grey')


                                   
