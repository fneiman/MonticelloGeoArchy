
library(rioja)
library(dplyr)
library(viridisLite)
library(maptools)


sb2Chem <- read.csv('data/slateBranch2-Wisc-Chem.csv', stringsAsFactors = F)


# check to see what it looks like
str(sb2Chem)
# what vars are in what rows?
cbind(1:ncol(sb2Chem), colnames(sb2Chem))


sb2Chem <- arrange(sb2Chem, Depth)

lnSB2Chem <- log(sb2Chem[,10:22])
                 
                 
elevation <- sb2Chem$Elevation
depth <- sb2Chem$Depth 


# add a dendrogram from constrained cluster analysis
diss <- dist(lnSB2Chem, method="euclidean" )
clust <- chclust(diss, method="coniss")
# broken stick model suggest 3 significant zones
bstick(clust)

title= 'Slate Branch: Unit 2'


x <- strat.plot(lnSB2Chem, 
           yvar = depth, 
           y.rev= T, 
           #y.tks= depth,
           ylabel="Depth (feet)", 
           cex.ylabel = 1.5, 
           cex.yaxis = 1,
           srt.xlabel = 45,
           cex.xlabel= 1.5,
           scale.percent=F, 
           
           plot.poly=TRUE, 
           col.poly= viridis(1), 
           col.poly.line=NA,
           xSpace=0.01, x.pc.lab=TRUE, x.pc.omit0=TRUE, las=2,
           clust=clust)

#lineCol <- adjustcolor(col="", alpha=.75)


addClustZone(x, clust,6, col="grey",
             lwd=2, lty=2)

 



# Classical PCA  #####

# first we check out the eigenvalues in a scree plot
# define a function for the broken stick model  
broken.stick <- function(p)
  # Compute the expected values of the broken-stick distribution for 'p' pieces.
  # Example: broken.stick.out.20 = broken.stick(20)
  #             Pierre Legendre, April 2007
{
  result = matrix(0,p,2)
  colnames(result) = c("j","E(j)")
  for(j in 1:p) {
    E = 0
    for(x in j:p) E = E+(1/x)
    result[j,1] = j
    result[j,2] = E/p
  }
  return(result)
}


pc1<-prcomp(~ Al+ Ba + Ca+ Fe + Mg +K + Mn + Na + P + Sr + Ti + Zn, 
            data=lnSB2Chem, retx = TRUE, center = TRUE, scale. = T)
# if scale.=TRUE is used, the PCs are extracted from the correlation matrix,
# vs. the covariance matrix 
#get a summary
summary(pc1)


# first we check out the eigenvalues in a scree plot
ProportionVariance<- pc1$sd^2/(sum(pc1$sd^2))

par(mfrow=c(1,1))
bs<-broken.stick(length(pc1$sd))
barplot(ProportionVariance, names.arg=1:length(pc1$sd),
        lwd=1, xlab="PC Order",
        pch=21, bg="light grey", cex=1.5,
        ylab= "Proportion of Variance", cex.axis=1.5, cex.lab=1.5)
lines(bs[,1],bs[,2], col=viridis(1), lwd=5, lty=1)


# plot the scores of the obs -- note we scale these to unit variance
# so scatterplot distances approximate Mahalaobis D
scores<-(predict(pc1))
pcaZScores <- apply(scores, 2, scale)


par(mfrow=c(1,2))
plot(pcaZScores[,1],pcaZScores[,2], 
     pch= 21, 
     bg= adjustcolor(viridis(14)[sb2Chem$Layer], alpha=.8),
     col="black", 
     #xlim=c(-2,4),
     cex= 3, cex.lab= 1.5, 
     xlab="PC 1 (40%)", 
     ylab="PC 2 (33%) ",
     asp=1)
abline(h=0,v=0,col="black", lwd=1, lty=2)
pointLabel(pcaZScores[,1],pcaZScores[,2] , as.character(sb2Chem$Layer),
           col="black")

layers <- unique(sb2Chem$Layer)
legend('topleft', pch = 21, 
       pt.bg= adjustcolor(viridis(length(layers)), alpha=.8),
       legend= layers,
       ncol=1,
       cex= 1.5,
       pt.cex = 2,
       bty='o')
#plot the variables on PC1 and PC2
plot(pc1$rotation[,1],pc1$rotation[,2], asp=1,type="n",
     #xlim= c(-.3,.3), ylim=c(-.3,.3),
     cex=2, cex.lab= 1.5, xlab="PC 1", ylab='PC 2')
arrows(0, 0, pc1$rotation[,1],pc1$rotation[,2], len=.2, lwd=2, col="grey")
abline(h=0,v=0,col="black", lwd=1, lty=2)
pointLabel(1.05* pc1$rotation[,1] , 1.05* pc1$rotation[,2] , names(pc1$rotation[,1]),
           col="black", cex=1.5)




par(mfrow=c(1,2))
plot(pcaZScores[,1],pcaZScores[,3], 
     pch= 21, 
     bg= adjustcolor(viridis(14)[sb2Chem$Layer], alpha=.8),
     col="black", 
     #xlim=c(-2,4),
     cex= 3, cex.lab= 1.5, 
     xlab="PC 1 (40%)", 
     ylab="PC 3 (13%) ",
     asp=1)
abline(h=0,v=0,col="black", lwd=1, lty=2)
pointLabel(pcaZScores[,1],pcaZScores[,3] , as.character(sb2Chem$Layer),
           col="black")

layers <- unique(sb2Chem$Layer)
legend('topleft', pch = 21, 
       pt.bg= adjustcolor(viridis(length(layers)), alpha=.8),
       legend= layers,
       ncol=1,
       cex= 1.5,
       pt.cex = 2,
       bty='o')


#plot the variables on PC1 and PC3
plot(pc1$rotation[,1],pc1$rotation[,3], asp=1,type="n",
     #xlim= c(-.3,.3), ylim=c(-.3,.3),
     cex=2, cex.lab= 1.5, xlab="PC 1", ylab='PC 3')
arrows(0, 0, pc1$rotation[,1],pc1$rotation[,3], len=.2, lwd=2, col="grey")
abline(h=0,v=0,col="black", lwd=1, lty=2)
pointLabel(1.05* pc1$rotation[,1] , 1.05* pc1$rotation[,3] , names(pc1$rotation[,1]),
           col="black", cex=1.5)








par(mfrow=c(1,1))

p1 <- strat.plot(lnSB2Chem, xRight =.7,
                yvar = depth, 
                y.rev= T, 
                #y.tks= depth,
                ylabel="Depth (feet)", 
                cex.ylabel = 1.5, 
                cex.yaxis = 1,
                srt.xlabel = 45,
                cex.xlabel= 1.5,
                scale.percent=F,
                plot.poly=TRUE, 
                col.poly= viridis(1), 
                col.poly.line=NA,
                xSpace=0.01, x.pc.lab=TRUE, x.pc.omit0=TRUE, las=2,
                clust=clust)

#lineCol <- adjustcolor(col="", alpha=.75)



p2 <- strat.plot(pcaZScores[,1:2], xLeft = 0.7, yvar = depth, y.rev=TRUE, 
           xRight=0.99,
           cex.ylabel = 1.5, 
           cex.yaxis = 1,
           srt.xlabel = 45,
           cex.xlabel= 1.5,
           plot.line=FALSE, 
           plot.poly=FALSE, 
           plot.bar=TRUE, 
           col.bar= viridis(14)[sb2Chem$Layer],
           lwd.bar=9,
           sep.bar=TRUE,
           las=2,
           y.axis=FALSE,  add=TRUE)




addClustZone(p1, clust,6, col="grey",
             lwd=2, lty=1)
addClustZone(p2, clust,6, col="grey",
             lwd=2, lty=1)







