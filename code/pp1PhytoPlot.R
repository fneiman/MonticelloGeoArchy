library(rioja)
library(dplyr)
library(stringr)
library(mgcv)

ppPhyt <- read.Tilia("P:/Sediment Analyses/Phytoliths/JonesPhytoData/Pawphyt.til")

baseTaxa <- ppPhyt$vars$Sums %in% c( 'A', 'B', 'C' , 'D' , 'E' , 'F', 'G')

ppPhytData <- ppPhyt$data[,baseTaxa]

colnames(ppPhytData) <- ppPhyt$vars$FullNames[baseTaxa]

# get rid of taxa with all 0s 
ppPhytData <- ppPhytData[, colSums(ppPhytData) > 0] 

# histogram of number of occurrences
hist(colSums(ppPhytData > 0), breaks=1:50, col='grey') 

# get rid of taxa with single occurrences
ppPhytData <- ppPhytData[,colSums(ppPhytData > 0) > 5]

str(ppPhytData)

 
ppPhytData<- rename(ppPhytData, Arboreal.Coarse = Arboreal ) 


names(ppPhytData)<-str_replace_all(names(ppPhytData), c(" " = ".") )
names(ppPhytData)<-str_replace_all(names(ppPhytData), c("/" = '_' ) )
names(ppPhytData)<-str_replace_all(names(ppPhytData), c("-" = '' ) )
names(ppPhytData)<-str_replace_all(names(ppPhytData), c('Fest.' = "Pooid" ) )

str(ppPhytData)

ppPhytData<- rename(ppPhytData, Pooid.Irregular =  Pooidcoid.Irregular  ) 

ppPhytData<- rename(ppPhytData, Chloridoid.Saddle = Saddle ) 
ppPhytData<- rename(ppPhytData, Panicoid.Cross = Cross ) 
ppPhytData<- rename(ppPhytData, Panicoid.Zea = Zea.mays ) 

ppPhytData<- rename(ppPhytData, Poaceae.Bilobate = Bilobate)              
ppPhytData<- rename(ppPhytData, Poaceae.Stipa.Type = Stipatype) 

str(ppPhytData)                                  

cbind(colnames(ppPhytData))

ppPhytData <- ppPhytData %>% select(
              Poaceae.Bulliform,     
              Poaceae.Elongate,      
              Poaceae.Hair_Edge.Cell,
              Poaceae.Bilobate,      
              Poaceae.Stipa.Type,    
              Pooid.Irregular,       
              Pooid.Keeled,          
              Pooid.Conical,         
              Pooid.Pyramidal,       
              Pooid.Crenate,         
              Chloridoid.Saddle,     
              Panicoid.Bilobate,     
              Panicoid.Cross,        
              Panicoid.Zea,          
              Bambusoidea,
              Hollow.spine,          
              Sabaltype.Palm,
              Bulliform,             
              Elongate.Rod,          
              Echinate.Irregular,    
              Sclereid,              
              Arboreal.Coarse,   
              Arboreal.Fine)         
  
  
  

Elevation <- as.numeric(as.character(ppPhyt$levels$Names))
Depth <- max(Elevation)-Elevation


ppPhytPcts <- prop.table(as.matrix(ppPhytData),1) *100


# add a dendrogram from constrained cluster analysis
diss <- dist(ppPhytPcts, method="euclidean" )
clust <- chclust(diss, method="coniss")
# broken stick model suggest 3 significant zones
bstick(clust)


title= 'Paw Paw: Unit 1'
xx <- 1:ncol(ppPhytPcts)
tColors <- ifelse(xx < 22, "darkgreen", "darkred")


x <- strat.plot(ppPhytPcts, 
                yvar = Depth, 
                y.rev=TRUE, 
                #y.tks= depth,
                ylabel="Depth (feet)", 
                cex.ylabel = 1.5, 
                cex.yaxis = 1,
                srt.xlabel = 45,
                cex.xlabel= 1.25,
                
                scale.percent=TRUE, 
                
                plot.poly=TRUE, 
                col.poly= tColors, 
                col.poly.line=NA,
                exag=TRUE,
                exag.mult =10,
                col.exag="auto",
                clust=clust)

addClustZone(x, clust, 4, col="black")


Poaceae <- rowSums(ppPhytData[,1:15])
Arboreal <-  rowSums(ppPhytData[,22:23])
Pooid <-  rowSums(ppPhytData[,6:10])
nonPooidPoaceae <- Poaceae - Pooid 


poaIndex <- data.frame(Poaceae, Arboreal, Depth = -Depth)
poaModel <- gam( cbind(Poaceae, Arboreal) ~ s(Depth),
                  link = logit,
                  family = quasibinomial,
                  method='REML',
                  data= poaIndex)

summary(poaModel)
anova(poaModel)
plot(poaModel)


acf(residuals(poaModel),ci.type = "ma")


invLogit <- function(x){ exp(x)/(1+exp(x)) }

pred <- predict(poaModel, type='link', se.fit=T) 
pHat <- invLogit(pred$fit)
upperCL <- invLogit(pred$fit + 2*pred$se.fit)
lowerCL <- invLogit(pred$fit - 2*pred$se.fit)

par(mar=c(6,6,2,2))
with(poaIndex, plot(Depth, Poaceae/(Poaceae + Arboreal),
                     pch=21, bg='green', cex=2,
                     cex.lab= 2,
                     cex.axis=1.5,
                    xlim= c(-5,0)
))

lines(poaIndex$Depth, pHat, lwd=3, col= 'darkgreen')
lines(poaIndex$Depth, upperCL, lwd=1, col= 'darkgreen', lty=2)
lines(poaIndex$Depth, lowerCL, lwd=1, col= 'darkgreen', lty=2)
abline( v=c(-.5, -2.5, -3.1), lty=1, lwd=2, col='grey')


##################################################
pooidIndex <- data.frame(Pooid, nonPooidPoaceae, Depth = -Depth)
pooidModel <- gam( cbind(Pooid, nonPooidPoaceae) ~ s(Depth),
                 link = logit,
                 family = quasibinomial,
                 #method='REML',
                 data= pooidIndex)

summary(pooidModel)
anova(pooidModel)
plot(pooidModel)


acf(residuals(poaModel),ci.type = "ma")
 

pred <- predict(pooidModel, type='link', se.fit=T) 
pHat <- invLogit(pred$fit)
upperCL <- invLogit(pred$fit + 2*pred$se.fit)
lowerCL <- invLogit(pred$fit - 2*pred$se.fit)

par(mar=c(6,6,2,2))
with(pooidIndex, plot(Depth, Pooid/(Pooid + nonPooidPoaceae),
                    pch=21, bg='green', cex=2,
                    cex.lab= 2,
                    cex.axis=1.5,
                    xlim= c(-5,0)
))

lines(pooidIndex$Depth, pHat, lwd=3, col= 'darkgreen')
lines(pooidIndex$Depth, upperCL, lwd=1, col= 'darkgreen', lty=2)
lines(pooidIndex$Depth, lowerCL, lwd=1, col= 'darkgreen', lty=2)
abline( v=c(-.5, -2.5, -3.1), lty=1, lwd=2, col='grey')



