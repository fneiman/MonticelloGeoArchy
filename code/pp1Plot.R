
library(rioja)
library(dplyr)

pp1 <- read.csv('data/pawpaw1Data.csv', stringsAsFactors = F)


# check to see what it looks like
str(pp1)
# what vars are in what rows?
cbind(1:ncol(pp1), colnames(pp1))

pp1 <- mutate(pp1,Fabaceae = Fabaceae.1 )


# leave out "unknown x", other, indeterminant, and sums
taxonCounts <- pp1[,c(9:117,123)]
# get rid of taxa with all NAs
taxonCounts <- taxonCounts[,!is.na(colSums(taxonCounts))]
# get rid of taxa with all 0s 
taxonCounts <- taxonCounts[, colSums(taxonCounts) > 0] 

# histgram of number of occurrences
hist(colSums(taxonCounts > 0), breaks=1:50, col='grey') 

# get rid of taxa with single occurrences
taxonCounts <- taxonCounts[,colSums(taxonCounts > 0) > 5]

# get rig of taxa with a max < 3 in a given sample



# create appropriately sized graphics window
windows(width=16, height=10)

taxonPcts <- prop.table(as.matrix(taxonCounts), 1)*100
rowSums(taxonPcts)

elevation <- pp1$Elevation
depth <- pp1$Depth 


# add a dendrogram from constrained cluster analysis
diss <- dist(taxonPcts, method="manhattan" )
clust <- chclust(diss, method="coniss")
# broken stick model suggest 3 significant zones
bstick(clust)

title= 'Paw Paw: Unit1'
xx <- 1:ncol(taxonPcts)
tColors <- ifelse(xx < 10, "darkgreen", "darkred")
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

addClustZone(x, clust, 4, col="black")









