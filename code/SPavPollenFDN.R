
# setwd("P:/Sediment Analyses/Pollen/Data/SouthPavilion")

library(dplyr)
library(ggplot2)
library(ggrepel)
library(viridis)
library(ca)

pp1 <- read.csv('https://raw.github.com/fneiman/MonticelloGeoArchy/master/data/pawPawUnit1Pollen.csv', 
                stringsAsFactors = F)

# read the first three rows and transpose them into a data frame
labels <- read.csv('https://raw.github.com/fneiman/MonticelloGeoArchy/master/data/montSouthPavilion2018.csv', 
                   header=F, skip=1, sep=",", nrows = 4,
                   stringsAsFactors=F,colClasses='character' )  
labels1 <- data.frame(t(labels),  stringsAsFactors=F)
labels2 <- labels1[2:nrow(labels1),] 
names(labels2) <- as.character(labels1[1,])
labels2$Elevation <- as.numeric(labels2$Elevation)

str(labels2)

# read the pollen counts rows and transpose them into a data frame
pdata<- read.table('https://raw.github.com/fneiman/MonticelloGeoArchy/master/data/montSouthPavilion2018.csv', 
                   header=F, skip=6, sep=",")
                    ############
pdata1 <- data.frame(t(pdata),  stringsAsFactors=F)
pdata2 <- pdata1[2:nrow(pdata1),] 
names(pdata2) <- as.character(pdata1[1,])
pdata3<-data.frame(lapply(pdata2,as.numeric))
pdata3[is.na(pdata3)] <-0

# glue the two together
SPav<-data.frame(labels2,pdata3)
# get rid of empty records
SPav <- filter(SPav, !is.na(Elevation) )
# check sample sizes
table(SPav$Total.Pollen.Sum)
# get rid of tiny samples
SPav1 <- filter(SPav, Total.Pollen.Sum > 100 )
# put the taxon counts in a DF with rownames as context and sample number
SPav1Mat <- SPav1[,5:59] 
rownames(SPav1Mat) <- paste(SPav1[,2], SPav1[,1], sep = "-")



# do the ca
caResult <-ca(SPav1Mat)

# plot the ca
# put the result in dataframes
inertia <- data.frame('Inertia' = prop.table(caResult$sv^2))
rowScores <- data.frame(caResult$rowcoord, rownames=caResult$rownames)
colScores <- data.frame(caResult$colcoord, rownames=caResult$colnames)


# Compute the broken stick model inertia
broken.stick <- function(p)
  # Compute the expected values of the broken-stick distribution for 'p' pieces.
  # Example: broken.stick.out.20 = broken.stick(20)
  #             Pierre Legendre, April 2007
{
  result = matrix(0,p,2)
  colnames(result) = c("Dim","Expected.Inertia")
  for(j in 1:p) {
    E = 0
    for(x in j:p) E = E+(1/x)
    result[j,1] = j
    result[j,2] = E/p
  }
  result <- result
  return(data.frame(result))
}

bs <- broken.stick(nrow(inertia))



# plot the proportion of inertia
theme_set(theme_classic(base_size = 18))
p <- ggplot(data=inertia , aes(x= 1:length(Inertia), y=Inertia)) +
  geom_bar(stat="identity", fill="grey") +
  xlab ('Dimension') + 
  ylab( "Proportion of Inertia") +
  geom_line(aes(y = bs[,2], x= bs[,1]), color = "black", linetype = "dashed", size=1)
p


#######################################################
# scatter plots of dim scores

set.seed(42)
ggplot(data = rowScores, aes(x=Dim1, y=Dim2, fill=Dim3))  + 
  geom_point(shape=21,  alpha = .75, size= 6) + 
  scale_fill_viridis() +
  labs(fill='Dimension 3') +
  xlab (paste('Dimension 1',':   ', format(inertia[1,]*100, digits=2), '%')) + 
  ylab (paste('Dimension 2',':   ', format(inertia[2,]*100, digits=2), '%')) +
  coord_fixed(ratio=1) +
  geom_text_repel(aes(label = rownames(rowScores)),
                  size = 4, point.padding = .5)
    

ggplot(data = colScores, aes(x=Dim1, y=Dim2, fill=Dim3))  + 
  geom_point(shape=21,  alpha = .75, size= 6) + 
  xlab (paste('Dimension 1',':   ', format(inertia[1,]*100, digits=2), '%')) + 
  ylab (paste('Dimension 2',':   ', format(inertia[2,]*100, digits=2), '%')) +
  scale_fill_viridis() +
  labs(fill='Dimension 3') +
  coord_fixed(ratio=1) +
  geom_text_repel(aes(label = rownames(colScores)),
                  size = 4, point.padding = .5)

# do it again without the three  outliers: samples 11 17 19
# define the function to  subset of the rows and cols
dumpOutliers <- function(dataframe, outliers) {
  df1 <- subset(dataframe, ! rownames(dataframe) %in% outliers)
  df2 <- df1[, colSums(df1) > 0]
  return (df2)

}

SPav2Mat <- dumpOutliers(SPav1Mat,c('2588F-11', '2588J-19', '2588F-17')) 


# do the ca
caResult <-ca(SPav2Mat)

# plot the ca
# put the result in dataframes
inertia <- data.frame('Inertia' = prop.table(caResult$sv^2))
rowScores <- data.frame(caResult$rowcoord, rownames=caResult$rownames)
colScores <- data.frame(caResult$colcoord, rownames=caResult$colnames)


# Compute the broken stick model inertia
broken.stick <- function(p)
  # Compute the expected values of the broken-stick distribution for 'p' pieces.
  # Example: broken.stick.out.20 = broken.stick(20)
  #             Pierre Legendre, April 2007
{
  result = matrix(0,p,2)
  colnames(result) = c("Dim","Expected.Inertia")
  for(j in 1:p) {
    E = 0
    for(x in j:p) E = E+(1/x)
    result[j,1] = j
    result[j,2] = E/p
  }
  result <- result
  return(data.frame(result))
}

bs <- broken.stick(nrow(inertia))



# plot the proportion of inertia
theme_set(theme_classic(base_size = 18))
p <- ggplot(data=inertia , aes(x= 1:length(Inertia), y=Inertia)) +
  geom_bar(stat="identity", fill="grey") +
  xlab ('Dimension') + 
  ylab( "Proportion of Inertia") +
  geom_line(aes(y = bs[,2], x= bs[,1]), color = "black", linetype = "dashed", size=1)
p


#######################################################
# scatter plots of dim scores

# get the indices of the taxa
rownames(colScores)

Group <- c(rep('Herbs and Cultigens',27 ), rep('Arboreal', 50-27), 'Indeterminant')
colScores$Group <- Group  
##### Working on this!
set.seed(42)
ggplot(data = rowScores, aes(x=Dim1, y=Dim2))  + 
  geom_point(shape=21,  fill= viridis(1), alpha = .75, size= 6) + 
  xlab (paste('Dimension 1',':   ', format(inertia[1,]*100, digits=2), '%')) + 
  ylab (paste('Dimension 2',':   ', format(inertia[2,]*100, digits=2), '%')) +
  coord_fixed(ratio=1) +
  geom_text_repel(aes(label = rownames(rowScores)),
                  size = 4, point.padding = .5)


ggplot(data = colScores, aes(x=Dim1, y=Dim2, fill=Group))  + 
  geom_point(shape=21,  alpha = .75, size= 6) + 
  xlab (paste('Dimension 1',':   ', format(inertia[1,]*100, digits=2), '%')) + 
  ylab (paste('Dimension 2',':   ', format(inertia[2,]*100, digits=2), '%')) +
  scale_fill_viridis(discrete = T) +
  labs(fill='Group') +
  coord_fixed(ratio=1) +
  geom_text_repel(aes(label = rownames(colScores)),
                  size = 4, point.padding = .5)






