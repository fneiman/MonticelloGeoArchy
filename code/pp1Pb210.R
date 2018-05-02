
library(mgcv)
library (dplyr)
library(ggplot2)

Pb210 <- read.csv('data/Pb210.csv', stringsAsFactors = F)

Pb210$depthFt <- 0.0328084* Pb210$Depth_cm

# check to see what it looks like
str(Pb210)

gy <- filter(Pb210, Site == 'Monticello Grave Yard' )



gy$depthFt2 <- gy$depthFt^2 

pp <- filter(Pb210, Site == 'PawPaw Alignment' )



fit <- lm(ExcessPb210 ~ depthFt + depthFt2,  data=gy)
summary(fit)

depthFt = seq(2.5 * .0328084, 17.5*.0328084,.01) 
depthFt2 <- depthFt^2

pred <- predict(fit, newdata=data.frame(depthFt, depthFt2) , interval='prediction')

par(mar=c(6,6,2,2))
plot(depthFt, pred[,1], ylim= c(-1,4), xlim=c(0,1.5), 
     type='l', col='blue', lwd=3,
     xlab='Depth (ft.)',
     ylab= 'Excess 210Pb',
     cex.lab= 2 , cex.axis=1.5)
lines(depthFt, pred[,2], col='blue')
lines(depthFt, pred[,3], col='blue')
points(pp$depthFt,pp$ExcessPb210, pch=21,bg=adjustcolor('gold',alpha=.75), cex=3 )

points(jitter(gy$depthFt),gy$ExcessPb210, pch=21,
       bg=ifelse(gy$Date == 1860, adjustcolor('blue',alpha=.5),adjustcolor('red',alpha=.5)),
       cex=3 )


legend('topright', legend= c('Alexander Garrett: 1860',
                             'Mary Jefferson Randolf: 1878',
                             'Paw Paw'),
        pch= c(21,21, 21),
        pt.bg = c(adjustcolor('blue',alpha=.5), 
                  adjustcolor('red',alpha=.5),
                  adjustcolor('gold',alpha=.5)), 
        pt.cex = 3,
        cex=1.5)



