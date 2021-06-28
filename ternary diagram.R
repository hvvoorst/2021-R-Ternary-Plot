'

ternary plot with calcualte height map based on Least-Quares fitting.
10 different blends of cement were compared in performance, specifically the unconfined compressive strength after firing to 815°C
Three types of cement were used, suitable for ternamry plotting.
The compositions were:
 - 3 pure cements
 - 6 2 cements in wt ratio 1:2 | 2:1
 - 1 3 cement in wt ratio 1:1:1


sources
https://cran.r-project.org/web/packages/Ternary/vignettes/Ternary.html

'


#-------------------- LOAD DATA --------------------
library(Ternary)

#all data (measured)
kdv815 <- c(34.7,43.7,49.6,50.4,44.6,47.7,50.9,54.1,45.5,47.7) #unconfined ocmpressive strength
plc815 <- c(-0.138,-0.154,-0.121,-0.133,-0.142,-0.117,-0.134,-0.143,-0.103,-0.149) #permanent lineair change
ds815 <- c(2145,2089,2050,2109,2145,2144,2112,2090,2133,2099) #apperent density
water <- c(350,375,400,330,325,330,340,340,320,350) #ml of water added (800g of cement, 2kg aggregate of size 0-4mm)

datapoints <- kdv815

#ternamry diagram positions resembling the weight fraction in cement blends
x <- c(1,0,0,1/3,2/3,1/3,0,0,1/3,2/3)
y <- c(0,1,0,1/3,1/3,2/3,2/3,1/3,0,0)
z <- c(0,0,1,1/3,0,0,1/3,2/3,2/3,1/3)

coordinates <- as.data.frame(cbind(x,y,z))

#-------------------- TIDY DATA --------------------
#-------------------- LOAD MODEL --------------------
library(dplyr)
#define the constants 
kvalues <- coordinates %>%
  mutate(xy=x*y,
         xz=x*z,
         yz=y*z,
         xyz=x*y*z)

#calculate solution wityh LSM
A <- as.matrix(kvalues)
b <- solve(t(A) %*% A) %*% t(A) %*% datapoints

#draw the solution
FunctionToContour <- function (x, y, z) {
  b[1]*x + b[2]*y + b[3]*z + b[4]*x*y + b[5]*x*z + b[6]*y*z + b[7]*x*y*z
}

calcvalues <- TernaryPointValues(FunctionToContour, resolution = 24L)

#-------------------- PLOT --------------------
#frame
TernaryPlot(point = 'up', atip = 'cem 40', btip = 'cem 50', ctip = 'cem 70',alab = '% cem 40', blab = '% cem 50', clab = '% cem 70')
#color map
ColourTernary(calcvalues)
#height contour lines
TernaryContour(FunctionToContour, resolution = 36L)
#yellow dots to write actual measurements on
TernaryPoints(coordinates[,c('x','y','z')],pch=21,cex=3,col = alpha("yellow", 0.5),bg=alpha("yellow",0.5))
#actual measurements 
TernaryText(coordinates[,c('x','y','z')],labels=as.character(datapoints))

legend('topleft', pch = NULL, 
       legend = "",
       title = 'KDV815°C cement blend', bty = 'n', cex = 1.2)

#-------------------- TERMINATE CODE --------------------

rm(list = ls())
cat("\014") #or cntrl+L
#if(!is.null(dev.list())) dev.off()
gc()
