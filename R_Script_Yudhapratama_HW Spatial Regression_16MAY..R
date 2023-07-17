################################################################################
#  Analisis Kuantitatif untuk Kebijakan Publik                           #
# Topik: Spatial Regression Analysis                                          #
# Yudhapratama Nugraha.                                              #
################################################################################
##Data Prep
#Panggil library
library(spdep)
library(maptools)
library(rgdal)
library(ggplot2)
library(car)
library(RColorBrewer)
library(rgeos)
library(rgeoda)
library(sf)
library(stargazer)
library(spatialreg)

#Set working directory
getwd()
setwd('/Volumes/Samsung T7 SSD/video Kuliah/Met Kuan Video/Precious File Pak Ade/PR_Spatial')
dir()

#Menggunakan data Spasial
nat_east <- readOGR(".","NAT_EAST")
head(nat_east)
summary(nat_east)
nat_east_data <-data.frame(nat_east)
head(nat_east_data)
nat_global <- nat_east
plot(nat_global)

#Subset biar lebih mengerucut
subset(nat_east, EAST=="1")
nat_east_only<-subset(nat_east, EAST=="1")
nat_east_only_data <- data.frame(nat_east_only)
plot(nat_east_only, main = "East Coast", axes = TRUE)
summary(nat_east_only)
summary(nat_east_only_data)
head(nat_east_only_data)
coords_east <- coordinates(nat_east_only)
head(coords_east)

###Menentukan Bobot Spasial dengan k-nearest neighborhood
#memasukan knn
ID_NAT_East <- row.names(as(nat_east_only, "data.frame"))
nat_east_nnb5 <- knn2nb(knearneigh(coords_east, k=5), row.names=ID_NAT_East)
nat_east_nnb10 <- knn2nb(knearneigh(coords_east, k=10), row.names=ID_NAT_East)
#Plot
par(mfrow=c(1,1))
plot(nat_east_only)
plot(nat_east_nnb5, coords_east, col = "blue", add = T)

par(mfrow=c(1,1))
plot(nat_east_only)
plot(nat_east_nnb10, coords_east, col = "green", add = T)
#memindahkan pada bendtuk list
nat_east_nnb5_wl <- nb2listw(nat_east_nnb5, style = "W")
nat_east_nnb10_wl <- nb2listw(nat_east_nnb10, style = "W")
summary(nat_east_nnb5_wl)

#mensimetriskan
nat_east_symnnb5 <- make.sym.nb(nat_east_nnb5)
nat_east_symnnb5_wl <- nb2listw(nat_east_symnnb5, style = "W")
nat_east_symnnb10 <- make.sym.nb(nat_east_nnb10)
nat_east_symnnb10_wl <- nb2listw(nat_east_symnnb10, style = "W")

############################################################
#1st Statistics, Global Moran's I
ID_NAT <- row.names(as(nat_global, "data.frame"))
coords <- coordinates(nat_global)
nat_global_10nnb <-knn2nb(knearneigh(coords, k=10), row.names=ID_NAT)
nat_global_10nnb_wl <- nb2listw(nat_global_10nnb, style = "W")

#############################################################

#1st Step, Global Moran's I East-Coast
moran.test(nat_east_only$HR60,nat_east_nnb10_wl, randomisation = FALSE, alternative="two.sided" )
moran.test(nat_east_only$HR70,nat_east_nnb10_wl, randomisation = FALSE, alternative="two.sided" )
moran.test(nat_east_only$HR80,nat_east_nnb10_wl, randomisation = FALSE, alternative="two.sided" )
moran.test(nat_east_only$HR90,nat_east_nnb10_wl, randomisation = FALSE, alternative="two.sided" )
#Plot Moran's I Test
par(mfrow=c(1,4)) 
moran.plot(nat_east_only$HR60, nat_east_nnb10_wl)
moran.plot(nat_east_only$HR70, nat_east_nnb10_wl)
moran.plot(nat_east_only$HR80, nat_east_nnb10_wl)
moran.plot(nat_east_only$HR90, nat_east_nnb10_wl)
##Dari Hasil tersebut dapat disimpulkan bahwa terjadi spetial Correlations.

#2nd Step, local Moran's I East-Coast 60s
lisa_hr60 <- localmoran_perm(nat_east_only$HR60, nat_east_nnb10_wl)
head(lisa_hr60)
DV.60 <- nat_east_only$HR60
quadrant <- vector(mode="numeric",length=nrow(lisa_hr60))
cDV.60 <- DV.60 - mean(DV.60) 
lagDV.60 <- lag.listw(nat_east_nnb10_wl, DV.60)
clagDV.60 <- lagDV.60 - mean(lagDV.60)

#Preparasi Quadran
quadrant <- 5
quadrant[cDV.60 >0 & clagDV.60>0 & lisa_hr60[,6]<=.05] <- 1 
quadrant[cDV.60 <0 & clagDV.60<0 & lisa_hr60[,6]<=.05] <- 2      
quadrant[cDV.60 <0 & clagDV.60>0 & lisa_hr60[,6]<=.05] <- 3
quadrant[cDV.60 >0 & clagDV.60<0 & lisa_hr60[,6]<=.05] <- 4

#Plot Lisa HR60
brks <- c(1,2,3,4,5)
colors <- c("#FF0000", "#0000FF", "#a7adf9", "#f4ada8", "#eeeeee")

par(mar=c(0,0,1,0))
plot(nat_east_only, border = "#333333", lwd=0.2, col=colors[findInterval(quadrant,brks,all.inside=FALSE)], axes = TRUE)
# box()
legend('bottomleft',legend=c("High-High","Low-Low","Low-High","High-Low", "Not Significant"),
       fill=colors, border = "#eeeeee")
title("LISA East-Coast Cluster Map, 1960 Homicides")


#Plot Significance
class.p <- 0

class.p[lisa_hr60[,6]>.01 & lisa_hr60[,6]<=.05] <- 1 
class.p[lisa_hr60[,6]>0.001 & lisa_hr60[,6]<=.01] <- 2 
class.p[lisa_hr60[,6]<=.001] <- 2 


brks.p <- c(0,1,2,3)
colors.p <- c("#eeeeee", "#84f576", "#53c53c", "#348124")

par(mar=c(0,0,1,0)) # sets margin parameters for plot space
plot(nat_east_only, border = "#333333", lwd=0.2, col=colors.p[findInterval(class.p,brks.p,all.inside=FALSE)], axes = TRUE)
# box()
legend('bottomleft',legend=c("Not significant", "p <= 0.05", "p <= 0.01", "p <= 0.001"),
       fill=colors.p , border = "#eeeeee")
title("East Coast Local Moran Map of HR60 (P-Values)")




#Local Moran's I East-Coast 70s
lisa_hr70 <- localmoran_perm(nat_east_only$HR70, nat_east_nnb10_wl)
head(lisa_hr70)
DV.70 <- nat_east_only$HR70
quadrant <- vector(mode="numeric",length=nrow(lisa_hr70))
cDV.70 <- DV.70 - mean(DV.70) 
lagDV.70 <- lag.listw(nat_east_nnb10_wl, DV.70)
clagDV.70 <- lagDV.70 - mean(lagDV.70)

#Preparasi Quadran
quadrant <- 5
quadrant[cDV.70 >0 & clagDV.70>0 & lisa_hr70[,6]<=.05] <- 1 
quadrant[cDV.70 <0 & clagDV.70<0 & lisa_hr70[,6]<=.05] <- 2      
quadrant[cDV.70 <0 & clagDV.70>0 & lisa_hr70[,6]<=.05] <- 3
quadrant[cDV.70 >0 & clagDV.70<0 & lisa_hr70[,6]<=.05] <- 4

#Plot Lisa HR70
brks <- c(1,2,3,4,5)
colors <- c("#FF0000", "#0000FF", "#a7adf9", "#f4ada8", "#eeeeee")

par(mar=c(0,0,1,0))
plot(nat_east_only, border = "#333333", lwd=0.2, col=colors[findInterval(quadrant,brks,all.inside=FALSE)], axes = TRUE)
# box()
legend('bottomleft',legend=c("High-High","Low-Low","Low-High","High-Low", "Not Significant"),
       fill=colors, border = "#eeeeee")
title("LISA East-Coast Cluster Map, 1970 Homicides")


#Plot Significance
class.p <- 0

class.p[lisa_hr70[,6]>.01 & lisa_hr70[,6]<=.05] <- 1 
class.p[lisa_hr70[,6]>0.001 & lisa_hr70[,6]<=.01] <- 2 
class.p[lisa_hr70[,6]<=.001] <- 2 


brks.p <- c(0,1,2,3)
colors.p <- c("#eeeeee", "#84f576", "#53c53c", "#348124")

par(mar=c(0,0,1,0)) # sets margin parameters for plot space
plot(nat_east_only, border = "#333333", lwd=0.2, col=colors.p[findInterval(class.p,brks.p,all.inside=FALSE)], axes = TRUE)
# box()
legend('bottomleft',legend=c("Not significant", "p <= 0.05", "p <= 0.01", "p <= 0.001"),
       fill=colors.p , border = "#eeeeee")
title("East Coast Local Moran Map of HR70 (P-Values)")


#Local Moran's I East-Coast 80 & 90s
lisa_hr80 <- localmoran_perm(nat_east_only$HR80, nat_east_nnb10_wl)
lisa_hr90 <- localmoran_perm(nat_east_only$HR90, nat_east_nnb10_wl)
DV.80 <- nat_east_only$HR80
DV.90 <- nat_east_only$HR90
quadrant80 <- vector(mode="numeric",length=nrow(lisa_hr80))
quadrant90 <- vector(mode="numeric",length=nrow(lisa_hr90))
cDV.80 <- DV.80 - mean(DV.80) 
cDV.90 <- DV.90 - mean(DV.90)
lagDV.80 <- lag.listw(nat_east_nnb10_wl, DV.80)
lagDV.90 <- lag.listw(nat_east_nnb10_wl, DV.90)
clagDV.80 <- lagDV.80 - mean(lagDV.80)
clagDV.90 <- lagDV.90 - mean(lagDV.90)

#Preparasi Quadran
quadrant80 <- 5
quadrant[cDV.80 >0 & clagDV.80>0 & lisa_hr80[,6]<=.05] <- 1 
quadrant[cDV.80 <0 & clagDV.80<0 & lisa_hr80[,6]<=.05] <- 2      
quadrant[cDV.80 <0 & clagDV.80>0 & lisa_hr80[,6]<=.05] <- 3
quadrant[cDV.80 >0 & clagDV.80<0 & lisa_hr80[,6]<=.05] <- 4

quadrant90 <- 6
quadrant[cDV.90 >0 & clagDV.90>0 & lisa_hr90[,6]<=.05] <- 1 
quadrant[cDV.90 <0 & clagDV.90<0 & lisa_hr90[,6]<=.05] <- 2      
quadrant[cDV.90 <0 & clagDV.90>0 & lisa_hr90[,6]<=.05] <- 3
quadrant[cDV.90 >0 & clagDV.90<0 & lisa_hr90[,6]<=.05] <- 4

#Plot Lisa HR80&90
brks80 <- c(1,2,3,4,5)
colors <- c("#FF0000", "#0000FF", "#a7adf9", "#f4ada8", "#eeeeee")

brks90 <- c(1,2,3,4,6)
colors <- c("#FF0000", "#0000FF", "#a7adf9", "#f4ada8", "#eeeeee")

par(mar=c(0,0,1,0))
plot(nat_east_only, border = "#333333", lwd=0.2, col=colors[findInterval(quadrant,brks80,all.inside=FALSE)], axes = TRUE)
# box()
legend('bottomleft',legend=c("High-High","Low-Low","Low-High","High-Low", "Not Significant"),
       fill=colors, border = "#eeeeee")
title("LISA East-Coast Cluster Map, 1980 Homicides")

par(mar=c(0,0,1,0))
plot(nat_east_only, border = "#333333", lwd=0.2, col=colors[findInterval(quadrant,brks90,all.inside=FALSE)], axes = TRUE)
# box()
legend('bottomleft',legend=c("High-High","Low-Low","Low-High","High-Low", "Not Significant"),
       fill=colors, border = "#eeeeee")
title("LISA East-Coast Cluster Map, 1990 Homicides")


#Plot Significance 1980
class.p <- 0

class.p[lisa_hr80[,6]>.01 & lisa_hr80[,6]<=.05] <- 1 
class.p[lisa_hr80[,6]>0.001 & lisa_hr80[,6]<=.01] <- 2 
class.p[lisa_hr80[,6]<=.001] <- 2 


brks80.p <- c(0,1,2,3)
colors.p <- c("#eeeeee", "#84f576", "#53c53c", "#348124")

par(mar=c(0,0,1,0)) # sets margin parameters for plot space
plot(nat_east_only, border = "#333333", lwd=0.2, col=colors.p[findInterval(class.p,brks.p,all.inside=FALSE)], axes = TRUE)
# box()
legend('bottomleft',legend=c("Not significant", "p <= 0.05", "p <= 0.01", "p <= 0.001"),
       fill=colors.p , border = "#eeeeee")
title("East Coast Local Moran Map of HR80 (P-Values)")


#Plot Significance 90
class90.p <- 0

class90.p[lisa_hr90[,6]>.01 & lisa_hr90[,6]<=.05] <- 1 
class90.p[lisa_hr90[,6]>0.001 & lisa_hr90[,6]<=.01] <- 2 
class90.p[lisa_hr90[,6]<=.001] <- 2 


brks90.p <- c(0,1,2,3)
colors90.p <- c("#eeeeee", "#84f576", "#53c53c", "#348124")

par(mar=c(0,0,1,0)) # sets margin parameters for plot space
plot(nat_east_only, border = "#333333", lwd=0.2, col=colors.p[findInterval(class.p,brks90.p,all.inside=FALSE)], axes = TRUE)
# box()
legend('bottomleft',legend=c("Not significant", "p <= 0.05", "p <= 0.01", "p <= 0.001"),
       fill=colors.p , border = "#eeeeee")
title("East Coast Local Moran Map of HR90 (P-Values)")





#####################
###RegresiOLS
#####################
#Persiapan Data
nat_east_dataR<-nat_east_only@data
nat_east_dataR
data_nat60 <-nat_east_dataR[,c("HR60","RD60","PS60","MA60","DV60","UE60")]
head(data_nat60)
data_nat70 <-nat_east_dataR[,c("HR70","RD70","PS70","MA70","DV70","UE70")]
head(data_nat70)
data_nat80 <-nat_east_dataR[,c("HR80","RD80","PS80","MA80","DV80","UE80")]
head(data_nat80)
data_nat90 <-nat_east_dataR[,c("HR90","RD90","PS90","MA90","DV90","UE90")]
head(data_nat90) 

#Regresi linear Tahun 60
nat_east60_lm <- lm(HR60~ RD60 + PS60 + MA60 + DV60 + UE60, data = data_nat60)
stargazer(nat_east60_lm, title = "Hasil OLS 1960", align = TRUE, type = "text", dep.var.labels = c("Tingkat Pembunuhan"), covariate.labels = c("Resource Deprivation", "Population Structure", "Median Age", "Percent of males divorced", "Percent of Male Unemployed"))

#Regresi linear Tahun 70
nat_east70_lm <- lm(HR70~ RD70 + PS70 + MA70 + DV70 + UE70, data = data_nat70)
stargazer(nat_east70_lm, title = "Hasil OLS 1970", align = TRUE, type = "text", dep.var.labels = c("Tingkat Pembunuhan"), covariate.labels = c("Resource Deprivation", "Population Structure", "Median Age", "Percent of males divorced", "Percent of Male Unemployed"))

#Regresi Linear Tahun 80
nat_east80_lm <- lm(HR80~ RD80 + PS80 + MA80 + DV80 + UE80, data = data_nat80)
stargazer(nat_east80_lm, title = "Hasil OLS 1980", align = TRUE, type = "text", dep.var.labels = c("Tingkat Pembunuhan"), covariate.labels = c("Resource Deprivation", "Population Structure", "Median Age", "Percent of males divorced", "Percent of Male Unemployed"))

#Regresi Linear Tahun 90
nat_east90_lm <- lm(HR90~ RD90 + PS90 + MA90 + DV90 + UE90, data = data_nat90)
stargazer(nat_east60_lm, title = "Hasil OLS 1990", align = TRUE, type = "text", dep.var.labels = c("Tingkat Pembunuhan"), covariate.labels = c("Resource Deprivation", "Population Structure", "Median Age", "Percent of males divorced", "Percent of Male Unemployed"))

#Uji Diagnostik, LM Test
lm.LMtests(nat_east60_lm, nat_east_symnnb10_wl, test = "all")
lm.LMtests(nat_east70_lm, nat_east_symnnb10_wl, test = "all")
lm.LMtests(nat_east80_lm, nat_east_symnnb10_wl, test = "all")
lm.LMtests(nat_east90_lm, nat_east_symnnb10_wl, test = "all")

#Spatial Durbin Model Test 1980
install.packages("sandwich")
library(spatialreg)
library(lmtest)
library(sandwich)
SDM_East <- lagsarlm(HR90~ RD90 + PS90 + MA90 + DV90 + UE90, data = data_nat90, listw = nat_east_symnnb10_wl, type = "mixed", method="Matrix", trs=T,control=list(), tol.solve = 1e-12 )
summary(SDM_East)

##Menentukan Direct & non-Direct Effect
#Memperbaiki Std Error
KoefSDM <- as.data.frame(SDM_East$coefficients)
lm.target <- lm(SDM_East$tary ~ SDM_East$tarX - 1)
robust.est <- coeftest(lm.target, vcov. = vcovHC(lm.target, type="HC0", df=Inf))
rse <- as.vector(robust.est[,2])
summary(rsdm)

#Memasukan pada Matriks
w.t <- as(as_dgRMatrix_listw(nat_east_symnnb10_wl), "CsparseMatrix")
trMC <- trW(w.t, type = "MC")
SDM_Impact <- summary(impacts(rsdm, tr=trMC, R = 1000), zstats = T)
SDM_Impact
#!!---- ada bukti dari Spillover

x.list <- names(SDM_Impact$res$direct)
x.list
#Partisi Impact berdasarkan Orde
SDM_Impact_par <- summary(impacts(rsdm, tr=trMC, R=1000 , Q=4), zstats= T, short = T, Q = 4 , reportQ = T)
SDM_Impact_par

#Visualisasikan Impact Berdasarkan Orde
n.w <- 4
x_direct <- c(0,1,2,3)
x_indirect <- c(0,1,2,3)

for(k in 1:length(x.list)){
  mean.val.d <- as.data.frame(SDM_Impact_par$Qdirect_sum$statistics[((k*n.w)-(n.w-1)):(k*n.w),1])
  colnames(mean.val.d) <-paste0("est_", x.list[k])
  sd.val.d <- as.data.frame(SDM_Impact_par$Qdirect_sum$statistics[((k*n.w)-(n.w-1)):(k*n.w),2])
  colnames(sd.val.d)<-paste0("sd_",x.list[k])
  
  x_direct <- cbind(x_direct, mean.val.d, sd.val.d)
  rownames(x_direct)<-seq(n.w)
                                                      
  mean.val.ind <-as.data.frame(SDM_Impact_par$Qindirect_sum$statistics[((k*n.w)-(n.w-1)):(k*n.w),1])
  colnames(mean.val.ind)<-paste0("est_", x.list[k])
  sd.val.ind <- as.data.frame(SDM_Impact_par$Qindirect_sum$statistics[((k*n.w)-(n.w-1)):(k*n.w),2])
  colnames(sd.val.ind) <- paste0("sd_",x.list[k])
                                                      
   x_indirect <- cbind(x_indirect,mean.val.ind, sd.val.ind)
   rownames(x_indirect) <- seq(1:n.w)}

colnames(x_direct)[1] <- "W.Orde"
colnames(x_indirect)[1] <- "W.Orde"

#Panggil Hasil
x_direct
x_indirect

#Pemberian Label
est_var_list <- grep('est_', colnames(x_indirect), value = TRUE)
sd_var_list <- grep('sd_', colnames(x_indirect), value = TRUE)

est_var_label <- as.data.frame(grep('est_', colnames(x_indirect), value = TRUE))
colnames(est_var_label) <- "Variabel"  
est_var_label$var_label <- c('Resource Deprivation (RD)', 'Population Structure (PS)', 'Median Age (MA)', 'Percent of Males Divorced(DV)', 'Percent of male Unemployed(UE)')


###Visualisasi dari direct dan indirect Effect
par(mfrow=c(2,5))
for(k in 1:length(x.list)){
  plot <- ggplot(data = x_direct, aes = W.Ord, y = x_direct[,est_var_list[k])) +
  geom_hline(yintercept =0, colour= "black")+
  geom ribbon(aes(ymin = x_direct[est_var_list[k]] - 2.576 * x_direct[,sd_var_list[k]],
                  ymax = x_direct[est_var_list[k]] + 2.576 * x_direct[,sd_var_list[k]],
                  alpha = 0,4, colour="grey", fill = 'lightgrey', linetype = 'dashed', size 0,6)
  plota <-plot +
    )}

