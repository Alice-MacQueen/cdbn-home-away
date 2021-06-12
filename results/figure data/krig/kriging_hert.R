#Heritability Kriging
############
#Kriging of heritability
library(sp)
library(raster)
library(gstat)
library(lme4)
library(lmerTest)
library(dplyr) # for "glimpse"
library(ggplot2)
library(scales) # for "comma"
library(magrittr)
biocLite(c("GO.db", "preprocessCore", "impute"))
#Make Grid
library(ggmap)
library(maps)
library(mapdata)
library(spData) # contains datasets used in this example
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(plotdap)
library(geofi)
library(sf)
library(dplyr)
library(ggplot2)
library(rgdal)
library(oce)
library(raster)
library(gridExtra)
library(fields)

libs = c('ggplot2', 
         'ggpmisc',     # annotate
         'wesanderson', 
         'here',        # local paths in projects
         'magrittr')

for (i in libs) require(i, character.only=TRUE)


in_dir = file.path('results', 'figure data', 'krig')

color_palette = 'FantasticFox1'
colors = wes_palette(color_palette)
pal <- colorRampPalette(colors)


hert2<-here::here(in_dir, 'Bean_heritability_by_location_V2.csv') %>% 
  read.csv(stringsAsFactors=FALSE)
names(hert2)

hert2<-na.omit(hert2)
hert2 %>% as.data.frame %>% 
  ggplot(aes(Longitude, Latitude)) + geom_point(aes(size=Heritability), color="blue", alpha=3/4) + 
  ggtitle("Hertitability") + coord_equal() + theme_bw()


###get summary for kriging
hert<-here::here(in_dir, 'Herittability_year_loc_beans.csv') %>% 
  read.csv(stringsAsFactors=FALSE)
head(hert)

#### BEGIN PME EDITS ####
hert2 %<>% merge(unique(hert[, c('Site', 'Latitude', 'Longitude')]), by=c('Latitude', 'Longitude')) %>% 
  subset(Site != 'PRIS')

# hert2<-na.omit(hert)
# str(hert2)
#####Get summary by Treatment
#Sumerize the four yield by treatment, here we are exploring the mean
# 
# comb_stats = group_by(hert2, Site) %>%  # <- remember to group by the important factor
#   summarize(occur=sum(MATCH),Lat=mean(Latitude), Lon=mean(Longitude),Means = mean(Heritability), SD= sd(Heritability), SE = sd(Heritability)/sqrt(n()), 
#             CV = cv(Heritability))
# 
# comb_stats<-subset(comb_stats, Site != "PRIS")
# x<-cbind(comb_stats$Lon, comb_stats$Lat)
# y<-comb_stats$Means

x = hert2[, c('Longitude', 'Latitude')]
y = hert2[, 'Heritability']

#### END PME EDITS ####


#fit<- Tps(x,y)
#intras <- interpolate(y, spline)
#intras 
# fits a GCV thin plate smoothing spline surface to ozone measurements.
# Hey, it does not get any easier than this!

fit<-Krig(x,y)
surface(fit, type="C")
plot(fit)
summary(fit) #diagnostic summary of the fit 

fit2<- spatialProcess(x,y)
#fit2<-MLESpatialProcess.fast(x, y)
summary(fit2)
predict(fit2)

#save(comb_stats, file = "C:/Users/Michael/Downloads/comb_stats.rda")
#save(fit2, file = "C:/Users/Michael/Downloads/fit.rda")

#plot(fit2, xlab="Longitude", ylab="Latitude", main="1983")
surface(fit2, xlim=c(-130,-60), ylim=c(25,55),xlab="Longitude", 
        ylab="Latitude", main ="Heritability", breaks=c(0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6, 0.65,0.7,0.75,0.85),
        col=pal(11),add=T)

world(add=T)
US(add=TRUE)

# predictions on a grid 
surface( fit2, xlab="longitude", ylab="latitude")
US( add=TRUE, col="grey", lwd=2)
#(see also predictSurface for more control on evaluation grid
# and plotting)

# prediction standard errors, note two steps now to generate and then plot surface
look<- predictSurfaceSE(fit2)
surface( look, xlim=c(-130,-60), ylim=c(25,55), xlab="longitude", ylab="latitude")
world(add=T)
US(add=TRUE)
points( x, col="magenta")
title("prediction standard errors (Heritability)")

###make raster
rnglat <- range(hert2$Latitude)*c(0.95, 1.05)
rnglon <- range(hert2$Longitude)*c(1.05, 0.95)
xvals <- seq(rnglon[1], rnglon[2], len=120)  # PME I also increased resolution an order of magnitude. 
yvals <- seq(rnglat[1], rnglat[2], len=120)
griddf <- expand.grid(xvals, yvals)
griddf$pred <- predict(fit2, x=as.matrix(griddf) )

names(griddf)[1]<- "x"
names(griddf)[2]<- "y"
names(griddf)[3]<- "z"

griddf<-data.frame(griddf)

r = rasterFromXYZ(griddf)
crs(r) = crs('+init=epsg:4326')  # PME added CRS
plot(r)
writeRaster(r, here::here('results', 'figure data', 'kriged_bean_heritability.tif'), overwrite=TRUE)
#writeRaster(r, "C:/Users/Michael/Downloads/test_kriged_bean_map.tif")





