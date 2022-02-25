################################################################################
# Analysis of all catchments for cryosphere-linkages
# 
# cryolink.R
#
# ReadMe: 
# Read and analyze discharge data from Langtang
#
#
#
# Created:          2022/02/05
# Latest Revision:  2022/02/05
#
#
# Jakob F Steiner| ICIMOD | jakob.steiner@icimod.org | x-hydrolab.org 
################################################################################
# clear entire workspace (excl. packages)
rm(list = ls())
gc()

# define &-sign for pasting string-elements
'&' <- function(...) UseMethod('&')
'&.default' <- .Primitive('&')
'&.character' <- function(...) paste(...,sep='')

library('geosphere')
library(leastcostpath)
library(raster)
library(rgdal)
library(rgeos)
library(geosphere)

projec_utm <- '+proj=utm +datum=WGS84'


data_path <- 'D:\\Work\\ICIMODProjects\\LivelihoodsCryosphere'
sitesFile <- 'Fieldsites_HKH.csv'

DEM_path <- 'D:\\Work\\GeospatialData\\HMA\\SRTM\\VoidFilled'
DEM_file <- 'SRTM_Corrected_Extended_HMA.tif'
DEM_HMA <- raster(DEM_path&'\\'&DEM_file)
#DEM_HMA <- projectRaster(DEM_HMA, crs=projec_utm)


snow_path <- 'D:\\Work\\GeospatialData\\HMA\\ICIMOD_modis_snow_eight_day_binary_2003_2020'

##########################
# Load RGI
##########################
path_RGI <- 'D:\\Work\\GeospatialData\\RGI60'                                     # Folder for RGI glacier outlines
RGI15_filename <- '15_rgi60_SouthAsiaEast.shp'
RGI14_filename <- '14_rgi60_SouthAsiaWest.shp'                                    # RGI filename
RGI13_filename <- '13_rgi60_CentralAsia.shp' 

ogrInfo(path_RGI&'\\'&RGI15_filename)
RGI60_15<-readOGR(dsn=path_RGI&'\\'&RGI15_filename)
#RGI60_15 <- spTransform(RGI60_15, projec)

ogrInfo(path_RGI&'\\'&RGI14_filename)
RGI60_14<-readOGR(dsn=path_RGI&'\\'&RGI14_filename)
#RGI60_14 <- spTransform(RGI60_14, projec)

ogrInfo(path_RGI&'\\'&RGI13_filename)
RGI60_13<-readOGR(dsn=path_RGI&'\\'&RGI13_filename)
#RGI60_13 <- spTransform(RGI60_13, projec)

##########################
# Load Permafrost Data
##########################
path_pf <- 'D:\\Work\\GeospatialData\\HMA\\Permafrost\\Obu2019\\UiO_PEX_PERPROB_5.0_20181128_2000_2016_NH' 
pf_filename <- 'OBU_PZI_HMA_WGS.tif'                                          # PZI filename

regionalPF <- raster(path_pf&'\\'&pf_filename)
#regionalPF <- projectRaster(regionalPF, crs = projec_utm)
regionalPF <- crop(regionalPF,extent(DEM_HMA))

##########################
# Load Precipitation Intensity Data
##########################
path_Precip <- 'D:\\Work\\ICIMODProjects\\TBWG\\LakeStandardization\\PPDIST\\' 
Precip_filename <- 'daily_intensity_point_scale.tif'                           # precipitation intensity with 10yr return period [mm d-1]

regionalPrecip <- raster(path_Precip&'\\'&Precip_filename)
#regionalPrecip <- crop(regionalPrecip,extent(regionalDEM))


outline_path <- 'D:\\Work\\ICIMODProjects\\LivelihoodsCryosphere\\Data\\Watersheds'

sitesDat <- read.csv(data_path&'\\'&sitesFile,header = T,sep=',')


##### Takmachik
outline_file <- paste(sitesDat$Site[1]&'_WS.shp')
ogrInfo(outline_path&'\\'&outline_file)
catch_outline<-readOGR(dsn=outline_path&'\\'&outline_file)
#catch_outline <- spTransform(catch_outline,CRSobj = projec_utm)

xy <- data.frame(X = c(sitesDat$Lon[1]), Y = c(sitesDat$Lat[1]))
coordinates(xy) <- ~ X + Y

glacMassLoss <- NA

elevSite <- extract(crop(DEM_HMA,catch_outline),xy)
elevTakmachik <- mask(crop(DEM_HMA,catch_outline),catch_outline)
areaTakmachik <- area(catch_outline)
PFTakmachik <- mask(crop(regionalPF,catch_outline),catch_outline)





##### Shishkat
outline_file <- paste(sitesDat$Site[2]&'_WS.shp')
ogrInfo(outline_path&'\\'&outline_file)
catch_outline<-readOGR(dsn=outline_path&'\\'&outline_file)
#catch_outline <- spTransform(catch_outline,CRSobj = projec_utm)

xy <- data.frame(X = c(sitesDat$Lon[2]), Y = c(sitesDat$Lat[2]))
coordinates(xy) <- ~ X + Y

glacMassLoss_1 <- raster('D:\\Work\\GeospatialData\\HMA\\GlacierChangeData\\Brun2017\\ASTER_shared_HMA\\ASTER_shared_HMA\\n37_e074\\dh_dt_2000-2016_ASTER_n37_e074.tif')
glacMassLoss_2 <- raster('D:\\Work\\GeospatialData\\HMA\\GlacierChangeData\\Brun2017\\ASTER_shared_HMA\\ASTER_shared_HMA\\n37_e075\\dh_dt_2000-2016_ASTER_n37_e075.tif')
glacMassLoss_3 <- raster('D:\\Work\\GeospatialData\\HMA\\GlacierChangeData\\Brun2017\\ASTER_shared_HMA\\ASTER_shared_HMA\\n36_e074\\dh_dt_2000-2016_ASTER_n36_e074.tif')
glacMassLoss_4 <- raster('D:\\Work\\GeospatialData\\HMA\\GlacierChangeData\\Brun2017\\ASTER_shared_HMA\\ASTER_shared_HMA\\n36_e075\\dh_dt_2000-2016_ASTER_n36_e075.tif')
glacMassLoss <- merge(glacMassLoss_1,glacMassLoss_2,glacMassLoss_3,glacMassLoss_4,tolerance=1)
glacMassLoss <- projectRaster(glacMassLoss,crs=projection(catch_outline))
glacMassLoss[abs(glacMassLoss[])>10] <- NA
glacMassLoss <- crop(glacMassLoss,catch_outline)

GlacSheeshkat <- rgeos::gIntersection(crop(RGI60_14,extent(catch_outline)),catch_outline)

glacMassLoss <- mask(glacMassLoss,GlacSheeshkat)
GLMasslossSheeshkat <- median( glacMassLoss[],na.rm=T )

elevSite <- extract(crop(DEM_HMA,catch_outline),xy)
elevSheeshkat <- mask(crop(DEM_HMA,catch_outline),catch_outline)
areaSheeshkat <- area(catch_outline)
PFSheeshkat <- mask(crop(regionalPF,catch_outline),catch_outline)



##### Lachung
outline_file <- paste(sitesDat$Site[3]&'_WS.shp')
ogrInfo(outline_path&'\\'&outline_file)
catch_outline<-readOGR(dsn=outline_path&'\\'&outline_file)
#catch_outline <- spTransform(catch_outline,CRSobj = projec_utm)

xy <- data.frame(X = c(sitesDat$Lon[2]), Y = c(sitesDat$Lat[2]))
coordinates(xy) <- ~ X + Y

glacMassLoss_1 <- raster('D:\\Work\\GeospatialData\\HMA\\GlacierChangeData\\Brun2017\\ASTER_shared_HMA\\ASTER_shared_HMA\\n37_e074\\dh_dt_2000-2016_ASTER_n37_e074.tif')
glacMassLoss_2 <- raster('D:\\Work\\GeospatialData\\HMA\\GlacierChangeData\\Brun2017\\ASTER_shared_HMA\\ASTER_shared_HMA\\n37_e075\\dh_dt_2000-2016_ASTER_n37_e075.tif')
glacMassLoss_3 <- raster('D:\\Work\\GeospatialData\\HMA\\GlacierChangeData\\Brun2017\\ASTER_shared_HMA\\ASTER_shared_HMA\\n36_e074\\dh_dt_2000-2016_ASTER_n36_e074.tif')
glacMassLoss_4 <- raster('D:\\Work\\GeospatialData\\HMA\\GlacierChangeData\\Brun2017\\ASTER_shared_HMA\\ASTER_shared_HMA\\n36_e075\\dh_dt_2000-2016_ASTER_n36_e075.tif')


elevSite <- extract(crop(DEM_HMA,catch_outline),xy)
elevSheeshkat <- mask(crop(DEM_HMA,catch_outline),catch_outline)
areaSheeshkat <- area(catch_outline)
PFSheeshkat <- mask(crop(regionalPF,catch_outline),catch_outline)

GlacSheeshkat <- rgeos::gIntersection(crop(RGI60_14,extent(catch_outline)),catch_outline)
mask(glacMassLoss,GlacSheeshkat)
 