library(Matrix)
library(tidyverse)
library(sf)

### All the needed shapefiles
US_ecorgns=st_read("GIS_Data/us_eco_l3_state_boundaries/us_eco_l3_state_boundaries.shp") 
US_shp=st_read("GIS_Data/gadm41_USA_shp/gadm41_USA_1.shp")
US_HUC8=st_read("GIS_Data/HUC8_CONUS/HUC8_US.shp")

### Subset ecoregions shapefile to select two key ecoregions:
US_ecorgns.sel=US_ecorgns%>%filter(US_L3CODE==50 | US_L3CODE==51)
US_ecorgns.sel.UTM=st_transform(US_ecorgns.sel,crs="+proj=utm +zone=15 +datum=WGS84 +units=m +no_defs")


### Subset the state of MN shape file
US_shp
MN_shp=US_shp%>%filter(NAME_1 %in% "Minnesota")
MN_shp.UTM=st_transform(MN_shp, crs="+proj=utm +zone=15 +datum=WGS84 +units=m +no_defs")
MN_ecorgns=st_intersection(MN_shp.UTM,US_ecorgns.sel.UTM)
ggplot(MN_shp.UTM)+geom_sf()+
  geom_sf(data=MN_ecorgns, fill="green")

###Final map of MN ecoregions
MN_ecorgns.unqL3=MN_ecorgns%>%distinct(US_L3CODE, .keep_all = TRUE)
ggplot(MN_shp.UTM)+geom_sf()+
geom_sf(data=MN_ecorgns, aes(fill=US_L3CODE))+geom_sf_text(data=MN_ecorgns.unqL3, aes(label=US_L3NAME))+
theme(legend.position = "NA")+xlab(" ")+ylab(" ")+theme(text=element_text(size=18))+ggtitle("Level-3 Ecoregions")
ggsave("Fig1.png", path="Figures/", device="png",width=9, height=12, dpi=900)


### Subset HUC8 by selecting only watersheds within the selected ECOREGIONS
US_HUC8
US_HUC8.UTM=st_transform(US_HUC8,crs="+proj=utm +zone=15 +datum=WGS84 +units=m +no_defs")

### overlapping polygons is tricky
### so converted HUC8 to points by estimating their centroids
US_HUC8_centr=st_centroid(US_HUC8.UTM)
HUC8.CENTR_IN_SEL.ECORGNS=st_join(US_HUC8_centr,US_ecorgns.sel.UTM, join=st_within)
HUC8.CENTR_IN_SEL.ECORGNS
### and then, spatially joined with selected ecoregions within MN
HUC8.CENTR_IN_SEL.ECORGNS2=HUC8.CENTR_IN_SEL.ECORGNS%>%filter(US_L3CODE !="NA")
HUC8.CENTR_IN_SEL.ECORGNS2
HUC8.CENTR_IN_MN_SEL.ECORGNS=st_join(HUC8.CENTR_IN_SEL.ECORGNS2,MN_shp.UTM, join=st_within)
HUC8.CENTR_IN_MN_SEL.ECORGNS2=HUC8.CENTR_IN_MN_SEL.ECORGNS%>%filter(NAME_1 != "NA")
### and finally, join the selected HUC8 centroids to HUC8 watershed shapefile to get the final subset of waterhsed
### shapefiles for further analyses
MN_allHUC8=HUC8.CENTR_IN_MN_SEL.ECORGNS2%>%pull(HUC8)
MN_allHUC8

HUC8.poly_MN_SEL.ECORGNS=US_HUC8.UTM%>%filter(HUC8 %in% MN_allHUC8)
HUC8.poly_MN_SEL.ECORGNS_bndry.corr=st_intersection(MN_shp.UTM,HUC8.poly_MN_SEL.ECORGNS)
st_write(HUC8.poly_MN_SEL.ECORGNS,"ProcessedData/HUC8.poly_MN_Sel.EcoRgns.shp") ### GIS FILE

### a quick visualzation of all the above steps
ggplot(MN_shp.UTM)+geom_sf()+
geom_sf(data=HUC8.poly_MN_SEL.ECORGNS, fill="green")+
geom_sf(data=HUC8.CENTR_IN_MN_SEL.ECORGNS2)+geom_sf_text(data=HUC8.CENTR_IN_MN_SEL.ECORGNS2, aes(label=HUC8))
###Final figure
ggplot(MN_shp.UTM)+geom_sf()+geom_sf(data=MN_ecorgns, col="blue", alpha=0.3)+
  geom_sf(data=HUC8.poly_MN_SEL.ECORGNS_bndry.corr, fill="green", alpha=0.6)+
    geom_sf(data=HUC8.CENTR_IN_MN_SEL.ECORGNS2)+ggtitle("Watersheds in Level-3 Ecoregions")+
      theme(text=element_text(size=18))
  ggsave("Fig2.png", path="Figures/", device="png",width=9, height=12, dpi=900)

#######################################################################################################################
#### Transform to Lat Lon, and convert the HUC8 shapefile to a csv file with lat lon
HUC8.CENTR_IN_MN_SEL.ECORGNS2_wgs=st_transform(HUC8.CENTR_IN_MN_SEL.ECORGNS2, crs=4326)
HUC8.CENTR_IN_MN_SEL.ECORGNS2_wgs
MN_HUC8.WTRSHDS= HUC8.CENTR_IN_MN_SEL.ECORGNS2_wgs%>%
    mutate(Lon = unlist(map(HUC8.CENTR_IN_MN_SEL.ECORGNS2_wgs$geometry,1)),
    Lat = unlist(map(HUC8.CENTR_IN_MN_SEL.ECORGNS2_wgs$geometry,2)))%>%select(c(HUC8,Lon, Lat))%>%
      st_drop_geometry()%>%arrange(Lat)
MN_HUC8.WTRSHDS
write_csv(MN_HUC8.WTRSHDS, "ProcessedData/All.MN_HUC8.wtrshds.csv")
MN_HUC8.WTRSHDS=read_csv("ProcessedData/All.MN_HUC8.wtrshds.csv")

### Also save all the centroid point shapefiles
st_write(HUC8.CENTR_IN_MN_SEL.ECORGNS2, "GIS_Data/HUC8_MN_centr.shp", delete_dsn = TRUE) ###GIS FILE

#######################################################################################################################
#######################################################################################################################




