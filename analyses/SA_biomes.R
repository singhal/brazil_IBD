library(sf)
library(raster)
library("tidyverse")

eco = read_sf("~/Dropbox/macroevolution/eco_IBD_oz/data/ecoregions/terr-ecoregions-TNC/tnc_terr_ecoregions.shp")

worldmap <- rworldmap::getMap(resolution = "low")
sa = worldmap[which(worldmap[["continent"]] == "South America"), ]
sa2 = crop(sa, extent(-85, -34, -56, 13))

insa = c()
econames = eco$ECO_NAME
for (i in 1:length(econames)) {
  r = eco[eco$ECO_NAME == econames[i], ]
  er = extent(r)
  
  eri = intersect(extent(er), extent(sa2))
  if (!is.null(eri)) {
    insa = c(insa, econames[i])
  }
}
eco2 = eco[eco$ECO_NAME %in% insa, ]

eco3 <-
  eco2 %>%
  group_by(WWF_MHTNAM) %>%
  summarise()

write_sf(eco3, '~/Desktop/SA_biomes.shp')
