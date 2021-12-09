
eco = readOGR("~/Dropbox (Personal)/macroevolution/eco_IBD_oz/data/ecoregions/terr-ecoregions-TNC/tnc_terr_ecoregions.shp")

x1 = read.csv("~/Dropbox (Personal)/brazil_gene_flow/data/brazil_samples_v8.csv",
              stringsAsFactors = F)
lins = table(x1$lineage)
lins = lins[lins > 3]
lins = lins[which(!names(lins) %in% c("Phrynonax_poecilonotus", "Tupinambis_teguixin"))]
x1 = x1[x1$lineage %in% names(lins), ]

pts = SpatialPoints(x1[,c("LON", "LAT")])
pts2 = raster::intersect(pts, eco)
eco2 = table(pts2$ECO_NAME)

worldmap <- rworldmap::getMap(resolution = "low")
sa = worldmap[which(worldmap[["continent"]] == "South America"), ]

# drop 
### lon < -62
### lat > 0
x1$keep = TRUE
x1[which(x1$LON < -62), "keep"] = FALSE
x1[which(x1$LAT > 0), "keep"] = FALSE
table(x1$keep)

par(mar = c(0,0, 0,0))
plot(NA, xlim = c(-85, -34), ylim = c(-56, 13), 
     axes = F, xlab = "", ylab = "")
plot(sa, col = "gray90", add = T, border = F)
points(x1[,c("LON", "LAT")], col = as.factor(x1$keep), pch = 16)

r = eco[eco$ECO_NAME == "Chiquitano Dry Forests", ]
plot(r, col = "red", border = F, add = T)

write.csv(x1, "~/Dropbox (Personal)/brazil_gene_flow/data/brazil_samples_v8_nonCerrado.csv",
          row.names = F)

