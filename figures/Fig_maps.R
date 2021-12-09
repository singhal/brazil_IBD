library(rgdal)
library(scales)

s = readRDS("~/Dropbox (Personal)/brazil_gene_flow/data/sample_map.Rds")
ll = read.csv("~/Dropbox (Personal)/brazil_gene_flow/data/brazil_samples_v8.csv",
              stringsAsFactors = F)
ll2 = ll[ll$lineage %in% s$lineage, ]


# Test beta_IBD vs range size correlation via PGLS
r = readOGR("~/Dropbox (Personal)/inornatus_gr/Meiri_etal_ranges/GARD1.1_dissolved_ranges/modeled_reptiles.shp")
rnames = r[[1]]
rnames = gsub(" ", "_", rnames)

spnames = unique(s$lineage)
names(spnames) = spnames
spnames[! spnames %in% rnames]

spnames["Xenodon_merremii"] = "Xenodon_merremi"
spnames["Trilepida_brasiliensis_2"] = "Trilepida_brasiliensis"
spnames["Taeniophallus_occipitalis_2"] = "Taeniophallus_occipitalis"
spnames["Oxyrhopus_trigeminus_2"] = "Oxyrhopus_trigeminus"
spnames["Spilotes_sp"] = "Spilotes_sulphureus"
spnames["Copeoglossum_nigropunctatum_1"] = "Copeoglossum_nigropunctatum"
spnames["Oxyrhopus_trigeminus_1"] = "Oxyrhopus_trigeminus"
spnames["Bachia_bresslaui_1"] = "Bachia_bresslaui"

worldmap <- rworldmap::getMap(resolution = "li")
sa = worldmap[which(worldmap[["continent"]] == "South America"), ]

pdf("~/Desktop/all_range_maps.pdf", width = 12, height = 27)
par(mfrow = c(10, 6))
for (i in 1:59) {
  par(mar = c(0, 0, 0, 0))
  plot(NA, xlim = c(-85, -34), ylim = c(-56, 13), 
       axes = F, xlab = "", ylab = "")
  plot(sa, col = "gray85", border = F, add = T)
  
  plot(r[ which(rnames == spnames[i]), ], add = T, 
       col = alpha("forestgreen", 0.5), border = F) 
  
  pts = ll[which(ll$lineage == names(spnames)[i]), c("LON", "LAT")]
  points(pts, pch = 21, bg = "white", cex = 2)
  text(x = -59, y = -50, labels = gsub("_", " ", names(spnames)[i]), font = 3, cex = 1.5)
}
dev.off()