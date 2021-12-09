library(ape)
library(rgdal)
library(rworldmap)
library(RColorBrewer)
library(scales)

get_dist <- function(row) {
  pt1 = x[which(x$sample == row[1]), c("LON", "LAT")]
  pt2 = x[which(x$sample == row[2]), c("LON", "LAT")]
  return(distHaversine(pt1, pt2))
}

get_gen_dist <- function(row) {
  ind1 = as.character(row[1])
  ind2 = as.character(row[2])
  
  d1 = d[d$ind1 == ind1 & d$ind2 == ind2, ]
  d2 = d[d$ind1 == ind2 & d$ind2 == ind1, ]
  dd = rbind(d1, d2)
  
  if (nrow(dd) < 1) {
    return(c(NA, NA))
  } else {
    return(c(dd$diff, dd$denom))
  }
}

# get divergence
d =  read.csv("~/Dropbox/brazil_gene_flow/data/divergence.concat_ind0.05_loci0.6_all_n4796.csv",
              stringsAsFactors = F)

worldmap <- getMap(resolution = "coarse")
cols = brewer.pal(8, "Dark2")

# get species data
x = read.csv("~/Dropbox/brazil_gene_flow/data/brazil_samples_v7.csv",
             stringsAsFactors = F)
lins = table(x$lineage)
lins = names(lins[lins > 3])
lins = lins[which(lins != "Tupinambis_teguixin")]

# get phylogeny
t = read.tree("~/Dropbox/brazil_gene_flow/data/prelim_tree/ExaML_result.concat_ind0.05_loci0.6_all_n4796.rooted.no_outs.tre")

# mtDNA tree
mt = read.tree("~/Dropbox/brazil_gene_flow/data/mtDNA/mtDNA.concatenated.treefile.tre")

# prepare maps
r = readOGR("~/Dropbox/inornatus_gr/Meiri_etal_ranges/GARD1.1_dissolved_ranges/modeled_reptiles.shp")
rnames = r[[1]]
rnames = gsub(" ", "_", rnames)

spnames = lins
names(spnames) = spnames
spnames[! spnames %in% rnames]

spnames["Xenodon_merremii"] = "Xenodon_merremi" 
spnames["Copeoglossum_nigropunctatum_1"] = "Copeoglossum_nigropunctatum" 
spnames["Bachia_bresslaui_1"] = "Bachia_bresslaui" 
spnames["Oxyrhopus_trigeminus_1"] = "Oxyrhopus_trigeminus" 
spnames["Oxyrhopus_trigeminus_2"] = "Oxyrhopus_trigeminus" 
spnames["Taeniophallus_occipitalis_2"] = "Taeniophallus_occipitalis" 
spnames["Trilepida_brasiliensis_2"] = "Trilepida_brasiliensis" 
spnames["Vanzosaura_rubricauda_2"] = "Vanzosaura_rubricauda" 


for (i in 1:length(lins)) {
  lin = lins[i]
  tsps = unique(x[which(x$lineage == lin), "species"])
  tlins = unique(x[which(x$species %in% tsps), "lineage"])
  x1 = x[which(x$species %in% tsps | x$lineage %in% tlins), ]
  
  # tt = getMRCA(t, intersect(x1$sample, t$tip.label))
  # tt = extract.clade(t, tt)
  tt = keep.tip(t, intersect(x1$sample, t$tip.label))
  
  # mtt = getMRCA(mt, intersect(x1$sample, mt$tip.label))
  # mtt = extract.clade(mt, mtt)
  mtt = keep.tip(mt, intersect(x1$sample, mt$tip.label))
  
  pdf(paste0("~/Desktop/species/", lin, ".pdf"), width = 10, height = 6)
  par(xpd = T, mar = c(0, 0, 0, 4))
  layout.matrix <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
  layout(mat = layout.matrix)
  
  # print tree with species colors
  tt1 = tt
  tcols = cols[1:length(tsps)]
  names(tcols) = tsps
  sptips = x[match(tt1$tip.label, x$sample), "species"]
  tt1$tip.label = paste(sptips, tt1$tip.label)
  plot(tt1, cex = 0.7, tip.color = tcols[sptips])
  
  # print tree with lin colors
  tt2 = tt
  tcols = cols[1:length(tlins)]
  names(tcols) = tlins
  lintips = x[match(tt2$tip.label, x$sample), "lineage"]
  tt2$tip.label = paste(lintips, tt2$tip.label)
  plot(tt2, cex = 0.7, tip.color = tcols[lintips])
  
  # print mt tree with species colors
  tt1 = mtt
  tcols = cols[1:length(tsps)]
  names(tcols) = tsps
  sptips = x[match(tt1$tip.label, x$sample), "species"]
  tt1$tip.label = paste(sptips, tt1$tip.label)
  plot(tt1, cex = 0.7, tip.color = tcols[sptips])
  
  # print mt tree with lin colors
  tt2 = mtt
  tcols = cols[1:length(tlins)]
  names(tcols) = tlins
  lintips = x[match(tt2$tip.label, x$sample), "lineage"]
  tt2$tip.label = paste(lintips, tt2$tip.label)
  plot(tt2, cex = 0.7, tip.color = tcols[lintips])
  
  # plot brazil
  par(mar = c(2, 2, 2, 2), xpd = F)
  plot(worldmap, col = "lightgrey", 
       border = "lightgrey",
       xlim = c(-100, -30), ylim = c(-22, 12),
       bg = "white",
       asp = 1)
  
  # plot range
  rnum = which(rnames == spnames[lin])
  if (length(rnum) > 0) {
    range = r[rnum, 1]
    plot(range, add = T, col = alpha("black", 0.5),
         border = F)
  }
  
  # add points
  points(x1$LON, x1$LAT, pch = 21, bg = tcols[x1$lineage])
  
  x2 = x[which(x$lineage == lin), ]
  dist = data.frame(t(combn(x2$sample, 2)), stringsAsFactors = F)
  names(dist) = c("sample1", "sample2")
  dist$geo_dist = apply(dist, 1, get_dist) / 1000
  vals = apply(dist, 1, get_gen_dist)
  dist$diff = vals[1, ]
  dist$denom = vals[2, ]
  dist$div = dist$diff / dist$denom
  
  plot(dist$geo_dist, dist$div, las = 1, pch = 21, bg ="forestgreen")
  
  dev.off()
}
