library(ape)
library(rgdal)
library(rworldmap)
library(RColorBrewer)
library(scales)
library(geosphere)

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
d =  read.csv("~/Dropbox (Personal)/brazil_gene_flow/data/divergence.concat_ind0.05_loci0.6_all_n4796.csv",
              stringsAsFactors = F)

worldmap <- getMap(resolution = "coarse")
cols = brewer.pal(8, "Set2")

# get species data
x = read.csv("~/Dropbox (Personal)/brazil_gene_flow/data/brazil_samples_v8.csv",
             stringsAsFactors = F)
lins = table(x$lineage)
lins = names(lins[lins > 3])
lins = lins[which(lins != "Tupinambis_teguixin")]

# get phylogeny
# t = read.tree("~/Dropbox (Personal)/brazil_gene_flow/data/prelim_tree/ExaML_result.concat_ind0.05_loci0.6_all_n4796.rooted.no_outs.tre")
# cal = makeChronosCalib(t, node = "root", age.min = 193, age.max = 193)
# t = chronos(t, 0.001, calibration = cal)
# write.tree(t, "~/Dropbox (Personal)/brazil_gene_flow/data/prelim_tree/ExaML_result.concat_ind0.05_loci0.6_all_n4796.rooted.no_outs.dated.tre")
t = read.tree("~/Dropbox (Personal)/brazil_gene_flow/data/prelim_tree/ExaML_result.concat_ind0.05_loci0.6_all_n4796.rooted.no_outs.dated.tre")

# prepare maps
r = readOGR("~/Dropbox (Personal)/inornatus_gr/Meiri_etal_ranges/GARD1.1_dissolved_ranges/modeled_reptiles.shp")
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

# mtDNA tree
mt = read.tree("~/Dropbox (Personal)/brazil_gene_flow/data/mtDNA/mtDNA.concatenated.treefile.tre")
mt = drop.tip(mt, "CHUNB56712")
rootfams = c("Phyllodactylidae", "Gekkonidae")
rootinds = x[x$family %in% rootfams, "sample"]
rootinds = rootinds[ rootinds %in% mt$tip.label ]
rootinds = rootinds[ rootinds != "GRCOLLI23532"]
library(phytools)
node = findMRCA(mt, rootinds)
tips = na.omit(mt$tip.label[getDescendants(mt, node)])
mt2 = root(mt, tips)


targetlins = c("Oxyrhopus_trigeminus_1")
letter = "A"
for (i in 1:length(targetlins)) {
  lin = targetlins[i]
  
  genus = strsplit(lin, "_")[[1]][1]
  
  tsps = unique(x[which(x$lineage == lin), "species_checked"])
  tlins = unique(x[which(x$species_checked %in% tsps), "lineage"])
  x1 = x[which(x$species_checked %in% tsps | x$lineage %in% tlins), ]
  
  tsps = unique(x1$species_checked)
  tlins = unique(x1$lineage)
    
  # tt = getMRCA(t, intersect(x1$sample, t$tip.label))
  # tt = extract.clade(t, tt)
  tt = keep.tip(t, intersect(x1$sample, t$tip.label))
  
  # mtt = getMRCA(mt, intersect(x1$sample, mt$tip.label))
  # mtt = extract.clade(mt, mtt)
  mtt = keep.tip(mt2, intersect(x1$sample, mt$tip.label))
  
  png(paste0("~/Desktop/", lin, ".png"), width = 8, height = 4, units = "in", res = 300)
  layout.matrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2, ncol = 4)
  # layout.matrix <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
  layout(mat = layout.matrix)
  
  par(xpd = T, mar = c(2, 0.5, 0, 0))
  # print tree with species colors
  tt1 = tt
  tcols = cols[1:length(tsps)]
  names(tcols) = tsps
  sptips = x[match(tt1$tip.label, x$sample), "species_checked"]
  tt1$tip.label = paste(sptips, tt1$tip.label)
  tt1$tip.label = gsub(genus, 
                       paste0(strsplit(genus, "")[[1]][1], "."),
                       tt1$tip.label)
  plot(tt1, cex = 0.7, tip.color = tcols[sptips])
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  text(lastPP$x.lim[2] * 0.01, lastPP$y.lim[2] * 0.95, letter, font = 2)
  axisPhylo(cex = 0.5)
  
  # print tree with lin colors
  tt2 = tt
  tcols = cols[1:length(tlins)]
  names(tcols) = tlins
  lintips = x[match(tt2$tip.label, x$sample), "lineage"]
  tt2$tip.label = paste(lintips, tt2$tip.label)
  tt2$tip.label = gsub(genus, 
                       paste0(strsplit(genus, "")[[1]][1], "."),
                       tt2$tip.label)
  plot(tt2, cex = 0.7, tip.color = tcols[lintips])
  axisPhylo(cex = 0.5)
  
  par(xpd = T, mar = c(1, 0, 0, 0))
  # print mt tree with species colors
  tt1 = mtt
  scols = cols[1:length(tsps)]
  names(scols) = tsps
  sptips = x[match(tt1$tip.label, x$sample), "species_checked"]
  tt1$tip.label = paste(sptips, tt1$tip.label)
  tt1$tip.label = gsub(genus, 
                       paste0(strsplit(genus, "")[[1]][1], "."),
                       tt1$tip.label)
  plot(tt1, cex = 0.7, tip.color = scols[sptips])
  
  # print mt tree with lin colors
  tt2 = mtt
  tcols = cols[1:length(tlins)]
  names(tcols) = tlins
  lintips = x[match(tt2$tip.label, x$sample), "lineage"]
  tt2$tip.label = paste(lintips, tt2$tip.label)
  tt2$tip.label = gsub(genus, 
                       paste0(strsplit(genus, "")[[1]][1], "."),
                       tt2$tip.label)
  plot(tt2, cex = 0.7, tip.color = tcols[lintips])
  
  # plot brazil
  par(mar = c(1, 1, 1, 1), xpd = F)
  plot(worldmap, col = "lightgrey", 
       border = "lightgrey",
       xlim = c(-90, -30), ylim = c(-22, 1),
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
  points(x1$LON, x1$LAT, pch = 21, bg = scols[x1$species_checked])
  

  # plot brazil
  plot(worldmap, col = "lightgrey", 
       border = "lightgrey",
       xlim = c(-90, -30), ylim = c(-22, 1),
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
  
  sp = gsub("_\\d+", "", lin)
  x2 = x[which(x$species_checked == sp), ]
  dist = data.frame(t(combn(x2$sample, 2)), stringsAsFactors = F)
  names(dist) = c("sample1", "sample2")
  dist$geo_dist = apply(dist, 1, get_dist) / 1000
  vals = apply(dist, 1, get_gen_dist)
  dist$diff = vals[1, ]
  dist$denom = vals[2, ]
  dist$div = dist$diff / dist$denom
  # dist$sp1 = x1[match(dist$sample1, x1$sample), "lineage"]
  # dist$sp2 = x1[match(dist$sample2, x1$sample), "lineage"]
  # dist$same = ifelse(dist$sp1 == dist$sp2, TRUE, FALSE)
  ylim2 = range(dist$div)
   
  par(mar = c(4, 4, 0.5, 0.5))
  plot(dist$geo_dist, dist$div, las = 1,
       pch = 21, bg = "gray", ylim = ylim2,
       xlab = "geo. dist., km", ylab = "nucl. dxy")
  
  x3 = x[which(x$lineage == lin), ]
  dist = data.frame(t(combn(x3$sample, 2)), stringsAsFactors = F)
  names(dist) = c("sample1", "sample2")
  dist$geo_dist = apply(dist, 1, get_dist) / 1000
  vals = apply(dist, 1, get_gen_dist)
  dist$diff = vals[1, ]
  dist$denom = vals[2, ]
  dist$div = dist$diff / dist$denom
  # dist$sp1 = x1[match(dist$sample1, x1$sample), "lineage"]
  # dist$sp2 = x1[match(dist$sample2, x1$sample), "lineage"]
  # dist$same = ifelse(dist$sp1 == dist$sp2, TRUE, FALSE)
  
  plot(dist$geo_dist, dist$div, las = 1,
       pch = 21, bg = "gray", ylim = ylim2,
       xlab = "geo. dist., km", ylab = "nucl. dxy")
  
  dev.off()
}
