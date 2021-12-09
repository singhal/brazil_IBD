library(ape)
library(adephylo)
library(RColorBrewer)
library(scales)
library(dplyr)
library(ggplot2)
library(geosphere)
library(cowplot)
library(vegan)
library(tidyverse)
theme_set(theme_cowplot())

get_dist <- function(row) {
  pt1 = x1[which(x1$sample == row[2]), c("LON", "LAT")]
  pt2 = x1[which(x1$sample == row[3]), c("LON", "LAT")]
  return(distHaversine(pt1, pt2))
}

get_climdist <- function(row) {
  return(climdist[ row[2], row[3] ])
}


inv_fst <- function(x) {
  inv = x / (1 - x)
  return(inv)
}

get_mt_dist <- function(dist) {
  inds = unique(c(dist$ind1, dist$ind2))
  indix = match(inds, dimnames(mt)[[1]])
  indix = indix[!is.na(indix)]
  mt2 = mt[indix, ] 
  mt3 = dist.dna(mt2, model = "raw", 
                 pairwise.deletion = T, as.matrix = T)
  mtinds = dimnames(mt2)[[1]]
  for (j in 1:nrow(dist)) {
    ind1 = dist[j, "ind1"]
    ind2 = dist[j, "ind2"]
    
    seq1 = as.character(mt2[match(ind1, mtinds), ])
    seq2 = as.character(mt2[match(ind2, mtinds), ])
    seqs = rbind(seq1, seq2)
    miss = apply(seqs, 2, function(x) {sum(x == '-')})
    denom = length(which( miss < 1))
    
    if (ind1 %in% mtinds & ind2 %in% mtinds) {
      dist[j, "mt_dist"] = mt3[ind1, ind2]
      dist[j, "mt_tree_dist"] = mtdist2[ind1, ind2]
      dist[j, "mt_denom"] = denom
    } else {
      dist[j, "mt_dist"] = NA
      dist[j, "mt_tree_dist"] = NA
      dist[j, "mt_denom"] = NA
    }
    
  }
  return(dist)
}

get_mt_dist2 <- function(dist) {
  
  mtinds = dimnames(mt)[[1]]

  for (j in 1:nrow(dist)) {
    ind1 = dist[j, "ind1"]
    ind2 = dist[j, "ind2"]
    
    if (ind1 %in% mtinds & ind2 %in% mtinds) {
      seq1 = as.character(mt[match(ind1, mtinds), ])
      seq2 = as.character(mt[match(ind2, mtinds), ])
      seqs = rbind(seq1, seq2)
      miss = apply(seqs, 2, function(x) {sum(x == '-')})
      seqs2 = seqs[ ,which(miss < 1)]
      
      denom = length(which(miss < 1))
      if (denom > 50) {
        mtdiff = apply(seqs2, 2, function(x) {x[1] == x[2]})
        
        dist[j, "mt_dist"] = sum(mtdiff) / denom
        dist[j, "mt_len"] = denom
      } else {
        dist[j, "mt_dist"] = NA
        dist[j, "mt_len"] = NA        
      }
    } else {
      dist[j, "mt_dist"] = NA
      dist[j, "mt_len"] = NA
    }
    
  }
  return(dist)
}

get_slope <- function(divtype, disttype) {
  dist2 = dist[ ,c("ind1", "ind2", divtype, disttype)]
  dist2 <- dist2[is.finite(rowSums(dist2[, c(divtype, disttype)])), ]
  num = nrow(dist2)
  
  dist3 = dist2 %>% select(ind2 = ind1, 
                           ind1 = ind2, 
                           all_of(divtype),
                           all_of(disttype))
  dist4 = rbind(dist2, dist3)
  
  a = dist4 %>% select(-all_of(disttype)) %>%
    spread(ind1, all_of(divtype), fill=0) %>% 
    column_to_rownames(var="ind2") %>% 
    as.matrix
  b = dist4 %>% select(-all_of(divtype)) %>%
    spread(ind1, all_of(disttype), fill=0) %>% 
    column_to_rownames(var="ind2") %>% 
    as.matrix
  
  mt = mantel(a, b)
  sig_mt = mt$signif
  r2_mt = mt$statistic ^ 2
  
  m = lm(dist2[, divtype] ~ dist2[, disttype])
  slope = m$coefficients[2]
  intercept = m$coefficients[1]
  r2 = summary(m)$r.squared
  sig = summary(m)$coefficients[2, 4]
  se = summary(m)$coefficients[2, 2]
  
  gg = ggplot(dist2, aes_string(disttype, divtype)) + 
    geom_point() + 
    stat_smooth(method = "lm", col = "black") +
    ggtitle(as.expression(bquote(r^2~"="~.(round(r2_mt, 2)) ~ "; p="~ .(round(sig_mt, 2)))))
  
  res = c(divtype, disttype, slope, intercept, se, num, r2, sig, r2_mt, sig_mt)
  return(list(res, gg))
}

# get samples
x1 = read.csv("~/Dropbox (Personal)/brazil_gene_flow/data/brazil_samples_v8.csv",
             stringsAsFactors = F)
x1[x1$lineage == "Vanzosaura_savanicola", "lineage"] = "Vanzosaura_rubricauda_2"
div = list.files("~/Dropbox (Personal)/brazil_gene_flow/data/pop_gen/",
                 pattern = "csv", full.names = T)
div = div[grep("teguixin", div, invert = T)]

# mt dna
mt = read.dna("~/Dropbox (Personal)/brazil_gene_flow/data/mtDNA/concatenated.fasta",format = "fasta")
mttree = read.tree("~/Dropbox (Personal)/brazil_gene_flow/data/mtDNA/mtDNA.concatenated.treefile.tre")
mtdist = distTips(mttree)
mtdist2 = as.matrix(mtdist)

# what to output
# 1. all divergences across species
# 2. all slopes: # comp, val, sig, r2
# 3. figure of all three
res1 = vector("list", length(dist))
res2 = vector("list", length(dist))

# climate distances
clim = list.files("~/Dropbox (Personal)/brazil_gene_flow/data/climate_variables/",
                  full.names = T)
clim2 = raster::stack(clim)
lins = lapply(div, read.csv, stringsAsFactors = F)
lins = unlist(lapply(lins, function(x) {x$lineage[1]}))
x2 = x1[x1$lineage %in% lins, ]
climpts = raster::extract(clim2, x2[, c("LON", "LAT")])
rownames(climpts) = x2$sample
climpca = prcomp(climpts, scale. = T, center = T)
climpca2 = climpca$x
climdist = as.matrix(dist(climpca2, diag = T, upper = T))

for (i in 1:length(div)) {

  dist = read.csv(div[i], stringsAsFactors = F)
  dist$geo_dist = apply(dist, 1, get_dist) / 1000
  dist$ln_geo_dist = log(dist$geo_dist)
  dist$inv_fst = inv_fst(dist$fst)
  dist = get_mt_dist(dist)
  dist$clim_dist = apply(dist, 1, get_climdist)
  
  res1[[i]] = dist 
  lineage = dist$lineage[1]
  
  a = get_slope('inv_fst', 'ln_geo_dist')
  b = get_slope('nuc_dxy', 'geo_dist')
  d = get_slope('fst', 'clim_dist')
  dist = dist[which(dist$mt_denom > 500), ]
  if (nrow(dist) > 1) {
    c = get_slope('mt_dist', 'geo_dist')
  } else {
    c = list(c('mt_dist', 'geo_dist', NA, NA, NA, NA, NA, NA),
             ggplot() + theme_void())
  }
  dd = data.frame(rbind(a[[1]], b[[1]], c[[1]], d[[1]]))
  dd$lineage = lineage
  gr = plot_grid(a[[2]], b[[2]], c[[2]], d[[2]], ncol = 2)
  # save_plot(paste0("~/Dropbox/brazil_gene_flow/results/IBD_plots/", lineage,
  #                  ".pdf"), gr, base_width = 4, base_height = 3, 
  #           nrow = 2, ncol = 2)
  res2[[i]] = dd
}


res1a = do.call("rbind", res1)
res2a = do.call("rbind", res2)
names(res2a) = c("div_type", "dist_type",
                 "slope", "intercept", "se", "n", "r2", "sig",
                 "r2_mt", "sig_mt",
                 "lineage")
for (i in c("slope", "n", "r2", "sig", "r2_mt", "sig_mt")) {
  res2a[, i] = as.numeric(res2a[, i])
}

write.csv(res1a, paste0("~/Dropbox (Personal)/brazil_gene_flow/results/divergence-",
                        Sys.Date(), ".csv"), quote = F, row.names = F)
write.csv(res2a, paste0("~/Dropbox (Personal)/brazil_gene_flow/results/IBD-",
                        Sys.Date(), ".csv"), quote = F, row.names = F)
