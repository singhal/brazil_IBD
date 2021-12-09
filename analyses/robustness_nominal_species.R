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

# get species data
x = read.csv("~/Dropbox (Personal)/brazil_gene_flow/data/brazil_samples_v7.csv",
             stringsAsFactors = F)
x = x[complete.cases(x$LAT), ]
lins = table(x$species_checked)
lins = names(lins[lins > 3])
lins = lins[which(lins != "Tupinambis_teguixin")]
lins = lins[which(lins != "Phrynonax_poecilonotus")]

res = vector("list", length(lins))
for (i in 1:length(lins)) {
  x2 = x[which(x$species_checked == lins[i]), ]
  dist = data.frame(t(combn(x2$sample, 2)), stringsAsFactors = F)
  names(dist) = c("sample1", "sample2")
  dist$geo_dist = apply(dist, 1, get_dist) / 1000
  vals = apply(dist, 1, get_gen_dist)
  dist$diff = vals[1, ]
  dist$denom = vals[2, ]
  dist$div = dist$diff / dist$denom
  
  a = lm(dist$div ~ dist$geo_dist)
  slope = coef(a)[2]
  rsq = summary(a)$adj.r.squared
  pval = summary(a)$coefficients[2, 4]
  n = nrow(x2)
  
  res[[i]] = c(lins[i], n, slope, pval, rsq)
}

res2 = as.data.frame(do.call("rbind", res))
names(res2) = c("species", "n", "slope", "pval", "rsq")
write.csv(res2, "~/Dropbox (Personal)/brazil_gene_flow/results/nominal_species_IBD-6-Apr-2021.csv",
          quote = F, row.names = F)
