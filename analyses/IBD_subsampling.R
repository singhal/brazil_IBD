library(ape)
library(adephylo)
library(RColorBrewer)
library(scales)
library(dplyr)
library(rgeos)
library(ggplot2)
library(geosphere)
library(cowplot)
library(vegan)
library(tidyverse)
library(alphahull)
theme_set(theme_cowplot())

source("~/Dropbox (Personal)/scripts/eco_IBD/pascal_scripts/ah2sp.R")

get_dist <- function(row) {
  pt1 = x1[which(x1$sample == row[2]), c("LON", "LAT")]
  pt2 = x1[which(x1$sample == row[3]), c("LON", "LAT")]
  return(distHaversine(pt1, pt2))
}

inv_fst <- function(x) {
  inv = x / (1 - x)
  return(inv)
}

get_slope <- function(dist, divtype, disttype) {
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
  if (nrow(dist2) > 2) {
    m = lm(dist2[, divtype] ~ dist2[, disttype])
    slope = m$coefficients[2]
    intercept = m$coefficients[1]
    r2 = summary(m)$r.squared
    sig = summary(m)$coefficients[2, 4]
    se = summary(m)$coefficients[2, 2]
    
    res = c(divtype, disttype, slope, intercept, se, num, r2, sig, r2_mt, sig_mt)
  } else {
    res = rep(NA, 10)
  }
  return(res)
}

# get samples
x1 = read.csv("~/Dropbox (Personal)/brazil_gene_flow/data/brazil_samples_v8.csv",
              stringsAsFactors = F)
x1[x1$lineage == "Vanzosaura_savanicola", "lineage"] = "Vanzosaura_rubricauda_2"
div = list.files("~/Dropbox (Personal)/brazil_gene_flow/data/pop_gen/",
                 pattern = "csv", full.names = T)
div = div[grep("teguixin", div, invert = T)]

all = vector("list", length(div))
for (i in 1:length(div)) {
  cat(i, "\n")
  dist = read.csv(div[i], stringsAsFactors = F)
  dist$geo_dist = apply(dist, 1, get_dist) / 1000
  dist$ln_geo_dist = log(dist$geo_dist)
  dist$inv_fst = inv_fst(dist$fst)
  
  lineage = dist$lineage[1]
  inds = unique(c(dist$ind1, dist$ind2))
  
  # pts = x1[match(inds, x1$sample), c("LON", "LAT")]
  # range = getAlphaHullRange(pts, percent=1, partCount = 3,
  #                            buff=50, 
  #                            coordHeaders=c('LON','LAT'))
  # area = gArea(range[[1]])
  
  if (length(inds) > 5) {
    combos = combn(inds, 5)
    res = vector("list", ncol(combos))
    for (j in 1:ncol(combos)) {
      inds2 = combos[, j]
      
      pts2 = x1[match(inds2, x1$sample), c("LON", "LAT")]
      pts3 = unique(pts2)
      # range2 = getAlphaHullRange(pts2, percent=1, partCount = 3,
      #                          buff=50, 
      #                          coordHeaders=c('LON','LAT'))
      # area2 = gArea(range2[[1]])
      if(nrow(pts3) >= 5) {
        # need to check if these are actually unique localities
        dist2 = dist[dist$ind1 %in% inds2, ]
        dist2 = dist2[dist2$ind2 %in% inds2, ]
        res[[j]] = get_slope(dist2, 'inv_fst', 'ln_geo_dist')
      } else {
        res[[j]] = rep(NA, 10)
      }
    }
    res2 = as.data.frame(do.call("rbind", res))
    res2$lineage = lineage
    names(res2) = c("div_type", "dist_type",
                     "slope", "intercept", "se", "n", "r2", "sig",
                     "r2_mt", "sig_mt",
                     "lineage")
    all[[i]] = res2
  } else {
    all[[i]] = NULL
  }
}

all2 = all[-which(sapply(all, is.null))]
all3 = do.call("rbind", all2)
vars = c("slope", "intercept", "se", "n",
         "r2", "sig", "r2_mt", "sig_mt")
for (var in vars) {
  all3[ , var] = as.numeric(all3[, var])
}
all3$log_slope = log(all3$slope)
all4 = all3[complete.cases(all3$log_slope), ]

xx = read.csv("~/Dropbox (Personal)/brazil_gene_flow/results/IBD-2021-09-15.csv")
xx1 = xx[xx$div_type == "inv_fst", ]
xx1$log_slope = log(xx1$slope)

dd = left_join(all4, xx1, by = c("lineage" = "lineage"))
dd2 = dd %>% group_by(lineage) %>% mutate(comps = n()) %>% 
  ungroup() %>% filter(comps > 10) %>% 
  filter(lineage != "Oxyrhopus_trigeminus_1", lineage != "Trilepida_koppesi")
dd2$lineage2 = gsub("_", " ", dd2$lineage)
ab = ggplot(dd2, aes(log_slope.x)) + geom_histogram(fill = "gray60") +
  geom_vline(aes(xintercept = log_slope.y), col = "red") +
  facet_wrap(~lineage2, scales = "free_y", ncol = 6) + theme_classic() +
  xlab(expression("log" ~ beta[IBD])) + 
  theme(strip.text = element_text(size = 7, face= "italic"))
save_plot("~/Desktop/IBD_bootstrapping.pdf", ab, base_height = 5,
          base_width = 10)  

dd3 = dd2 %>% group_by(lineage) %>% 
  summarise(num = mean(n.x), full_log_slope = mean(log_slope.y),
            se = mean(se.y), r2 = mean(r2.y), 
            sample_sd = sd(log_slope.x), sample_se = mean(se.x),
            sample_r2 = mean(r2.x), sample_log_slope = mean(log_slope.x),) %>% ungroup()

a = ggplot(dd3, aes(full_log_slope, sample_log_slope)) + 
  geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
  geom_point(shape = 21, fill = "dodgerblue") + 
  xlab(expression("log" ~ beta[IBD] * ", all inds")) +
  ylab(expression("mean subsampled log" ~ beta[IBD]))
b = ggplot(dd3, aes(r2, sample_r2)) + 
  geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
  geom_point(shape = 21, fill = "dodgerblue") + 
  xlab(expression(r^2 ~ "all inds")) +
  ylab(expression("mean subsampled" ~ r^2)) +
  ylim(0, 0.7) + xlim(0, 0.7)
c = ggplot(dd3, aes(r2, sample_sd)) + 
  geom_point(shape = 21, fill = "dodgerblue") + 
  xlab(expression(r^2 ~ "all inds")) +
  ylab(expression("mean subsampled SD")) 
abc = plot_grid(a, b, c, labels = c("A", "B", "C"), ncol = 3)
save_plot("~/Desktop/subsampling.pdf", abc, base_height = 4,
          base_width = 12)
