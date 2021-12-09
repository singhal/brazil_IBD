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

inv_fst <- function(x) {
  inv = x / (1 - x)
  return(inv)
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
x1 = read.csv("~/Dropbox (Personal)/brazil_gene_flow/data/brazil_samples_v8_nonCerrado.csv",
              stringsAsFactors = F)
x1[x1$lineage == "Vanzosaura_savanicola", "lineage"] = "Vanzosaura_rubricauda_2"
x2 = x1[x1$keep == TRUE,]
div = list.files("~/Dropbox (Personal)/brazil_gene_flow/data/pop_gen/",
                 pattern = "csv", full.names = T)
div = div[grep("teguixin", div, invert = T)]


# what to output
# 1. all divergences across species
# 2. all slopes: # comp, val, sig, r2
# 3. figure of all three
res1 = vector("list", length(dist))
res2 = vector("list", length(dist))


for (i in 1:length(div)) {
  dist = read.csv(div[i], stringsAsFactors = F)
  
  dist = dist[dist$ind1 %in% x2$sample, ]
  dist = dist[dist$ind2 %in% x2$sample, ]
  
  if (nrow(dist) > 2) {
    dist$geo_dist = apply(dist, 1, get_dist) / 1000
    dist$ln_geo_dist = log(dist$geo_dist)
    dist$inv_fst = inv_fst(dist$fst)
   
   
    res1[[i]] = dist 
    lineage = dist$lineage[1]
    
    a = get_slope('inv_fst', 'ln_geo_dist')
    dd = c(a[[1]], lineage)
  } else {
    dd = c(rep(NA, 10), lineage)
  }
  res2[[i]] = dd
}

res1a = do.call("rbind", res1)
res2a = as.data.frame(do.call("rbind", res2))
names(res2a) = c("div_type", "dist_type",
                 "slope", "intercept", "se", "n", "r2", "sig",
                 "r2_mt", "sig_mt",
                 "lineage")
for (i in c("slope", "n", "r2", "sig", "r2_mt", "sig_mt")) {
  res2a[, i] = as.numeric(res2a[, i])
}

d = res2a[complete.cases(res2a$slope), ]
d1 = d[d$n > 6, ]
d1$log_slope = log(d1$slope)
d2 = d1[complete.cases(d1$log_slope), ]
d2[which(d2$lineage == "Vanzosaura_rubricauda_2"), "lineage"] = "Vanzosaura_savanicola"

e = read.csv("~/Dropbox (Personal)/brazil_gene_flow/results/IBD-2021-09-15.csv")
e$log_slope = log(e$slope)
e1 = e %>% filter(div_type == "inv_fst")
e2 = e1[complete.cases(e1$log_slope), ]

f = left_join(d2, e2, by = c("lineage" = "lineage"))
cor.test(f$slope.x, f$slope.y)
ggplot(f, aes(log_slope.y, log_slope.x)) + geom_point() +
  xlab(expression("log" ~ beta[IBD] ~ "all samples")) + 
  ylab(expression("log" ~ beta[IBD] ~ "Cerrado samples"))


s = readRDS("~/Dropbox (Personal)/brazil_gene_flow/data/sample_map.Rds")
crates = readRDS("~/Dropbox (Personal)/brazil_gene_flow/data/speciation_rates/clads_tiprates.Rdata")

t = read.tree("~/Dropbox (Personal)/brazil_gene_flow/data/phylogeny/Tonini.squam_shl_new_Consensus_9755.targeted_species.tre")
t$tip.label[which(t$tip.label == "Vanzosaura_rubricauda_2")] = "Vanzosaura_savanicola"

# Calculate phylogenetic signal of beta_IBD
t2 = keep.tip(t, d2$lineage)
d3 = d2[match(t2$tip.label, d2$lineage), ]
physig1 = phylosig(t2, d3$log_slope, test = TRUE,
                   method = "lambda", nsim = 1000)

# Test beta_IBD vs ClaDS speciation rate correlation via PGLS
crates2 = crates[ s[match(d3$lineage, s$lineage), "full_tree"] ]
d3$log_crates = log(crates2)
rownames(d3) = d3$lineage
d4 = d3[t2$tip.label, ]
gls1 = gls(log_crates ~ log_slope, 
           correlation = corBrownian(phy = t2),
           data = d4, method = "ML")

obj = phyl.vcv(as.matrix(d4[, c("log_slope", "log_crates")]), vcv(t2), 1)
corvar = cov2cor(obj$R)[1, 2]
t.xy = corvar * sqrt((Ntip(t2)-2)/(1-corvar^2))
pval = 2*pt(abs(t.xy),df=Ntip(t2)-2,lower.tail=F)

