library(nlme)
library(ape)
library(dplyr)
library(rgdal)
library(phytools)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

s = readRDS("~/Dropbox (Personal)/brazil_gene_flow/data/sample_map.Rds")
ll = read.csv("~/Dropbox (Personal)/brazil_gene_flow/data/brazil_samples_v8.csv",
              stringsAsFactors = F)
crates = readRDS("~/Dropbox (Personal)/brazil_gene_flow/data/speciation_rates/clads_tiprates.Rdata")

d = read.csv("~/Dropbox (Personal)/brazil_gene_flow/data/speciation_rates/tonini-dr-rates-100.csv", stringsAsFactors = F, header = F)
drates = apply(d[2:101], 1, mean)
names(drates) = d$V1

t = read.tree("~/Dropbox (Personal)/brazil_gene_flow/data/phylogeny/Tonini.squam_shl_new_Consensus_9755.targeted_species.tre")
t$tip.label[which(t$tip.label == "Vanzosaura_rubricauda_2")] = "Vanzosaura_savanicola"

# Calculate phylogenetic signal of beta_IBD
d = read.csv("~/Dropbox (Personal)/brazil_gene_flow/results/IBD-2021-09-14.csv", stringsAsFactors = F)
d1 = d %>% dplyr::filter(div_type == 'inv_fst')
d1$log_slope = log(d1$slope)
d2 = d1[complete.cases(d1$log_slope), ]
d2[which(d2$lineage == "Vanzosaura_rubricauda_2"), "lineage"] = "Vanzosaura_savanicola"

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

# # Test beta_IBD vs range size correlation via PGLS
# r = readOGR("~/Dropbox/inornatus_gr/Meiri_etal_ranges/GARD1.1_dissolved_ranges/modeled_reptiles.shp")
# rnames = r[[1]]
# rnames = gsub(" ", "_", rnames)
# 
# spnames = d4$lineage
# names(spnames) = spnames
# spnames[! spnames %in% rnames]
# 
# spnames["Vanzosaura_rubricauda_2"] = "Vanzosaura_rubricauda" 
# spnames["Bachia_bresslaui_1"] = "Bachia_bresslaui" 
# spnames["Trilepida_brasiliensis_2"] = "Trilepida_brasiliensis" 
# spnames["Xenodon_merremii"] = "Xenodon_merremi" 
# spnames["Taeniophallus_occipitalis_2"] = "Taeniophallus_occipitalis"
# spnames["Oxyrhopus_trigeminus_2"] = "Oxyrhopus_trigeminus" 
# spnames["Spilotes_sp"] = "Spilotes_sulphureus" 
# 
# ranges = rep(NA, length(spnames))
# for (i in 1:length(ranges)) {
#   rnum = which(rnames == spnames[i])
#   range = r[rnum, 1]
#   ranges[i] = rgeos::gArea(range)
# }
# d4$area = ranges
# d4$log_area = log(d4$area)
# gls_area = gls(log_area ~ log_slope, 
#            correlation = corBrownian(phy = t2),
#            data = d4, method = "ML")
# areaplt = ggplot(d4) + geom_point(aes(log_slope, log_area),
#                         shape = 21, 
#                         fill = "dodgerblue") +
#   xlab("log IBD slope") + ylab("log range area")
# save_plot("~/Dropbox/brazil_gene_flow/figures/area_IBD.pdf",
#           areaplt)

# Test beta_IBD vs ClaDS speciation rate correlation via ES-SIM
source("~/Dropbox (Personal)/scripts/brazil/essim.R")
trait = d4$log_slope
names(trait) = rownames(d4)
es = d4$log_crates
names(es) = rownames(d4)
essim(t2, trait, nsim = 1000, es, return.es = FALSE, logrates = FALSE)

# Test beta_IBD vs DR speciation rate correlation via PGLS
drates2 = drates[ s[match(d3$lineage, s$lineage), "full_tree"] ]
d3$log_drates = log(drates2)
rownames(d3) = d3$lineage
d4 = d3[t2$tip.label, ]
gls2 = gls(log_drates ~ log_slope, 
           correlation = corBrownian(phy = t2),
           data = d4, method = "ML")

# Test beta_IBD vs DR speciation rate correlation via ES-SIM
es2 = d4$log_drates
names(es2) = rownames(d4)
essim(t2, trait, nsim = 1000, es2, return.es = FALSE, logrates = FALSE)

# Remove species with only 4 sampled individuals
# 4 sampled inds = 6 comparisons
d5 = d4[d4$n > 6, ]
t3 = keep.tip(t2, d5$lineage)
d5 = d5[t3$tip.label, ]
gls3 = gls(log_crates ~ log_slope, 
           correlation = corBrownian(phy = t3),
           data = d5, method = "ML")
physig3 = phylosig(t3, d5$log_slope, test = TRUE,
         method = "lambda", nsim = 1000)

# Remove species with non-significant beta_IBD slopes
d5 = d4[d4$sig <= 0.05, ]
t3 = keep.tip(t2, d5$lineage)
d5 = d5[t3$tip.label, ]
gls4 = gls(log_crates ~ log_slope, 
           correlation = corBrownian(phy = t3),
           data = d5, method = "ML")
physig4 = phylosig(t3, d5$log_slope, test = TRUE,
                   method = "lambda", nsim = 1000)

obj = phyl.vcv(as.matrix(d5[, c("log_slope", "log_crates")]), vcv(t3), 1)
corvar = cov2cor(obj$R)[1, 2]
t.xy = corvar * sqrt((Ntip(t3)-2)/(1-corvar^2))
pval = 2*pt(abs(t.xy),df=Ntip(t3)-2,lower.tail=F)

# Use nuclear d_xy and geographic distance
ndy = d[d$div_type == "nuc_dxy", ] 
ndy$log_slope = log(ndy$slope)
ndy2 = ndy[complete.cases(ndy$log_slope), ]
rownames(ndy2) = ndy2$lineage
td = keep.tip(t, ndy2$lineage)
ndy3 = ndy2[td$tip.label, ]
physig5 = phylosig(td, ndy3$log_slope, test = TRUE,
         method = "lambda", nsim = 1000)
crates3 = crates[ s[match(ndy3$lineage, s$lineage), "full_tree"] ]
ndy3$log_crates = log(crates3)
gls5 = gls(log_crates ~ log_slope, 
           correlation = corBrownian(phy = td),
           data = ndy3, method = "ML")

# Use mtDNA d_xy and geographic distance
ndy = d[d$div_type == "mt_dist", ] 
ndy$log_slope = log(ndy$slope)
ndy2 = ndy[complete.cases(ndy$log_slope), ]
rownames(ndy2) = ndy2$lineage
td = keep.tip(t, ndy2$lineage)
ndy3 = ndy2[td$tip.label, ]
physig6 = phylosig(td, ndy3$log_slope, test = TRUE,
         method = "lambda", nsim = 1000)
crates3 = crates[ s[match(ndy3$lineage, s$lineage), "full_tree"] ]
ndy3$log_crates = log(crates3)
gls6 = gls(log_crates ~ log_slope, 
           correlation = corBrownian(phy = td),
           data = ndy3, method = "ML")

# Use F_ST and environmental distance
ndy = d[d$div_type == "fst", ] 
ndy$log_slope = log(ndy$slope)
ndy2 = ndy[complete.cases(ndy$log_slope), ]
rownames(ndy2) = ndy2$lineage
td = keep.tip(t, ndy2$lineage)
ndy3 = ndy2[td$tip.label, ]
physig7 = phylosig(td, ndy3$log_slope, test = TRUE,
         method = "lambda", nsim = 1000)
crates3 = crates[ s[match(ndy3$lineage, s$lineage), "full_tree"] ]
ndy3$log_crates = log(crates3)
gls7 = gls(log_crates ~ log_slope, 
           correlation = corBrownian(phy = td),
           data = ndy3, method = "ML")

# Nominal species
dd = read.csv("~/Dropbox (Personal)/brazil_gene_flow/results/nominal_species_IBD-6-Apr-2021.csv",
              stringsAsFactors = F)
t = read.tree("~/Dropbox (Personal)/brazil_gene_flow/data/speciation_rates/Squamata_Tree_Methods_Data/squam_shl_new_Consensus_9755.tre")
dd$species[ ! dd$species %in% t$tip.label ]
dd$tree_species = dd$species

dd[dd$tree_species == "Cercosaura_argula", "tree_species"] = "Cercosaura_argulus"       
dd[dd$tree_species == "Spilotes_sp", "tree_species"] = "Spilotes_sulphureus"     
dd[dd$tree_species == "Xenodon_merremii", "tree_species"] = "Xenodon_merremi"      
dd$log_slope = log(dd$slope)
dd2 = dd[complete.cases(dd$log_slope), ]
t2 = keep.tip(t, dd2$tree_species)
dd3 = dd2[match(t2$tip.label, dd2$tree_species), ]
physig8 = phylosig(t2, dd3$log_slope, test = TRUE,
         method = "lambda", nsim = 1000)
dd3$log_crates = log(crates[ dd3$tree_species ])
gls8 = gls(log_crates ~ log_slope, 
           correlation = corBrownian(phy = t2),
           data = dd3, method = "ML")


# calculate IBD accounting for intercept at 1000 km
d1$fst_1000 = d1$slope * log(1000) + d1$intercept
d1$log_fst_1000 = log(d1$fst_1000)
d2 = d1[complete.cases(d1$log_fst_1000), ]
d2[which(d2$lineage == "Vanzosaura_rubricauda_2"), "lineage"] = "Vanzosaura_savanicola"
d2 = d2[d2$lineage != "Phrynonax_poecilonotus", ]
t3 = keep.tip(t, d2$lineage)
d3 = d2[match(t3$tip.label, d2$lineage), ]
d3$log_crates = log(crates[ s[match(d3$lineage, s$lineage), "full_tree"] ])

gls10 = gls(log_crates ~ log_fst_1000, 
           correlation = corBrownian(phy = t3),
           data = d3, method = "ML")
physig10 = phylosig(t3, d3$log_fst_1000, test = TRUE,
                   method = "lambda", nsim = 1000)

summarize_results <- function(name, physig, pgls) {
  lambda = round(physig$lambda, 2)
  pval1 = round(physig$P, 3)
  
  n = summary(pgls)$dims$N
  pval2 = round(summary(pgls)$tTable[2, 4], 3)
  corr = summary(pgls)$tTable[2,1]
  
  corrval = "+"
  if (corr < 0) {
    corrval = '-'
  }
  
  return(c(name, n, lambda, pval1, corrval, pval2))
}

res = list(
summarize_results("main", physig1, gls1),
summarize_results(">4inds", physig3, gls3),
summarize_results("sig_IBD", physig4, gls4),
summarize_results("nuc_dxy", physig5, gls5),
summarize_results("mt_dxy", physig6, gls6),
summarize_results("env_dist", physig7, gls7),
summarize_results("nominal", physig8, gls8))
summarize_results("ibd at 1000km", physig10, gls10)
res2 = as.data.frame(do.call("rbind", res))
write.csv(res2, "~/Dropbox (Personal)/brazil_gene_flow/results/robustness.csv",
          quote = F, row.names = F)
