library(mvtnorm)
library(ape)
library(phytools)
library(nlme)

t <- read.tree("~/Dropbox (Personal)/brazil_gene_flow/data/phylogeny/Tonini.squam_shl_new_Consensus_9755.targeted_species.tre")
diff = read.csv("~/Dropbox (Personal)/brazil_gene_flow/results/IBD-2020-12-02.csv",  stringsAsFactors = F)
crates = readRDS("~/Dropbox (Personal)/brazil_gene_flow/data/speciation_rates/clads_tiprates.Rdata")
ldf = readRDS("~/Dropbox (Personal)/brazil_gene_flow/data/sample_map.Rds")

t$tip.label[which(t$tip.label == "Vanzosaura_rubricauda_2")] = "Vanzosaura_savanicola"

diff = diff[diff$dist_type == 'ln_geo_dist' & diff$div_type == 'inv_fst',]
diff$log_slope = log(diff$slope)
diff = diff[complete.cases(diff$log_slope), ]
ibd = diff$log_slope
names(ibd) = diff$lineage
names(ibd)[which(names(ibd) == "Vanzosaura_rubricauda_2")] = "Vanzosaura_savanicola"

t2 = keep.tip(t, names(ibd))

vv <- vcv.phylo(t2, corr = TRUE)
# let xvals be log-transformed measures of speciation rate
#    in same order as tip labels, so
#    if mytraits is vector of rates with names

rates <- crates[ ldf[ match(t2$tip.label, ldf$lineage), "full_tree"] ]
xvals <- log(rates)
names(xvals) = t2$tip.label

# Define a target correlation 
rcorrs = runif(100000, min = -1, max = 1)
res = data.frame(actualrho = rcorrs,
                 pval = rep(NA, length(rcorrs)),
                 estrho = rep(NA, length(rcorrs)))
for (i in 1:length(rcorrs)) {
  # Generate traits under Brownian motion:
  d1 <- as.vector(rmvnorm(1, sigma=vv))
    
  # rotate them by the target correlation with the speciation rates:
  newvals <- (rcorrs[i] * xvals) + sqrt(1 - rcorrs[i]^2)*d1
  # The newvals will have the approximate correlation desired 
  #   and will also have the right correlation due to phylogeny
    
  # test significance
  d = data.frame(xvals[t2$tip.label], newvals[t2$tip.label])
  names(d) = c("xvals", "newvals")
  g2 = gls(xvals ~ newvals, correlation=corBrownian(phy=t2), 
             data=d, method="ML")
  
  obj = phyl.vcv(as.matrix(cbind(newvals, xvals)), vcv(t2), 1)
  
  res[i, "pval"] = summary(g2)$tTable[2, 4]
  res[i, "estrho"] = cov2cor(obj$R)[1, 2]
  if (i %% 1000 == 0) {
    cat(i, "\n")
  }
}

write.csv(res, "~/Desktop/power.csv", row.names = F)

res = read.csv("~/Dropbox (Personal)/brazil_gene_flow/results/power.csv", stringsAsFactors = F)
res2 = res[res$estrho > 0.02 & res$estrho < 0.06, ]
hist(res2$actualrho, breaks = 30)
