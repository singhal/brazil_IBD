library(mvtnorm)
library(ape)
library(phytools)
library(nlme)

source("~/Dropbox (Personal)/scripts/brazil/essim.R")

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
rcorrs = seq(0, 1, by=0.1)
nrep = c(100)
res = vector('list', length(rcorrs))
for (i in 1:length(res)) {
  res[[i]] = list(rep(NA, nrep), 
                  rep(NA, nrep), 
                  rep(NA, nrep),
                  rep(NA, nrep),
                  rep(NA, nrep))
  names(res[[i]]) = c("lambda", "sig_gls", "cor_gls",
                      "sig_es", "cor_es")
}
names(res) = rcorrs
for (j in 1:length(rcorrs)) {
  cat(j, "\n")
  for (i in 1:nrep) {
    # Generate traits under Brownian motion:
    d1 <- as.vector(rmvnorm(1, sigma=vv))

    # rotate them by the target correlation with the speciation rates:
    newvals <- (rcorrs[j] * xvals) + sqrt(1 - rcorrs[j]^2)*d1
    # The newvals will have the approximate correlation desired 
    #   and will also have the right correlation due to phylogeny
    
    # test lambda
    res[[j]]$lambda[i] = phylosig(t2, newvals[t2$tip.label], 
                                  method="lambda")$lambda
  
    # test significance
    d = data.frame(xvals[t2$tip.label], newvals[t2$tip.label])
    names(d) = c("xvals", "newvals")
    g1 = essim(t2, newvals[t2$tip.label], nsim = 1000,
               xvals[t2$tip.label], logrates = F)
    g2 = gls(xvals ~ newvals, 
            correlation=corBrownian(phy=t2), 
            data=d, method="ML")
    
    obj = phyl.vcv(as.matrix(cbind(newvals, xvals)), vcv(t2), 1)
    
    res[[j]]$sig_gls[i] = summary(g2)$tTable[2, 4]
    res[[j]]$cor_gls[i] = cov2cor(obj$R)[1, 2]
    res[[j]]$sig_es[i] = g1[2]
    res[[j]]$cor_es[i] = g1[1]
  }
}

xx = lapply(res, function(x) {do.call("cbind", x)})
xx2 = data.frame(do.call("rbind", xx))
xx2$cor = rep(seq(0, 1, by=0.1), each = nrep)

power_es = unlist(lapply(res, function(x) {length(x$sig_es[x$sig_es < 0.05])})) / 100
power_gls = unlist(lapply(res, function(x) {length(x$sig_gls[x$sig_gls < 0.05])})) / 100
res1 = data.frame(cor = rcorrs, power_es, power_gls)

powerres = list(xx2, res1)
saveRDS(powerres, "~/Dropbox (Personal)/brazil_gene_flow/results/power.Rds")

d = readRDS("~/Dropbox (Personal)/brazil_gene_flow/results/power.Rds")
