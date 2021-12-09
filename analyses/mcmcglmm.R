require(MCMCglmm)

# check_and_fix_ultrametric
#     checks trees to see if they pass ape ultrametricity test.
# If not, it computes the differential root-to-tip distance across all tips.
# It adds the appropriate quantity to each terminal branch length to ensure that 
# tree passes ultrametric test.
# Note: this is only a valid method of making trees ultrametric when the 
# 	non-ultrametricity is due to small numerical discrepancies, e.g., 
#   rounding or other floating point issues during phylogeny construction.
# 

check_and_fix_ultrametric <- function(phy){
  
  if (!is.ultrametric(phy)){
    
    vv <- vcv.phylo(phy)
    dx <- diag(vv)
    mxx <- max(dx) - dx
    for (i in 1:length(mxx)){
      phy$edge.length[phy$edge[,2] == i] <- phy$edge.length[phy$edge[,2] == i] + mxx[i]
    }
    if (!is.ultrametric(phy)){
      stop("Ultrametric fix failed\n")
    }	
  }
  
  return(phy)
}

# from https://ourcodingclub.github.io/tutorials/mcmcglmm/
trace.plots <- function(x) {
  n <- dim(x)[2]
  par(mfrow=c(ceiling(n/2),2), mar=c(0,0.5,1,0.5))
  for (i in 1:n) {
    plot(as.numeric(x[,i]), t="l", main=colnames(x)[i], xaxt="n", yaxt="n")
  }
}

s = readRDS("~/Dropbox (Personal)/brazil_gene_flow/data/sample_map.Rds")
ll = read.csv("~/Dropbox (Personal)/brazil_gene_flow/data/brazil_samples_v8.csv",
              stringsAsFactors = F)
crates = readRDS("~/Dropbox (Personal)/brazil_gene_flow/data/speciation_rates/clads_tiprates.Rdata")

t = read.tree("~/Dropbox (Personal)/brazil_gene_flow/data/phylogeny/Tonini.squam_shl_new_Consensus_9755.targeted_species.tre")
t$tip.label[which(t$tip.label == "Vanzosaura_rubricauda_2")] = "Vanzosaura_savanicola"

# Calculate phylogenetic signal of beta_IBD
d = read.csv("~/Dropbox (Personal)/brazil_gene_flow/results/IBD-2021-09-15.csv", stringsAsFactors = F)
d1 = d %>% dplyr::filter(div_type == 'inv_fst')
d1$log_slope = log(d1$slope)
d2 = d1[complete.cases(d1$log_slope), ]
d2[which(d2$lineage == "Vanzosaura_rubricauda_2"), "lineage"] = "Vanzosaura_savanicola"
crates2 = crates[ s[match(d2$lineage, s$lineage), "full_tree"] ]
d2$log_crates = log(crates2)

t2 = keep.tip(t, d2$lineage)
t3 = check_and_fix_ultrametric(t2)
d3 = d2[match(t3$tip.label, d2$lineage), ]

inv.phylo = inverseA(t3, nodes="TIPS", scale=TRUE)
prior1 = list(G=list(G1=list(V=1,nu=0.02)),
             R=list(V=1,nu=0.02))
x1 = MCMCglmm(log_crates ~ log_slope,
         data = d3,
         random = ~lineage,
         family = "gaussian",
         ginverse=list(lineage=inv.phylo$Ainv),
         prior=prior1,
         thin = 1,
         burnin = 300,
         nitt = 4000)

prior2 = list(G=list(G1=list(V=1,nu=0.02),
                     G1=list(V=1,nu=0.02)),
              R=list(V=1,nu=0.02))
x2 = MCMCglmm(log_crates ~ log_slope,
              data = d3,
              random = ~lineage + idh(se):units,
              family = "gaussian",
              ginverse=list(lineage=inv.phylo$Ainv),
              prior=prior2,
              thin = 10,
              burnin = 3000,
              nitt = 30000)
trace.plots(x2$Sol)
plot(x2$VCV)
summary(x2)

lambda <- x2$VCV[,'lineage']/
  (x2$VCV[,'lineage']+x2$VCV[,'units'])
mean(lambda)
posterior.mode(lambda)
HPDinterval(lambda)
