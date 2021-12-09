inv_fst <- function(x) {
  inv = x / (1 - x)
  return(inv)
}

# setwd("~/Dropbox (Personal)/brazil_gene_flow/data/IBD_simulations/summary/")
setwd("~/Desktop/sims/summary/")
dd = list.files(pattern = "*csv")

res = vector("list", length(dd))

for (i in 1:length(dd)) {
  d = read.csv(dd[i])
  d$log_dist = log(d$distance)
  d$inv_fst = inv_fst(d$fst)
  
  
  sigma = d$sigma[1]
  rep = d$rep[1]
  
  inds = seq(0, 100)
  ninds = c(5, 10, 15, 20)
  reps = 10
  dslopes = c()
  fslopes = c()
  indstrack = c()
  for (x in 1:length(ninds)) {
    nind = ninds[x]
    currep = reps[x]
    
    for (y in 1:reps) {
      curinds = sample(inds, nind)
      
      d1 = d[d$ind1 %in% curinds, ]
      d1 = d1[d1$ind2 %in% curinds, ]
      
      f1 = f[f$ind1 %in% curinds, ]
      f1 = f1[f1$ind2 %in% curinds, ]
      
      m1 = lm(nuc_dxy ~ log_dist, data = d1)
      
      dslopes = c(dslopes, coef(m1)[2])
      indstrack = c(indstrack, nind)
    }
  }
  
  res[[i]] = data.frame(sigma, rep, dslopes, indstrack)
}

res2 = do.call("rbind", res)
library(ggplot2)
library(RColorBrewer)
ggplot(res2, aes(as.factor(indstrack), log(dslopes))) + 
  geom_boxplot(aes(fill = as.factor(sigma)), outlier.shape = NA) + 
  facet_wrap(~sigma, ncol = 4) +
  theme_classic() + xlab("# of sampled inds") +
  ylab(expression("log" ~ beta[IBD])) +
  scale_fill_manual(values = brewer.pal(4, "BuPu")) +
  theme(legend.position = "none")
save_plot("~/Desktop/ibd_sims.png", xx, base_width = 6.5,
          base_height = 2.5)

ggplot(res2, aes(as.factor(indstrack), dslopes)) + 
  geom_boxplot() + facet_wrap(~sigma, ncol = 4) +
  theme_classic()
  