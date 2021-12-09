library(ape)

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
d3 = d2[match(t2$tip.label, d2$lineage), ]

r2 = c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
res = list("vector", length(r2))
for (i in 1:length(r2)) {
  dd = d3[d3$r2 > r2[i], ]
  t3 = keep.tip(t, dd$lineage)
  dd1 = dd[match(t3$tip.label, dd$lineage), ]
  
  phyd = phylosig(t3, dd1$log_slope, test = TRUE,
                     method = "lambda", nsim = 1000)
  
  # Test beta_IBD vs ClaDS speciation rate correlation via PGLS
  glsd = gls(log_crates ~ log_slope, 
             correlation = corBrownian(phy = t3),
             data = dd1, method = "ML")
  
  obj = phyl.vcv(as.matrix(dd1[, c("log_slope", "log_crates")]), vcv(t3), 1)
  corvar = cov2cor(obj$R)[1, 2]
  t.xy = corvar * sqrt((Ntip(t3)-2)/(1-corvar^2))
  pval = 2*pt(abs(t.xy),df=Ntip(t3)-2,lower.tail=F)
  
  # r2, sample size, correlation, pval
  n = summary(glsd)$dims$N
  pval = round(summary(glsd)$tTable[2, 4], 3)
  res[[i]] = c(r2[i], n, corvar, pval)
}

res2 = data.frame(do.call("rbind", res))
names(res2) = c("r2", "n", "cor", "pval")
seplot = ggplot(res2, aes(r2, cor, size = n, fill = pval)) +
  geom_point(shape = 21) +
  xlab(expression("minimum" ~ r^2 ~ "of IBD slope")) +
  ylab("correlation")  + 
  scale_fill_continuous(type = "viridis")

lins = table(ll$lineage)
d3$num = lins[d3$lineage]

ss = c(4, 5, 6, 7, 8, 9, 10)
res = list("vector", length(ss) - 1)
for (i in 1:(length(ss) - 1)) {
  dd = d3[d3$num >= ss[i], ]
  t3 = keep.tip(t, dd$lineage)
  dd1 = dd[match(t3$tip.label, dd$lineage), ]
  
  phyd = phylosig(t3, dd1$log_slope, test = TRUE,
                  method = "lambda", nsim = 1000)
  
  # Test beta_IBD vs ClaDS speciation rate correlation via PGLS
  glsd = gls(log_crates ~ log_slope, 
             correlation = corBrownian(phy = t3),
             data = dd1, method = "ML")
  
  obj = phyl.vcv(as.matrix(dd1[, c("log_slope", "log_crates")]), vcv(t3), 1)
  corvar = cov2cor(obj$R)[1, 2]
  t.xy = corvar * sqrt((Ntip(t3)-2)/(1-corvar^2))
  pval = 2*pt(abs(t.xy),df=Ntip(t3)-2,lower.tail=F)
  
  # r2, sample size, correlation, pval
  n = summary(glsd)$dims$N
  pval = round(summary(glsd)$tTable[2, 4], 3)
  res[[i]] = c(ss[i], n, corvar, pval)
}

res2 = data.frame(do.call("rbind", res))
names(res2) = c("ss", "n", "cor", "pval")
ssplot = ggplot(res2, aes(ss, cor, size = n, fill = pval)) +
  geom_point(shape = 21) +
  xlab(expression("minimum sample size")) +
  ylab("correlation")  + 
  scale_fill_continuous(type = "viridis")

prow <- plot_grid(
  ssplot + theme(legend.position="none"),
  seplot + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B"),
  hjust = -1,
  nrow = 1
)
legend <- get_legend(
  # create some space to the left of the legend
  ssplot + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
ab = plot_grid(prow, legend, rel_widths = c(3, .4))
save_plot("~/Desktop/r2_ss_plot.png", ab, base_height = 3, base_width = 8)
