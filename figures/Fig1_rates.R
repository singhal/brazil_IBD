library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)

setwd("~/Dropbox (Personal)/brazil_gene_flow/cooney-data/")

d = list.files(".", pattern = "csv", full.names = T)
d1 = lapply(d, read.csv)
d2 = do.call("rbind", d1)
xlim = range(d2$mean.dr.full.tree)

sp = gsub(".*rates.", "", d)
sp = gsub(".csv", "", sp)

vals = rep(NA, length(d1))
for (i in 1:length(d1)) {
  vals[i] = quantile(d1[[i]]$mean.dr.full.tree, probs = c(0.975)) / 
    quantile(d1[[i]]$mean.dr.full.tree, probs = c(0.025)) 
}



d3 = vector("list", length = length(d))
plts = vector("list", length = length(d))
cols = brewer.pal(length(d), "Paired")
for (i in 1:length(d)) {
  x = read.csv(d[i], stringsAsFactors = F)
  x$clade = sp[i]
  plts[[i]] = ggplot(x, aes(mean.dr.full.tree)) + 
    geom_histogram(bins = 200, fill = cols[i]) + 
    theme_classic() +
    scale_x_log10(limits = xlim) +
    xlab("speciation rate")
  d3[[i]] = x
}
d4 = do.call("rbind", d3)
cols = brewer.pal(6, "Greys")
a = ggplot(d4, aes(mean.dr.full.tree)) + 
  geom_histogram(aes(fill = clade), bins = 150) +
  scale_x_log10() + xlab("speciation rate") + 
  scale_fill_manual(values = cols[2:6]) + 
  theme_classic()
save_plot("~/Desktop/speciation_rates_verts2.pdf", a,
          base_width = 6, base_height = 4)

pltA = ggplot(d4, aes(bamm.mcc.tree)) + 
  geom_histogram(aes(fill = clade), bins = 50) +
  scale_x_log10() + xlab(expression("speciation rate (" * lambda[BAMM] * ")")) + 
  scale_fill_manual(values = cols[2:6]) + 
  theme_classic()
pltB = ggplot(d4, aes(dr.mcc.tree, bamm.mcc.tree)) + 
  geom_point(alpha = 0.1, shape = 16) + theme_classic() + 
  scale_x_log10() + scale_y_log10() +
  xlab(expression("speciation rate (" * lambda[DR] * ")")) +
  ylab(expression("speciation rate (" * lambda[BAMM] * ")"))
ab = plot_grid( pltA, pltB, ncol = 2, labels = c("A", "B"))
save_plot("~/Desktop/vertebrate_speciation_BAMM_DR.png", 
          ab, base_height = 3.5, base_width = 8)

cort = cor.test(log(d4$dr.mcc.tree), log(d4$bamm.mcc.tree))
vals2 = rep(NA, length(d1))
for (i in 1:length(d1)) {
  vals2[i] = quantile(d1[[i]]$bamm.mcc.tree, probs = c(0.99)) / 
    quantile(d1[[i]]$bamm.mcc.tree, probs = c(0.01)) 
}


layout <- "
ABC
DE#
"
aa = plts[[1]] + plts[[2]] + plts[[3]] + plts[[4]] + plts[[5]]  + plot_layout(design = layout)
save_plot("~/Desktop/speciation_rates_verts.pdf", aa,
          base_width = 9, base_height = 5)
