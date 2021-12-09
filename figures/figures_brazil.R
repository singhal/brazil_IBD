library(ape)
library(cowplot)
library(classInt)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(rgdal)
library(tidyr)
library(jpeg)
library(phytools)
library(grid)
library(phangorn)
library(BAMMtools)
library(plotfunctions)
library(patchwork)
library(raster)
library(viridis)
theme_set(theme_cowplot())

source("~/Dropbox (Personal)/scripts/brazil//colorMap.R")

maincol2 = '#e41a1c'
maincol = '#377eb8'

###################
# keep lineages
###################

ll = read.csv("~/Dropbox (Personal)/brazil_gene_flow/data/brazil_samples_v8.csv",
              stringsAsFactors = F)
ll2 =  ll[complete.cases(ll$LAT), ]
lins = table(ll2$lineage)
lins = names(lins[ lins > 3])
lins = lins[lins != "Phrynonax_poecilonotus"]
inds = ll2[ll2$lineage %in% lins, "sample"]
ll3 =  ll2[ll2$lineage %in% lins, ]

# update lins to match names in trees
ldf = data.frame(lineage = lins,
                 full_tree = lins,
                 gene_tree = lins,
                 stringsAsFactors = F)
ldf[which(ldf$lineage == "Bachia_bresslaui_1"), 
    c("full_tree", "gene_tree")] = 
  "Bachia_bresslaui"
ldf[which(ldf$lineage == "Copeoglossum_nigropunctatum_1"), 
    c("full_tree", "gene_tree")] = 
  "Copeoglossum_nigropunctatum"
ldf[which(ldf$lineage == "Spilotes_sp"), 
    c("full_tree", "gene_tree")] = 
  "Spilotes_sulphureus"
ldf[which(ldf$lineage == "Taeniophallus_occipitalis_2"), 
    c("full_tree", "gene_tree")] = 
  "Taeniophallus_occipitalis"
ldf[which(ldf$lineage == "Trilepida_brasiliensis_2"), 
    c("full_tree", "gene_tree")] = 
  "Trilepida_brasiliensis"
ldf[which(ldf$lineage == "Vanzosaura_savanicola"), 
    c("full_tree", "gene_tree")] = 
  "Vanzosaura_rubricauda"
ldf[which(ldf$lineage == "Xenodon_merremii"), 
    c("full_tree", "gene_tree")] = 
  "Xenodon_merremi"
ldf[which(ldf$lineage == "Oxyrhopus_trigeminus_1"), 
    c("full_tree", "gene_tree")] = 
  "Oxyrhopus_trigeminus"
ldf[which(ldf$lineage == "Oxyrhopus_trigeminus_2"), 
    c("full_tree", "gene_tree")] = 
  "Oxyrhopus_trigeminus"

# based on our tree
ldf[which(ldf$lineage == "Apostolepis_ammodites"), 
    c("gene_tree")] = "Apostolepis_cearensis"
ldf[which(ldf$lineage == "Helicops_modestus"), 
    c("gene_tree")] = "Helicops_carinicaudus"
ldf[which(ldf$lineage == "Taeniophallus_occipitalis_2"), 
    c("gene_tree")] = "Taeniophallus_nicagus"
ldf[which(ldf$lineage == "Tantilla_melanocephala"), 
    c("gene_tree")] = "Tantilla_armillata"
ldf[which(ldf$lineage == "Trilepida_brasiliensis_2"), 
    c("gene_tree")] = "Trilepida_macrolepis"
ldf[which(ldf$lineage == "Trilepida_koppesi"), 
    c("gene_tree")] = "Trilepida_macrolepis"
# saveRDS(ldf, "~/Dropbox (Personal)/brazil_gene_flow/data/sample_map.Rds")

ll = read.csv("~/Dropbox (Personal)//brazil_gene_flow/data/brazil_samples_v7.csv",
              stringsAsFactors = F)
ll2 =  ll[complete.cases(ll$LAT), ]
ll3 = ll2[ll2$lineage %in% lins, ]

#####################
# ddRAD vs sqcl
#####################

# ddRAD data, pi
dpi1 = list.files("~/Dropbox (Personal)/macroevolution/eco_IBD_oz/data/pop_gen/diversity/",
                full.names = T)
dpi2 = lapply(dpi1, read.csv, stringsAsFactors = F)
dpi3 = do.call("rbind", dpi2)

# ddRAD data, div
ddiv1 = list.files("~/Dropbox (Personal)/macroevolution/eco_IBD_oz/data/pop_gen/divergence_cov5_results/divergence/", full.names = T)
ddiv2 = lapply(ddiv1, read.csv, stringsAsFactors = F)
ddiv3 = do.call("rbind", ddiv2)

# sqcl data, pi
spi1 = list.files("~/Dropbox (Personal)/macroevolution/SqCL_sequencing/rapid_test/compare_to_SqCL/pop_gen/", 
                  pattern = "diversity.csv",
                  full.names = T)
spi2 = lapply(spi1, read.csv, stringsAsFactors = F)
spi3 = do.call("rbind", spi2)
spi4 = spi3[spi3$type == 'IND', ]

# sqcl data, div
sdiv1 = list.files("~/Dropbox (Personal)/macroevolution/SqCL_sequencing/rapid_test/compare_to_SqCL/pop_gen/", 
                   pattern = "divergence",
                   full.names = T)
sdiv2 = lapply(sdiv1, read.csv, stringsAsFactors = F)
sdiv3 = do.call("rbind", sdiv2)

############################
# comparison of pi
############################

pi = inner_join(dpi3, spi4, by = c("individual"= "ind"))
pi2 = pi %>% group_by(individual) %>% slice_max(n_inds) %>% ungroup()

cor.test(pi2$pi.x, pi2$pi.y)

pigraph = ggplot(pi2, aes(pi.x, pi.y)) + 
  geom_point() + 
  xlab(expression(pi ~ ", ddRAD")) +
  ylab(expression(pi ~ ", SqCL")) 

############################
# comparison of fst
############################

fst1 = inner_join(ddiv3, sdiv3, by = c("ind1" = "ind1", "ind2" = "ind2"))
fst2 = inner_join(ddiv3, sdiv3, by = c("ind2" = "ind1", "ind1" = "ind2"))
fst = rbind(fst1, fst2)

cor.test(fst$fst.x, fst$fst.y)

fstgraph = ggplot(fst, aes(fst.x, fst.y)) + 
  geom_point() + 
  xlab(expression(F[ST] ~ ", ddRAD")) +
  ylab(expression(F[ST] ~ ", SqCL")) 

############################
# comparison of IBD
############################

# returns the inverse of Fst
inv_fst = function(x) {
  if (x == 1) {
    return(NA)
  } else {
    return(x / (1 - x))
  }
}

# returns log distance
# note that this drops all comparisons
# between points at the same place
log_dist = function(x) {
  if (x == 0) {
    return(NA)
  } else {
    return(log(x))
  }
}

get_slopes <- function(x, fstval) {
  inv_fst1 = sapply(x[, fstval], inv_fst)
  log_dist1 = sapply(x[, "geo_dist"], log_dist)
  return(coef(lm(inv_fst1 ~ log_dist1))[2])
}

ibd = split(fst, fst$lineage)
ibd = ibd[ which(unlist(lapply(ibd, nrow)) > 2 )]
ibd1 = unlist(lapply(ibd, get_slopes, "fst.x"))
ibd2 = unlist(lapply(ibd, get_slopes, "fst.y"))

ibd3 = data.frame(ibd_dd = ibd1, ibd_sq = ibd2)
cor.test(ibd3$ibd_dd, ibd3$ibd_sq)

ibdgraph = ggplot(ibd3, aes(ibd_dd, ibd_sq)) + 
  geom_point() + 
  xlab(expression(beta[IBD] ~ ", ddRAD")) +
  ylab(expression(beta[IBD] ~ ", SqCL")) + 
  xlim(0, 1.5) +
  ylim(0, 0.75)

dd_sq = plot_grid(pigraph, fstgraph, ibdgraph, labels = c("A", "B", "C"),
          ncol = 3)
save_plot("~/Dropbox (Personal)/brazil_gene_flow/figures/ddRAD_to_SqCL.png", 
          dd_sq, ncol = 3, base_height = 3, base_width = 4)


############################
# map
############################

eco = readOGR("~/Dropbox/macroevolution/eco_IBD_oz/data/ecoregions/terr-ecoregions-TNC/tnc_terr_ecoregions.shp")

pts = SpatialPoints(ll3[,c("LON", "LAT")])
pts2 = raster::intersect(pts, eco)
eco2 = table(pts2$ECO_NAME)

worldmap <- rworldmap::getMap(resolution = "high")
sa = worldmap[which(worldmap[["continent"]] == "South America and the Caribbean"), ]

# based on this, 
# plot most common biomes
ecor = c("Caatinga", "Cerrado",
         "Mato Grosso Seasonal Forests", 
         "Atlantic Dry Forests")
cols = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")

png('~/Dropbox/brazil_gene_flow/figures/map.png', height = 7,
    width = 5, res = 200, units = "in")
par(mar = c(0,0, 0,0))
plot(NA, xlim = c(-85, -34), ylim = c(-56, 13), 
     axes = F, xlab = "", ylab = "")
plot(sa, col = "gray85", border = F, add = T)
for (i in 1:length(ecor)) {
  r = eco[eco$ECO_NAME == ecor[i], ]
  plot(r, col = alpha(cols[i], 0.5), border = F, add = T)
}
points(ll3$LON, ll3$LAT, pch = 21, bg = "white", cex = 0.75)
legend(x = -60, y = -40, legend = ecor, 
       fill = sapply(cols, alpha, 0.5),
       ncol = 1, bty = "n", cex = 0.8)
dev.off()


#####################
# comparison of divergence
#####################

d = read.csv("~/Dropbox/brazil_gene_flow/results/divergence-2020-12-02.csv", stringsAsFactors = F)
d = d[d$lineage %in% lins, ]

dist1 = ggplot(d, aes(fst, nuc_dxy)) + 
  geom_point(alpha = 0.5) + xlab(expression(F[ST])) +
  ylab(expression("nuc " ~ d[xy]))

dist2 = ggplot(d %>% filter(mt_denom > 500), aes(fst, mt_dist)) + 
  geom_point(alpha = 0.5) + xlab(expression(F[ST])) +
  ylab(expression("mt " ~ d[xy]))

dist3 = ggplot(d %>% filter(mt_denom > 500), aes(nuc_dxy, mt_dist)) + 
  geom_point(alpha = 0.5) + xlab(expression("nuc " ~ d[xy])) +
  ylab(expression("mt " ~ d[xy]))

d2 = d %>% filter(mt_denom > 500)
cor.test(d$fst, d$nuc_dxy)
cor.test(d2$fst, d2$mt_dist)
cor.test(d2$nuc_dxy, d2$mt_dist)

distcor = plot_grid(dist1, dist2, dist3, labels = c("A", "B", "C"),
                  ncol = 3)
save_plot("~/Dropbox/brazil_gene_flow/figures/genetic_divergence.png", 
          distcor, ncol = 3, base_height = 3, base_width = 4)

#####################
# comparison of slopes
#####################

d = read.csv("~/Dropbox (Personal)/brazil_gene_flow/results/IBD-2020-12-02.csv", stringsAsFactors = F)
d = d[d$lineage %in% lins, ]
dd = d %>% select(-n, -dist_type) %>% 
  pivot_wider(names_from = c(div_type), 
                  values_from = c(r2, sig, slope, r2_mt, sig_mt))

cor.test(dd$slope_inv_fst, dd$slope_nuc_dxy)
cor.test(dd$slope_inv_fst, dd$slope_mt_dist)
cor.test(dd$slope_inv_fst, dd$slope_fst)

ibda = ggplot(dd, aes(slope_inv_fst, slope_nuc_dxy)) + 
  geom_point() + 
  xlab(expression("slope of inv." ~ F[ST] ~ "and log(dist)")) +
  ylab(expression("slope of nuc." ~ d[xy] ~ "and dist."))
ibdb = ggplot(dd, aes(slope_inv_fst, slope_mt_dist)) + 
  geom_point() + 
  xlab(expression("slope of inv." ~ F[ST] ~ "and log(dist)")) +
  ylab(expression("slope of mt." ~ d[xy] ~ "and dist."))
ibdc = ggplot(dd, aes(slope_inv_fst, slope_fst)) + 
  geom_point() + 
  xlab(expression("slope of inv." ~ F[ST] ~ "and log(dist)")) +
  ylab(expression("slope of" ~ F[ST] ~ "and eco dist."))

ibdcor = plot_grid(ibda, ibdb, ibdc, labels = c("A", "B", "C"),
                    ncol = 3)
save_plot("~/Dropbox/brazil_gene_flow/figures/IBD_slopes.png", 
          ibdcor, ncol = 3, base_height = 4, base_width = 5)


#####################
# data quality
#####################

q = read.csv("~/Dropbox (Personal)//brazil_gene_flow/data/metadata/site_counts.csv",
             stringsAsFactors = F)
q = q[q$individual %in% inds, ]
q$avg_cov = q$tot_coverage / q$num_sites

q1 = ggplot(q, aes(num_loci)) + geom_histogram() +
  xlab("# of loci")
q2 = ggplot(q, aes(num_sites / 1e6)) + geom_histogram() +
  xlab("# of sites (Mb)")
q3 = ggplot(q, aes(avg_cov)) + geom_histogram() + 
  xlab("average coverage")
d = read.csv("~/Dropbox (Personal)/brazil_gene_flow/results/divergence-2020-12-02.csv",
             stringsAsFactors = F)
d1 = d[d$lineage %in% lins, ]
q4 = ggplot(d1, aes(fst_denom / 1000)) +
  geom_histogram() + 
  xlab(expression("# of sites used for" ~ F[ST] ~ "(kb)"))

mean(q$num_loci)
mean(q$num_sites) / 1e6
mean(q$avg_cov) 
mean(d$fst_denom) 

xx = ggplot(q %>% filter(num_loci > 2000), aes(num_loci, num_sites / 1e6)) + 
  geom_point(col = "#1f78b4", alpha = 0.5) + 
  xlab("number of loci") +
  ylab("number of sites (Mb)") +
  theme_classic()
save_plot("~/Desktop/quality.pdf", 
          xx, base_height = 3, base_width = 3)

qualplot = plot_grid(q1, q2, q3, q4,  labels = c("A", "B", "C", "D"),
                   ncol = 2)
save_plot("~/Dropbox (Personal)/brazil_gene_flow/figures/data_quality.png", 
          qualplot, ncol = 2, base_height = 4, base_width = 5, nrow = 2)

#####################
# ibd on phylogeny
#####################

t = read.tree("~/Dropbox/brazil_gene_flow/data/phylogeny/Tonini.squam_shl_new_Consensus_9755.targeted_species.tre")
diff = read.csv("~/Dropbox/brazil_gene_flow/results/IBD-2020-12-02.csv", stringsAsFactors = F)
diff = diff[diff$div_type == 'inv_fst' & diff$dist_type == 'ln_geo_dist', ]

rownames(diff) = diff$lineage
t4 = read.tree(text = write.tree(ladderize(t)))
diff = diff[t4$tip.label,]

png("~/Dropbox/brazil_gene_flow/figures/amount_explained.png",
    width = 5, height = 4, units = "in", res = 200)
par(mar=c(0, 0, 0, 0))
layout(matrix(c(1,2),1,2),c(0.5,0.5))
par(mar=c(4.1,1.1,1.1,0))
plot.phylo(t4, show.tip.label=F, cex=0.6)
par(mar=c(4.1,0,1.1,1.1))
cols = rep("gray25", nrow(diff))
cols[diff$sig_mt > 0.05] = "gray75"
barplot(diff$r2_mt, horiz=T, names="", space=0, col=cols, width=1, ylim=c(1, nrow(diff))-0.5, xaxt="n", border=F, xlim=c(0,1))
ticks = c(0, 0.5, 1.0)
axis(1, ticks, tck=-0.02, labels=NA)
axis(1, ticks, labels=ticks, lwd = 0, line = -0.6, las = 1)
mtext("amount explained", side=1, line=1.5)
dev.off()
mean(diff$r2_mt)

############################
# correlation btn DR & CLaDS
############################

crates = readRDS("~/Dropbox (Personal)/brazil_gene_flow/data/speciation_rates/clads_tiprates.Rdata")

d = read.csv("~/Dropbox (Personal)/brazil_gene_flow/data/speciation_rates/tonini-dr-rates-100.csv", stringsAsFactors = F, header = F)
drates = apply(d[2:101], 1, mean)
names(drates) = d$V1

df = data.frame(crates = crates,
                drates = drates[names(crates)],
                stringsAsFactors = F)

t = read.tree("~/Dropbox (Personal)/brazil_gene_flow/data/speciation_rates/tonini_100trees.txt")[[1]]
t = keep.tip(t, unique(ldf$full_tree))
t = read.tree(text = write.tree(ladderize(t)))
### full tree correlation

a = ggplot(df, aes(crates, drates)) + 
  geom_point(aes(color = log(crates))) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_gradientn(colors = gplots::rich.colors(20)) +
  xlab(expression(lambda[CLaDS])) + 
  ylab(expression(lambda[DR])) +
  theme(legend.position = "none")

#### plot clads results

ecc = rep("black", length(t$edge.length))
crates = df[ldf[match(t$tip.label, ldf$full_tree), "full_tree"], "crates"]
cols = old_colorMap(log(crates), "temperature", 100)
ecc[ match(1:Ntip(t), t$edge[,2]) ] = cols

ecd = rep("black", length(t$edge.length))
drates = df[ldf[match(t$tip.label, ldf$full_tree), "full_tree"], "drates"]
cols = old_colorMap(log(drates), "temperature", 100)
ecd[ match(1:Ntip(t), t$edge[,2]) ] = cols

plottrees <- function() {
  par(mfrow = c(1, 2), mar = c(0, 1.5, 0, 0), xpd = T)
  plot(t, edge.color = ecc, show.tip.label = F)
  pal = gplots::rich.colors(64)
  bks = assignColorBreaks(crates, method = "linear")
  addBAMMlegend_sonal(bks, pal, "horizontal", 
                      location = c(10, 50, 2, 3),
                      cex.axis = 0.5, nTicks = 0)
  
  par(mar = c(0, 0, 0, 1.5))
  plot(t, edge.color = ecd, 
       show.tip.label = F, direction = "leftwards")
  bks = assignColorBreaks(drates, method = "linear")
  addBAMMlegend_sonal(bks, pal, "horizontal", 
                      location = c(125, 165, 2, 3),
                      cex.axis = 0.5, nTicks = 0)
}

df2 = df[ldf$full_tree, ]
b = ggplot(df2, aes(crates, drates)) + 
  geom_point(aes(color = log(crates))) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_gradientn(colors = gplots::rich.colors(20)) +
  xlab(expression(lambda[CLaDS])) + 
  ylab(expression(lambda[DR])) +
  theme(legend.position = "none")
   
abtt = plot_grid(a, plottrees, b, rel_widths = c(0.25, 0.5, 0.25), nrow = 1,
                 labels = c("A", "B", "C"))
save_plot("~/Dropbox/brazil_gene_flow/figures/speciation_rate_comp.png", 
          abtt, base_width = 10, base_height = 2.5)

##############
# main figure 1
##############

b = readRDS("~/Dropbox (Personal)//brazil_gene_flow/data/speciation_rates/clads_as_bammData.Rdata")
b2 = subtreeBAMM(b, unique(ldf$gene_tree))
brates = BAMMtools::getTipRates(b)$lambda.avg

b_plot = plot.bammdata(b, tau=0.002, breaksmethod = "jenks",
                       method = "polar", lwd = 1, vtheta = 180)
cols = b_plot$palette[ cut(brates[b$tip.label], 
                           breaks = b_plot$colorbreaks,
                           include.lowest = T, right = FALSE) ]

bb = data.frame(brates = brates,
                colors = cols,
                in_br = FALSE)
bb[rownames(bb) %in% ldf$gene_tree, "in_br"] = TRUE
bb2 = bb %>% arrange(brates)
bb2$ix = seq(1:nrow(bb2))

col <- as.character(bb2$colors)
names(col) <- as.character(bb2$ix)

rateplt = ggplot(bb2, aes(ix, brates)) + 
  geom_line(col = "black", linetype = "dashed") +
  geom_point(data = bb2 %>% filter(in_br == TRUE), 
             color = "black",
             alpha = 1, size = 4.5, shape = 16) +
  geom_point(data = bb2 %>% filter(in_br == TRUE), 
             aes(color = as.character(ix)), shape = 16,
             alpha = 1, size = 3.5) +
  scale_color_manual(values=col) +
  ylab("speciation rate") + scale_y_log10() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") 
save_plot("~/Desktop/speciation_rate1.pdf", rateplt,
          base_width = 6, base_height = 2)
bb2 = bb %>% filter(in_br == TRUE)
bb2$ix = seq(1:nrow(bb2))
rateplt2 = ggplot(bb2, aes(ix, brates)) + 
  geom_point(col = "#1f78b4") +
  ylab("speciation rate") + scale_y_log10() +
  geom_hline(yintercept = min(bb$brates), linetype = "dotted") +
  geom_hline(yintercept = max(bb$brates), linetype = "dotted") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())
save_plot("~/Desktop/speciation_rate2.pdf", rateplt2,
          base_width  = 6, base_height = 2)

####################
# figure parts
####################

# focal species
pdf("~/Desktop/speciation_rates1.pdf", width = 2,
    height = 6)
par(mar = c(0, 0, 0, 0))
xx = BAMMtools::plot.bammdata(b2, lwd = 2,
                              breaksmethod = "jenks",
                   tau = 0.002, direction = "rightwards")
addBAMMlegend(xx, direction = "horizontal",
              location = c(10, 70, 1, 1.5), 
              cex.axis = 0.8, nTicks = 0, 
              tck = -0.01, padj = 1)
dev.off()

# data amount & quality
qq = q[q$num_loci > 2000, ]
qq$num_sites_Mb = qq$num_sites /  1e6
qq$num_loci_K = qq$num_loci /  1e3
qplt = ggplot(qq, aes(num_loci, num_sites_Mb)) +
  geom_point(pch = 21, fill = maincol) +
  xlab("number of loci") + 
  ylab("number of sites (Mb)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
save_plot("~/Desktop/speciation_rates2.pdf", qplt,
          base_height = 3, base_width = 3)

# map
worldmap <- rgdal::readOGR("~/Dropbox/brazil_gene_flow/data/ne_10m_land/ne_10m_land.shp")
sa2 = crop(worldmap, extent(-85, -34, -56, 13))
sa3 <- rgeos::gSimplify(sa2, 0.05)

gmap = stack("~/Dropbox/brazil_gene_flow/data/GRAY_HR_SR_OB_DR/GRAY_HR_SR_OB_DR.tif")
gmap2 = crop(gmap, extent(-85, -34, -56, 13))
gmap3 = mask(gmap2, sa3)

bc = raster("/Users/Sonal/Dropbox/brazil_gene_flow/data/envirem/current_30arcsec_climaticMoistureIndex.tif")
bc1 = crop(bc, extent(-85, -34, -56, 13))
bc2 = mask(bc1, sa3)

pdf("~/Desktop/speciation_rates3.pdf", height = 4, width = 2.7)
plot.new()
par(mar = c(0, 0, 0, 0))
plot.window(xlim = c(-85, -34), 
            ylim = c(-56, 13))
plot(gmap3,
        maxpixels = ncell(gmap3),
        col = gray.colors(1000, start = 0.9, 
                  end = 0.1, gamma = 2.2,
                  alpha = 0.9),
     legend = F, add = T)
plot(bc2, alpha = 0.5, 
     add = T, legend = FALSE)
points(ll3$LON, ll3$LAT, pch = 21, bg = "white", 
       cex = 0.6, lwd = 0.5)
dev.off()

# big overview
pdf("~/Desktop/speciation_rates4.pdf", 
    width = 4, height = 6)
mm = matrix(c(1, 2, 3, 4), 
       nrow = 4, ncol = 1, byrow = T)
layout(mm, heights = c(0.39, 0.07, 0.39, 0.15))
par(mai = c(0, 0, 0, 0), xpd = T)
b_plot = plot.bammdata(b, tau=0.002, breaksmethod = "jenks",
               method = "polar", lwd = 0.5)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
addBAMMlegend(b_plot, direction = "horizontal", 
              location = c(0.3, 0.7, 0, 0.2),
              cex.axis = 0.8, nTicks = 0, 
              tck = -0.06, padj = 0.6)
text(labels = "speciation rate", x = 0.5,
     y = 0.8, adj = 0.5)

par(mai = c(0, 0, 0, 0))
plot.bammdata(b, mask.color = alpha("gray70", 0.1),
              mask = seq(1, nrow(b$edge)),
              method = "polar",
              lwd = 0.5)
# get edges of interest
tips = match(unique(ldf$gene_tree), b$tip.label)
# edges to mask
me = which(!1:length(b$edge.length) %in% tips)
phyd = plot.bammdata_sonal(b, mask.color = alpha("black", 0),
                    mask = me, breaksmethod = "jenks",
                    tau = 0.002,
                    method = "polar",
                    show = FALSE)
plotvals = phyd$plotvals
plotcols = plotvals[[1]]
plotcols[which(plotcols != "#00000000")] = maincol
segments(plotvals[[7]]$x0,
          plotvals[[7]]$y0,
          plotvals[[7]]$x1,
          plotvals[[7]]$y1,
          # col = plotvals[[1]],
          col = plotcols,
          lwd = 2.5, lend = 2);
plotcols2 = plotvals[[6]]
plotcols2[which(plotcols2 != "#00000000")] = maincol
BAMMtools:::arc(0, 0, plotvals[[2]]$arcs[, 1],
                 plotvals[[2]]$arcs[, 2],
                 c(plotvals[[3]],
                   plotvals[[3]] + plotvals[[4]]/plotvals[[5]]),
                   border = plotcols2, lwd = 2.5)
# plot.new()

par(mar = c(1, 7, 0, 5), xpd = T)
vals = c(0.01, 0.1, 0.5)
pos = log(vals)
plot.new()
plot.window(xlim = range(bb$ix), 
            ylim = range(c(log(bb$brates), pos)))
points(bb$ix, log(bb$brates), pch = 16, 
       col = alpha("gray80", 0.1))
bb2 = bb[bb$in_br == TRUE, ]
points(bb2$ix, log(bb2$brates), pch = 16, 
       col = alpha(maincol, 0.5))
axis(2, at = pos, labels = NA, tck = -0.02)
mtext(side = 2, l = 0.3, 
      text = vals, at = pos, las = 2, cex = 0.5)
mtext(side = 2, l = 1.5, 
      text = "speciation rate", cex = 0.5)
dev.off()


##############
# main figure 2 - ibd
##############

sps = c("Bothrops_moojeni", 
        "Micrablepharus_atticolus")
r = c("~//Dropbox/brazil_gene_flow/data/distributions/try2/Bothrops_moojeni.shp",
      "~//Dropbox/brazil_gene_flow/data/distributions/try2/Micrablepharus_atticolus.shp")
r2 = lapply(r, readOGR)
div = read.csv("~/Dropbox/brazil_gene_flow/results/divergence-2020-12-02.csv", stringsAsFactors = F)
div2 = div[div$lineage %in% sps, ]
div2 = div2[is.finite(div2$ln_geo_dist), ]

hyp = stack("~/Dropbox/brazil_gene_flow/data/HYP_HR_SR_W_DR/HYP_HR_SR_W_DR.tif")
hyp1 = crop(hyp, extent(c(-85, -34, -56, 13)))
hyp2 = raster::mask(hyp1, sa3)

ibdplts = vector("list", length = length(sps))

pdf("~/Desktop/ibd1.pdf", height = 5, width = 3)
par(mar = c(0.5, 0.5, 0.5, 0.5), mfrow = c(2, 1))
for (i in 1:length(sps)) {
  plot(NA, xlim = c(-85, -34), ylim = c(-30, 13), 
       axes = F, xlab = "", ylab = "")
  # plot(sa, col = "gray90", border = F, add = T)
  plotRGB(hyp2, add = T,maxpixels = ncell(hyp))

  ll4 = ll3[ll3$lineage == sps[i], ]
  
  range = r2[[i]]
  # range = r[which(rnames == sps[i]), 1]
  plot(range, add = T, col = alpha("white", 0.35),
         border = F)
  lat = jitter(ll4$LAT)
  lon = jitter(ll4$LON)
  points(lon, lat, bg = maincol, 
         pch = 21, cex = 1)
  
  div3 = div2 %>% filter(lineage == sps[i])
  ibdplts[[i]] = ggplot(div3, 
         aes(ln_geo_dist, inv_fst)) +
    geom_smooth(method = "lm", fill = NA, 
                col = "gray50", lwd = 0.7) +
    geom_point(shape = 21, fill = maincol, size = 2) +
    xlab("log(geo dist, km)") +
    ylab(expression(F[ST] / (1 - F[ST]) )) +
    xlim(range(div2$ln_geo_dist)) +
    ylim(range(div2$inv_fst)) +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}
dev.off()  

ibdg = ibdplts[[1]] / ibdplts[[2]] 
save_plot("~/Desktop/ibd2.pdf", ibdg, ncol = 1,
          base_height = 5, base_width = 3)


##############
# main figure 3 - correlation
##############

d = read.csv("~/Dropbox (Personal)//brazil_gene_flow/results/IBD-2020-12-02.csv", stringsAsFactors = F)
d = d[d$lineage %in% lins,  ]
d2 = d[d$div_type == 'inv_fst', ]
d2$log_slope = log(d2$slope)

crates = readRDS("~/Dropbox (Personal)//brazil_gene_flow/data/speciation_rates/clads_tiprates.Rdata")
crates2 = crates[ ldf[ match(d2$lineage, ldf$lineage), "full_tree"] ]
  
d3 = cbind(d2, crates = crates2)
d3$log_crate = log(d3$crates)
d4 = d3[complete.cases(d3$log_slope), ]

t = read.tree("~/Dropbox (Personal)/brazil_gene_flow/data/phylogeny/Tonini.squam_shl_new_Consensus_9755.targeted_species.tre")
t4 = keep.tip(t, d4$lineage)
t4 = read.tree(text=write.tree(ladderize(t4)))

snake = c("Trilepida_koppesi", "Tantilla_melanocephala")
snake = getMRCA(t4, snake)
snakes = t4$tip.label[ Descendants(t4, snake, type = "tips")[[1]] ]
d4$lizard = "lizard"
d4[d4$lineage %in% snakes, "lizard"] = "snake"
d4$color = maincol2
d4[d4$lizard == 'snake', "color"] = maincol

slopelog = d4[match(t4$tip.label, d4$lineage), "log_slope"]
names(slopelog) = t4$tip.label

map = contMap(t4, slopelog, plot=FALSE)
map = setMap(map, colors=brewer.pal(11, "YlGn"),
             space="Lab")

d5 = d4[match(t4$tip.label, d4$lineage), ]

pdf("~/Desktop/correlation1.pdf", height = 3, width = 4)
mm2 = matrix(c(1, 2, 3), byrow = F, ncol = 3)
layout(mm2, widths = c(0.4, 0.3, 0.3))
plot(map, fsize=c(0,0.9), legend = F,
     lwd=2.1, ftype="off", outline=F, mar=c(1,0,0,0),
     ylim=c(1-0.09*(Ntip(map$tree)-1), Ntip(map$tree)), xpd=NA)
nodelabels("", snake, frame = "none", pch = 16, 
           cex = 2, col = "gray70")
nodelabels("S", snake, frame = "none", cex = 0.7)
add.color.bar(80, map$cols, title=expression(log ~ beta[IBD]),
              lims=map$lims, digits=1, prompt=FALSE, x=10,
              y=1, lwd=4, fsize = 1.3,
              subtitle="")
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
ticks = pretty(d5$log_slope, 3)
par(mar = c(1, 0.8, 0, 0.8))
plot(NULL, xlim=range(ticks), 
     ylim=lastPP$y.lim, 
     xlab="", ylab="", axes=F)
for (i in 1:nrow(d5)) {
  x1 = min(ticks)
  x2 = d5[i, "log_slope"]
  lines(c(x1, x2), c(i, i), lwd=0.3, col="gray60")
}
points(d5$log_slope, 1:nrow(d5), 
       pch= 21, 
       bg = d5$color, cex = 0.8)
axis(1, ticks, tck=-0.02, labels=NA, line = -1.8)
mtext(side = 1, l = -1.7, text = ticks, 
      at = ticks, las = 1, cex = 0.5)
mtext(expression(log ~ beta[IBD]), 
      side = 1, line = 0, cex = 0.8)
ticks = pretty(d5$log_crate, 3)
plot(NULL, xlim=range(ticks), 
     ylim=lastPP$y.lim, 
     xlab="", ylab="", axes=F)
for (i in 1:nrow(d5)) {
  x1 = min(ticks)
  x2 = d5[i, "log_crate"]
  lines(c(x1, x2), c(i, i), lwd=0.3, col="gray60")
}
points(d5$log_crate,
       1:nrow(d5), 
       pch=21, bg = d5$color, cex = 0.8)
axis(1, ticks, tck=-0.02, labels=NA, line = -1.8)
mtext(side = 1, l = -1.7, text = ticks, 
      at = ticks, las = 1, cex = 0.5)
mtext(expression("log speciation rate"), side = 1, 
      line = 0, cex = 0.8)
dev.off()


crplt = ggplot(d4, aes(log_slope, log_crate)) + 
  geom_point(pch = 21, size = 2, aes(fill = lizard)) +
  ylab("log speciation rate") +  
  xlab(expression(log ~ beta[IBD])) + 
  scale_fill_manual(values = c(maincol2, maincol)) +
  theme_classic()
save_plot("~/Desktop/correlation2.pdf", crplt, 
         base_height = 2.6, base_width = 4)


##############
# main figure 4 - snakes vs lizards
##############

a = ggplot(d4, aes(lizard, log_slope)) + 
  geom_boxplot(aes(fill = lizard)) + 
  ylab(expression(log ~ beta[IBD])) + xlab("") + 
  scale_fill_manual(values = c(maincol2, maincol)) +
  theme_classic() +
  theme(legend.position = "none")
b = ggplot(d4, aes(lizard, log_crate)) + 
  geom_boxplot(aes(fill = lizard)) + 
  ylab("log speciation rate") + xlab("") + 
  scale_fill_manual(values = c(maincol2, maincol)) +
  theme_classic() +
  theme(legend.position = "none")

library(readxl)
m = read_xlsx("~/Dropbox (Personal)/brazil_gene_flow/data/Meiri_data/feldman_et_al._2015_lepidosaur_body_sizes_appendix_s1.xlsx")
tips = gsub("_", " ",
            ldf[match(d4$lineage, ldf$lineage), "full_tree"])
d4$length = as.numeric(pull(m[match(tips, m$binomial), "max length (mm)"]))
d4$mass = as.numeric(pull(m[match(tips, m$binomial), "mass (g)"]))

er <- function(L, M){
  return(0.5*L^(3/2)*sqrt(pi/M))
}

d4$ratio = log(d4$mass) / log(d4$length)
d4$er = er(d4$length, d4$mass)
d4$log_er = log(d4$er)

c = ggplot(d4, aes(log_er, log_slope)) + 
  geom_point(shape = 21, aes(fill = lizard)) +
  xlab("log body elongation") +
  ylab(expression(log ~ beta[IBD])) +
  scale_fill_manual(values = c(maincol2, maincol)) +
  theme_classic() +
  theme(legend.position = "none")

d = ggplot(d4, aes(x = log_slope, fill = lizard)) + 
  geom_density(alpha = 0.6, lwd = 0.2) + 
  theme_void() + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c(maincol2, maincol)) +
  coord_flip()

e = ggplot(d4, aes(x = log_er, fill = lizard)) + 
  geom_density(alpha = 0.6, lwd = 0.2) + 
  theme_void() + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c(maincol2, maincol)) 

half1 = e + plot_spacer() + c + d + 
  plot_layout(ncol = 2, nrow = 2, 
              widths = c(4, 1), 
              heights = c(1, 4))

cab = (c | plot_spacer() | (a / b)) + 
  plot_layout(widths = c(1.5, 0.1, 1))
save_plot("~/Desktop/snakes_vs_lizards.pdf", cab,
          base_height = 3.2, base_width = 4.2)

library(nlme)
d4t = d4[match(t4$tip.label, d4$lineage), ]
xx = gls(log_slope ~ log_er, correlation = corBrownian(phy = t4),
    data = d4t, method = "ML")

##############
# power
##############

p = readRDS("~/Dropbox (Personal)/brazil_gene_flow/results/power.Rds")
p1 = p[[2]] %>% gather("test", "power", -cor) %>%
  filter(cor >= 0.1, cor < 1)
p1[p1$test == "power_es", "test"] = "ES-sim"
p1[p1$test == "power_gls", "test"] = "PGLS"

a = ggplot(p1, aes(cor, power)) +
  geom_line(aes(color = test)) +
  geom_point(aes(fill = test), shape = 21) +
  ylab("power") +
  xlab("simulated correlation") +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = c("#2b83ba", "#d7191c")) +
  scale_colour_manual(values = c("#2b83ba", "#d7191c"))
  
p2 = p[[1]]
p3 = p2 %>% dplyr::select(cor_gls, cor_es, cor) %>%
  gather("type", "correlation", -cor) %>% 
  filter(cor >= 0.1, cor < 1)

p3[p3$type == "cor_es", "type"] = "ES-sim"
p3[p3$type == "cor_gls", "type"] = "PGLS"

b = ggplot(p3, aes(as.character(cor), correlation, col = type)) +
  geom_point(position = position_jitterdodge(), 
             alpha=0.3, size= 0.5) +
  geom_boxplot(outlier.shape = NA, fill = NA) + 
  xlab("simulated correlation") +
  ylab("estimated correlation") +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = c("#2b83ba", "#d7191c")) +
  scale_colour_manual(values = c("#2b83ba", "#d7191c"))

d = read.csv("~/Dropbox (Personal)/brazil_gene_flow/results/power.csv",
             stringsAsFactors = F)
d2 = d[d$estrho > 0.03 & d$estrho < 0.05,]

c = ggplot(d2, aes(actualrho)) + 
  geom_histogram(binwidth = 0.025) +
  geom_vline(xintercept = median(d2$actualrho), col = "red") +
  xlab("value of simulated correlation")

ab = plot_grid(a, c, labels = c("A", "B"))
save_plot("~/Dropbox (Personal)/brazil_gene_flow/figures/power.png",ab,
          base_height = 3, base_width = 10)
