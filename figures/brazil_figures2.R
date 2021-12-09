library(BAMMtools)

willeerd.tiplabels <- function (text, tip, adj = c(0.5, 0.5), radial.adj=1, frame = "rect", pch = NULL, thermo = NULL, pie = NULL, piecol = NULL, col = "black", bg = "yellow", horiz = FALSE, width = NULL, height = NULL, ...) 
{
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  if (missing(tip)) 
    tip <- 1:lastPP$Ntip
  lastPP$xx <- lastPP$xx * radial.adj
  lastPP$yy <- lastPP$yy * radial.adj
  XX <- lastPP$xx[tip]
  YY <- lastPP$yy[tip]
  BOTHlabels(text, tip, XX, YY, adj, frame, pch, thermo, pie, 
             piecol, col, bg, horiz, width, height, ...)
}

b = readRDS("~/Dropbox (Personal)/brazil_gene_flow/data/speciation_rates/clads_as_bammData.Rdata")
b2 = subtreeBAMM(b, unique(ldf$gene_tree))

randtips = sample(b$tip.label, round(0.5 * length(b$tip.label), 0))
randtips2 = unique(c(randtips, ldf$gene_tree))

b_sub = subtreeBAMM(b, randtips2)
brates = BAMMtools::getTipRates(b)$lambda.avg
bb = data.frame(ix = seq(1, length(brates)),
                brates = brates[order(brates)],
                in_br = FALSE)
bb[rownames(bb) %in% ldf$gene_tree, "in_br"] = TRUE

tips = b_sub$tip.label
sizes = rep(0, length(tips))
sizes[ which(tips %in% ldf$gene_tree) ] = 1

fam = read.csv("~/Dropbox (Personal)/squamate_tree/species_family.csv", 
               stringsAsFactors = F)
fams = fam[match(b_sub$tip.label, fam$Species), "Family"]
names = unlist(lapply(fams, function(x) { substr(x, 1, 1) }))
names2 = names
names2[sample(1:length(names), 800)] = ""
# names = unlist(lapply(b_sub$tip.label, function(x) { substr(x, 1, 1) }))
# names2 = gsub("\\S+_", "", b_sub$tip.label)
# names2 = unlist(lapply(names2, function(x) { substr(x, 1, 1) }))
# names3 = mapply(function(x, y) {paste(x, y, sep = "")}, names, names2)
# names3[ which(!tips %in% ldf$gene_tree) ] = "" 

pdf("~/Desktop/speciation_rates.pdf", width = 8.5, height = 8.5)
par(mai = c(0, 0, 0, 0), xpd = T)
b_plot = plot.bammdata(b_sub, tau=0.002, breaksmethod = "jenks",
                       method = "polar", lwd = 1, vtheta = 180)
cols = b_plot$palette[ cut(brates[b_sub$tip.label], 
                           breaks = b_plot$colorbreaks,
                           include.lowest = T, right = FALSE) ]
# willeerd.tiplabels("", seq(1, length(tips)), frame = "none",
#                    radial.adj = 1.03,
#                    pch= 21, cex = sizes,
#                    bg = cols)
willeerd.tiplabels(names2, seq(1, length(tips)), frame = "none",
                   radial.adj = 1.03, cex = 0.3)
# plot.new()
# plot.window(xlim = c(0, 1), ylim = c(0, 1))
# addBAMMlegend(b_plot, direction = "horizontal", 
#              location = c(0.3, 0.7, 0, 0.2),
#              cex.axis = 0.8, nTicks = 0, 
#              tck = -0.06, padj = 0.6)
# text(labels = "speciation rate", x = 0.5,
#     y = 0.8, adj = 0.5)
dev.off()



worldmap <- rgdal::readOGR("~/Dropbox (Personal)/brazil_gene_flow/data/ne_10m_land/ne_10m_land.shp")
sa2 = crop(worldmap, extent(-85, -34, -56, 13))
sa3 <- rgeos::gSimplify(sa2, 0.05)

hyp = stack("~/Dropbox (Personal)/brazil_gene_flow/data/HYP_HR_SR_W_DR/HYP_HR_SR_W_DR.tif")
hyp1 = crop(hyp, extent(c(-85, -34, -56, 13)))
hyp2 = raster::mask(hyp1, sa3)

pdf("~/Desktop/speciation_rates2.pdf", height = 4, width = 2.7)
par(mar = c(0, 0, 0, 0))
plot(NA, xlim = c(-85, -34), ylim = c(-56, 13), 
     axes = F, xlab = "", ylab = "")
plotRGB(hyp2, add = T,maxpixels = ncell(hyp))
points(ll3$LON, ll3$LAT, pch = 21, bg = "white", 
       cex = 0.6, lwd = 0.5)
dev.off()