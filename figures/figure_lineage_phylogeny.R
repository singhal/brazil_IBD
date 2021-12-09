library(ape)
library(RColorBrewer)

#######
#	called by colorClade(...)
#
zip = function(x,y) {
  ret = numeric(2*length(x));
  i = 1; j =1;
  while(i <= length(ret)) {
    ret[i] = x[j];
    i = i+1;
    ret[i] = y[j];
    i = i+1;
    j = j+1;
  }
  return(ret);
}

########
# polar tree plotting function that returns coordinates of edges
#
#
polar.phylo = function(phy, tip_names, vtheta = 1, rbf = 0.001, labels = FALSE, lwd = 1, edge.color = 1, arc.color = 1, cex = 1)
{
  phy = BAMMtools:::getStartStopTimes(phy)
  tH = max(phy$end);
  rb = tH * rbf;
  ret = BAMMtools:::setPolarTreeCoords(phy, vtheta, rbf);
  
  x0 = ret$segs[, 1];
  x1 = ret$segs[, 3];
  y0 = ret$segs[, 2];
  y1 = ret$segs[, 4];
  ofs = 0;
  if (labels) {
    ofs = max(nchar(tip_names) * 0.03 * cex);
  }
  
  plot.new();
  plot.window(xlim = c(-1, 1) + c(-rb, rb) + c(-ofs, ofs), ylim = c(-1, 1) + c(-rb, rb) + c(-ofs, ofs), asp = 1);
  segments(x0, y0, x1, y1, col = edge.color, lwd = lwd, lend = 2);
  BAMMtools:::arc(0, 0, ret$arcs[, 1], ret$arcs[, 2], c(rb, rb + phy$end/tH), border = rep(arc.color, nrow(ret$arcs)), lwd = lwd);
  
  invisible(ret);
}

#######
#	Arguments:
#		phy     -  phylo object 
#		p		-  plotting coordinates of phylogeny. returned invisibly from plot.bammdata
#		node1   -  left bound tip node (may be character string) 
#		node2   -  right bound tip node (may be character string)
#		col     -  fill and border color of polygon
#		alpha   -  transparency of color
colorClade = function(phy, p, nodes, col, border, alpha, name, cex) {
  n2 = max(nodes);
  n1 = min(nodes);
  seq.nod = .Call("seq_root2tip",phy$edge,phy$Nnode+1,phy$Nnode,PACKAGE="ape");
  nmrca = BAMMtools:::getmrca(phy, n1, n2);
  left = seq.nod[[n2]];
  left = left[which(left == nmrca):length(left)];
  right = seq.nod[[n1]];
  right = right[which(right == nmrca):length(right)];
  if (n1 == n2) {
    tips = seq.int(n1, n2-1);	
  }
  else if (n1 == n2+1) {
    tips = seq.int(n1, n2-1);
  }
  else {
    tips = seq.int(n1, n2-1);
  }
  
  lc = p[rownames(p) %in% left,];
  rc = p[rownames(p) %in% right,];
  tc = p[rownames(p) %in% tips,];
  
  if (is.null(dim(tc))) {
    tc = matrix(tc, nrow = 1);
  }
  
  xv = c(zip(lc[-1,1],lc[-1,3]), tc[nrow(tc):1,3], zip(rc[nrow(rc):2,3], rc[nrow(rc):2,1]));
  yv = c(zip(lc[-1,2],lc[-1,4]), tc[nrow(tc):1,4], zip(rc[nrow(rc):2,4], rc[nrow(rc):2,2]));
  polygon(xv, yv, col=BAMMtools:::transparentColor(col,alpha), border= border);
  
  angle = (mean(nodes) * 2 * 3.14159) / length(phy$tip.label)
  text(x=1.05 * cos(angle),
       y = 1.05 * sin(angle), name, srt=(angle * 360 / (2 * 3.14159)), 
       adj=c(0,0.5), cex=cex, font = 3)
}


ldf = readRDS("~/Dropbox (Personal)/brazil_gene_flow/data/sample_map.Rds")

ll = read.csv("~/Dropbox (Personal)//brazil_gene_flow/data/brazil_samples_v8.csv",
              stringsAsFactors = F)
ll2 =  ll[complete.cases(ll$LAT), ]
lins = table(ll2$lineage)
lins = names(lins[ lins > 3])
ll3 = ll2[ll2$lineage %in% lins, ]

t = read.tree("~/Dropbox (Personal)/brazil_gene_flow/data/prelim_tree/ExaML_result.concat_ind0.05_loci0.6_all_n4796.rooted.no_outs.tre")
ll4 = ll3[ll3$sample %in% t$tip.label, ]
t2 = keep.tip(t, ll4$sample)
t3 = chronopl(t2, 0.01)
t3 = read.tree(text = write.tree(ladderize(t3)))

colors2 = rep(brewer.pal(12, "Paired"), 20)

pdf("~/Desktop/lineage_phylogeny.pdf", width=3, height=3)
par(mar=c(0, 0, 0, 0))

cexval = 0.4
ret = polar.phylo(t3, lwd=0.2, vtheta=0, labels=FALSE, cex=cexval)

lins2 = unique(ll4[match(t3$tip.label, ll4$sample), "lineage"])
groups = split(ll4, ll4$lineage)
for (i in 1:length(lins)) {
  group = lins2[i]
  inds = groups[[group]]$sample
  nodes = match(inds, t3$tip.label)
  colorClade(t3, ret$seg, nodes, colors2[i], NA, 0.5, "", cexval)
}
dev.off()



pdf("~/Desktop/lineage_phylogeny2.pdf", width=6, height=6)
par(mar=c(5, 5, 5, 5), xpd = TRUE)

cexval = 0.4
ret = polar.phylo(t3, lwd=0.7, vtheta=0, labels=FALSE, cex=cexval)

lins2 = unique(ll4[match(t3$tip.label, ll4$sample), "lineage"])
groups = split(ll4, ll4$lineage)
for (i in 1:length(lins)) {
  group = lins2[i]
  inds = groups[[group]]$sample
  nodes = match(inds, t3$tip.label)
  group2 = gsub("_", " ", group)
  colorClade(t3, ret$seg, nodes, colors2[i], NA, 0.5, group2, cexval)
}
dev.off()