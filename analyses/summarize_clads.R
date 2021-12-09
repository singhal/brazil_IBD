library(ape)
library(BAMMtools)

source("~/Dropbox/scripts/brazil/clads-to-bamm-v2.R")

### get list of tips that are gene only
t = read.tree("~/Dropbox/brazil_gene_flow/data/speciation_rates/Squamata_Tree_Methods_Data/squam_shl_new.tre")
keep = t$tip.label

# t2 = read.tree("~/Dropbox/brazil_gene_flow/data/speciation_rates/tonini_100trees.txt")
# t3 = lapply(t2, keep.tip, keep)
# t4 = lapply(t3, check_and_fix_ultrametric)
# t5 = unique.multiPhylo(t4)

### read in CLaDS objects
clads = list.files("~/Dropbox/brazil_gene_flow/data/speciation_rates/ClaDS/", full.names = T, pattern = "Rdata")

clads2 = vector("list", length(clads))
for (i in 1:length(clads)) {
   load(clads[i])
   clads2[[i]] = CladsOutput
}

crates = lapply(clads2, function(x) {return(x$lambdatip_map)})
tips = sort(clads2[[1]]$tree$tip.labels)
crates2 = vector("list", length(crates))
for (i in 1:length(crates)) {
  names(crates[[i]]) = clads2[[i]]$tree$tip.labels
  crates2[[i]] =  crates[[i]][ tips ]
}
crates3 = do.call(rbind, crates2)
crates4 = apply(crates3, 2, mean)
saveRDS(crates4, "~/Dropbox/brazil_gene_flow/data/speciation_rates/clads_tiprates.Rdata")

#### double check that keep tips match to tree tips
table(keep %in% clads2[[1]]$tree$tip.labels)

#### read in event data
####### subsample
####### make it look like an eventdata dataframe
####### add in a generation
####### write to file
for (i in 1:length(clads2)) {
  gen = gsub("^.*tree", "", clads[i])
  gen = as.numeric(gsub(".Rdata", "", gen))
  outf = paste0("~/Dropbox/brazil_gene_flow/data/speciation_rates/ClaDS/bamm_event_data", gen, ".csv")
  if (!file.exists(outf)) {
    cat(i, "\n")
    ed = clads_to_eventdata(clads2[[i]])
    ed_sub = subtreeBAMM(ed, keep)
    eddf = ed_sub$eventData[[1]]
    tax = assignSpanningSet(ed_sub)
    
    eddf2 = data.frame(generation = gen,
                       leftchild = tax[eddf$node, 1],
                       rightchild = tax[eddf$node, 2],
                       abstime = eddf$time,
                       lambdainit = eddf$lam1,
                       lambdashift = eddf$lam2,
                       muinit = eddf$mu1,
                       mushift = eddf$mu2,
                       stringsAsFactors = F)
    write.csv(eddf2, outf, quote = F, row.names = F) 
  }
}

#### read in all event data
ed = clads_to_eventdata(clads2[[99]])
ed_sub = subtreeBAMM(ed, keep)

d = list.files("~/Dropbox/brazil_gene_flow/data/speciation_rates/ClaDS/", full.names = T, pattern = "csv")
d1 = lapply(d, read.csv, stringsAsFactors = F)
d2 = do.call(rbind, d1)
d3 = d2 %>% arrange(generation, abstime)
ed2 = getEventData(ed_sub, eventdata = d3, burnin=0)
saveRDS(ed2, "~/Dropbox/brazil_gene_flow/data/speciation_rates/clads_as_bammData.Rdata")

#####
# confirm results are same
#####

b = readRDS("~/Dropbox/brazil_gene_flow/data/speciation_rates/clads_as_bammData.Rdata")
brates = getTipRates(b)$lambda.avg

dd = data.frame(brates = brates,
                crates = crates4[names(brates)])
cor.test(dd$crates, dd$brates)
