library(ape)

ll = read.csv("~/Dropbox/brazil_gene_flow/data/brazil_samples_v7.csv",
              stringsAsFactors = F)
ll2 =  ll[complete.cases(ll$LAT), ]
lins = table(ll2$lineage)
lins = names(lins[ lins > 3])
lins = lins[lins != "Phrynonax_poecilonotus"]

t = read.tree("~/Dropbox/macroevolution/published_data/Tonini2016_Squamata_Tree_Methods_Data/squam_shl_new_Consensus_9755.tre")
t2 = phytools::bind.tip(t, "Oxyrhopus_trigeminus_2", 
                        where = which(t$tip.label=="Oxyrhopus_trigeminus"),
                        position = 1,
                        edge.length = 1)

names(lins) = lins
lins[!lins %in% t2$tip.label]

# clean up tips not in tree
names(lins)[which(lins == "Bachia_bresslaui_1")] = "Bachia_bresslaui"
names(lins)[which(lins == "Copeoglossum_nigropunctatum_1")] = "Copeoglossum_nigropunctatum"
# assign to this tip (for now?)
names(lins)[which(lins == "Spilotes_sp")] = "Spilotes_sulphureus"
names(lins)[which(lins == "Taeniophallus_occipitalis_2")] = "Taeniophallus_occipitalis"
names(lins)[which(lins == "Trilepida_brasiliensis_2")] = "Trilepida_brasiliensis"
names(lins)[which(lins == "Vanzosaura_rubricauda_2")] = "Vanzosaura_rubricauda"
names(lins)[which(lins == "Xenodon_merremii")] = "Xenodon_merremi"
names(lins)[which(lins == "Oxyrhopus_trigeminus_1")] = "Oxyrhopus_trigeminus"

tt = keep.tip(t2, names(lins))
tt$tip.label = lins[tt$tip.label]

t3 = t2
t3$tip.label = ifelse(t2$tip.label %in% names(lins),
                      lins[t2$tip.label], t2$tip.label)

write.tree(tt, "~/Dropbox/brazil_gene_flow/data/phylogeny/Tonini.squam_shl_new_Consensus_9755.targeted_species.tre")
write.tree(t3, "~/Dropbox/brazil_gene_flow/data/phylogeny/Tonini.squam_shl_new_Consensus_9755.renamed.tre")
