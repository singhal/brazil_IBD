calc_het <- function(h) {
  het = pull(h %>% filter(type == "s") %>%
               summarise(sum(var) / sum(denom)))
  denom = pull(h %>% filter(type == "s") %>%
                 summarise(sum(denom)))
  ind = h[1, "individual"]
  return(c(ind, het, denom))
}

maincol2 = '#e41a1c'
maincol = '#377eb8'

d = read.csv("~/Dropbox (Personal)/brazil_gene_flow/results/IBD-2020-12-02.csv", stringsAsFactors = F)
d1 = d %>% dplyr::filter(div_type == 'inv_fst')
d1$log_slope = log(d1$slope)
d1[which(d1$lineage == "Vanzosaura_rubricauda_2"), "lineage"] = "Vanzosaura_savanicola"
d1 = d1[d1$lineage != 'Phrynonax_poecilonotus', ]

ll1 = ll[which(ll$lineage %in% d1$lineage), ]
hetfiles = paste0("~/Dropbox (Personal)/squamate_tree_analyses/diversity/",
                  ll1$sample, "_diversity.csv")
h1 = lapply(hetfiles, read.csv)
h2 = lapply(h1, calc_het)
h3 = data.frame(do.call("rbind", h2))
names(h3) = c("sample", "pi", "denom")
h3$lineage = ll1[match(h3$sample, ll1$sample), "lineage"]
h3$pi = as.numeric(h3$pi)
h3$denom = as.numeric(h3$denom)
h4 = h3 %>% filter(denom > 4e4) %>% group_by(lineage) %>%
  summarize(pi = mean(pi)) %>% ungroup()
d1$pi = pull(h4[match(d1$lineage, h4$lineage), "pi"])

t = read.tree("~/Dropbox (Personal)/brazil_gene_flow/data/phylogeny/Tonini.squam_shl_new_Consensus_9755.targeted_species.tre")
t$tip.label[which(t$tip.label == "Vanzosaura_rubricauda_2")] = "Vanzosaura_savanicola"
t2 = keep.tip(t, d1$lineage)
d2 = d1[match(t2$tip.label, d1$lineage), ]
phytools::phylANOVA(t2, d2$critter, d2$pi)

d1$critter = "lizard"
t = read.tree("~/Dropbox (Personal)/brazil_gene_flow/data/phylogeny/Tonini.squam_shl_new_Consensus_9755.targeted_species.tre")
snake = c("Trilepida_koppesi", "Tantilla_melanocephala")
snake = getMRCA(t, snake)
snakes = t$tip.label[ phangorn::Descendants(t, snake, type = "tips")[[1]] ]
d1[d1$lineage %in% snakes, "critter"] = "snake"
a = ggplot(d1, aes(critter, pi)) + geom_boxplot(aes(fill = critter))  +
  xlab("") + ylab(expression(pi)) +
  scale_fill_manual(values = c(maincol2, maincol)) + 
  theme_classic() + 
  theme(legend.position = "none")


d = readxl::read_xlsx("~/Dropbox (Personal)/brazil_gene_flow/response/dispersal.xlsx")
d[d$critter == "lizards", "critter"] = "lizard"
d[d$critter == "snakes", "critter"] = "snake"
d$dispersal = as.numeric(d$`dispersal (km)`)
b = ggplot(d, aes(critter, dispersal)) + 
  geom_boxplot(aes(fill = critter))  +
  xlab("") + ylab("dispersal length") +
  scale_y_log10() +
  scale_fill_manual(values = c(maincol2, maincol)) + 
  theme_classic() + 
  theme(legend.position = "none")

ab = plot_grid(a, b, labels = c("A", "B"))
save_plot("~/Desktop/elongation_IBD.png", ab, base_height = 3, 
          base_width = 8)

# d[d$sp1 == "Eumeces skiltonianus", "sp1"] = "Plestiodon skiltonianus"
# d[d$sp1 == "Cnemidophorus punctilinealis", "sp1"] = "Aspidoscelis tigris"
# d[d$sp1 == "Hoplodactylus maculatus", "sp1"] = "Woodworthia maculatus"
# d[d$sp1 == "Rhinoplocephalus nigrescens", "sp1"] = "Cryptophis nigrescens"
# d[d$sp1 == "Lacerta vivipara", "sp1"] = "Zootoca vivipara"
# 
# er <- function(L, M){
#   return(0.5*L^(3/2)*sqrt(pi/M))
# }
# 
# m = read_xlsx("~/Dropbox (Personal)/brazil_gene_flow/data/Meiri_data/feldman_et_al._2015_lepidosaur_body_sizes_appendix_s1.xlsx")
# d1 = d
# d1$body_length = as.numeric(pull(m[match(d1$sp1, m$binomial), "max length (mm)"]))
# d1$body_mass = as.numeric(pull(m[match(d1$sp1, m$binomial), "mass (g)"]))
# 
# d1$er = er(d1$body_length, d1$body_mass)
# plot(log(d1$er), log(d1$dispersal))
