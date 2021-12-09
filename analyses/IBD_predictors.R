library(nlme)
library(ape)
library(dplyr)
library(rgdal)
library(phytools)
library(ggplot2)
library(cowplot)
library(readxl)
library(patchwork)
theme_set(theme_cowplot())

# basic info
ldf = readRDS("~/Dropbox (Personal)/brazil_gene_flow/data/sample_map.Rds")
ll = read.csv("~/Dropbox (Personal)/brazil_gene_flow/data/brazil_samples_v8.csv",
              stringsAsFactors = F)

# speciation rates
crates = readRDS("~/Dropbox (Personal)/brazil_gene_flow/data/speciation_rates/clads_tiprates.Rdata")

# tree
t = read.tree("~/Dropbox (Personal)/brazil_gene_flow/data/phylogeny/Tonini.squam_shl_new_Consensus_9755.targeted_species.tre")
t$tip.label[which(t$tip.label == "Vanzosaura_rubricauda_2")] = "Vanzosaura_savanicola"

# ibd
d = read.csv("~/Dropbox (Personal)/brazil_gene_flow/results/IBD-2020-12-02.csv", stringsAsFactors = F)
d1 = d %>% dplyr::filter(div_type == 'inv_fst')
d1$log_slope = log(d1$slope)
d1[which(d1$lineage == "Vanzosaura_rubricauda_2"), "lineage"] = "Vanzosaura_savanicola"
d1 = d1[d1$lineage != 'Phrynonax_poecilonotus', ]

#############
# get ranges
#############
r = readOGR("~/Dropbox (Personal)/inornatus_gr/Meiri_etal_ranges/GARD1.1_dissolved_ranges/modeled_reptiles.shp")
rnames = r[[1]]
rnames = gsub(" ", "_", rnames)

spnames = d1$lineage
names(spnames) = spnames
spnames[! spnames %in% rnames]

spnames["Xenodon_merremii"] = "Xenodon_merremi"
spnames["Trilepida_brasiliensis_2"] = "Trilepida_brasiliensis"
spnames["Taeniophallus_occipitalis_2"] = "Taeniophallus_occipitalis"
spnames["Oxyrhopus_trigeminus_2"] = "Oxyrhopus_trigeminus"
spnames["Spilotes_sp"] = "Spilotes_sulphureus"
spnames["Copeoglossum_nigropunctatum_1"] = "Copeoglossum_nigropunctatum"
spnames["Oxyrhopus_trigeminus_1"] = "Oxyrhopus_trigeminus"
spnames["Bachia_bresslaui_1"] = "Bachia_bresslaui" 
 
ranges = rep(NA, length(spnames))
for (i in 1:length(ranges)) {
   rnum = which(rnames == spnames[i])
   range = r[rnum, 1]
   range = spTransform(range, CRS("+init=epsg:3395"))
   ranges[i] = rgeos::gArea(range) / 1000 / 1000
}
d1$range_area = ranges

#############
# body size & mass & er
#############
m = read_xlsx("~/Dropbox (Personal)/brazil_gene_flow/data/Meiri_data/feldman_et_al._2015_lepidosaur_body_sizes_appendix_s1.xlsx")
tips = gsub("_", " ",
            ldf[match(d1$lineage, ldf$lineage), "full_tree"])
d1$body_length = as.numeric(pull(m[match(tips, m$binomial), "max length (mm)"]))
d1$body_mass = as.numeric(pull(m[match(tips, m$binomial), "mass (g)"]))

er <- function(L, M){
  return(0.5*L^(3/2)*sqrt(pi/M))
}

d1$er = er(d1$body_length, d1$body_mass)

#############
# heterozygosity
#############

calc_het <- function(h) {
  het = pull(h %>% filter(type == "s") %>%
               summarise(sum(var) / sum(denom)))
  denom = pull(h %>% filter(type == "s") %>%
               summarise(sum(denom)))
  ind = h[1, "individual"]
  return(c(ind, het, denom))
}

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

# rates
d1$speciation_rate = crates[ ldf[match(d1$lineage, ldf$lineage), "full_tree"] ]
d1$log_speciation_rate = log(d1$speciation_rate)

snake = c("Trilepida_koppesi", "Tantilla_melanocephala")
snake = getMRCA(t, snake)
snakes = t$tip.label[ phangorn::Descendants(t, snake, type = "tips")[[1]] ]
d1$clade = "lizard"
d1[d1$lineage %in% snakes, "clade"] = "snake"

# this range area is clearly in error
d1[which(d1$lineage == "Vanzosaura_savanicola"), "range_area"] = NA
lins = table(ll$lineage)
d1$inds = as.vector(lins[d1$lineage])

# for each explanatory variable
vars = c("pi", "body_mass", "range_area", "er", "inds")
xlabs = c(expression(pi), "log body mass", expression("log range area (" * km^2 * ")"),
          "log elongation", "log inds")
res1 = vector("list", length(vars))
plt1 = vector("list", length(vars))
res2 = vector("list", length(vars))
plt2 = vector("list", length(vars))

t2 = keep.tip(t, d1$lineage)
d2 = d1[match(t2$tip.label, d1$lineage), ]
xx = phytools::phylANOVA(t2, d2$clade, d2$log_speciation_rate)

for (i in 1:length(vars)) {
  # is it correlated with IBD?
  newvar = paste0("log_", vars[i])
  d1[ , newvar] = log(d1[, vars[i]])
  
  d2 = d1[complete.cases(d1[ ,c("log_slope", newvar)]),]
  
  t2 = keep.tip(t, d2$lineage)
  d3 = d2[match(t2$tip.label, d2$lineage), ]
  
  fmla <- as.formula(paste("log_slope", " ~ ", newvar))
  res1[[i]] = gls(fmla, correlation = corBrownian(phy = t2),
              data = d3, method = "ML")

  plt1[[i]] = ggplot(d3, aes_string(newvar, "log_slope")) + 
              geom_point(shape = 21, fill = "dodgerblue") + 
              xlab(xlabs[i]) + ylab(expression("log" ~ beta[IBD]))
  
  # is it correlated with speciation rate?
  fmla2 <- as.formula(paste("log_speciation_rate", " ~ ", newvar))
  res2[[i]] = gls(fmla2, 
              correlation = corBrownian(phy = t2),
              data = d3, method = "ML")
  
  if (i != 2) {
    plt2[[i]] = ggplot(d3, aes_string(newvar, "log_speciation_rate")) + 
      geom_point(shape = 21, fill = "dodgerblue") +
      xlab(xlabs[i]) + ylab("log speciation rate")
  } else {
    plt2[[i]] = ggplot(d3, aes_string(newvar, "log_speciation_rate")) + 
      geom_point(shape = 21, aes(fill = clade)) +
      xlab(xlabs[i]) + ylab("log speciation rate")
  }

}

gls1 = gls(log_speciation_rate ~ log_pi + log_slope, 
    correlation = corBrownian(phy = t2),
    data = d3, method = "ML")

summary(gls(log_slope ~ log_inds, 
    correlation = corBrownian(phy = t2),
    data = d3, method = "ML"))


abcd = (plt1[[1]] + plt1[[2]]) / (plt1[[3]] + plt1[[4]]) +
  plot_annotation(tag_levels = "A")
save_plot("~/Desktop/predictors_IBD.pdf", abcd, base_height = 6, base_width = 8)

save_plot("~/Desktop/bodymass_speciation.png", plt2[[2]], base_height = 3, base_width = 4)
save_plot("~/Desktop/pi_speciation.png", plt2[[1]], base_height = 3, base_width = 4)
