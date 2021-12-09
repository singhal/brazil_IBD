library(dplyr)
library(ggplot2)
library(cowplot)

ll = read.csv("~/Dropbox (Personal)/brazil_gene_flow/data/brazil_samples_v8.csv",
              stringsAsFactors = F)
ll2 =  ll[complete.cases(ll$LAT), ]
lins = table(ll2$lineage)
lins = names(lins[ lins > 3])
lins = lins[lins != "Phrynonax_poecilonotus"]

dd = read.csv("~/Dropbox (Personal)/brazil_gene_flow/results/IBD-2020-12-02.csv",
              stringsAsFactors = F)
dd = dd[dd$div_type == 'inv_fst', ]
dd$sig = ifelse(dd$sig < 0.05, TRUE, FALSE)
dd[dd$lineage == 'Vanzosaura_rubricauda_2', "lineage"] = "Vanzosaura_savanicola"

div = read.csv("~/Dropbox (Personal)//brazil_gene_flow/results/divergence-2020-12-02.csv", stringsAsFactors = F)
div[div$lineage == 'Vanzosaura_rubricauda_2', "lineage"] = "Vanzosaura_savanicola"
div2 = div[div$lineage %in% lins, ]
div2 = div2[is.finite(div2$ln_geo_dist), ]
div2 = left_join(div2, dd %>% dplyr::select(lineage, sig))
div2$lineage = gsub("_", " ", div2$lineage)

lins2 = unique(div2$lineage)
linsxx = split(lins2, ceiling(seq_along(lins2)/15))

for (i in 1:length(linsxx)) {
xx = ggplot(div2 %>% filter(lineage %in% linsxx[[i]]), 
            aes(ln_geo_dist, inv_fst)) + 
  geom_point(color = "gray70") +
  geom_smooth(method = "lm", fill = NA, aes(col = sig),
              lwd = 0.7) +
  facet_wrap(~lineage, nrow = 5, ncol = 3) +
  ylim(-0.5, 5) + xlim(-2, 10) + 
  theme_classic() +
  theme(strip.text = element_text(face = "italic")) +
  xlab("log(geo dist, km)") +
  ylab(expression(F[ST] / (1 - F[ST]) )) +
  scale_colour_manual(values = c( "black" , "red")) +
  theme(legend.position = "none")
save_plot(paste0("~/Desktop/ibd_lineages", i, ".png"), 
          xx, base_height = 10, base_width = 8)
}
