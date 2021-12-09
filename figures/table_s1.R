# list of all samples, species identity, lat & long, SRA ID, basic quality info

ll = read.csv("~/Dropbox (Personal)/brazil_gene_flow/data/brazil_samples_v7.csv",
              stringsAsFactors = F)
ll2 =  ll[complete.cases(ll$LAT), ]
lins = table(ll2$lineage)
lins = names(lins[ lins > 3])
lins = lins[lins != "Phrynonax_poecilonotus"]
ll3 = ll2[ll2$lineage %in% lins, ]

ll4 = ll3 %>% dplyr::select(sample, species_checked, lineage,
               LAT, LON)
q = read.csv("~/Dropbox (Personal)//brazil_gene_flow/data/metadata/site_counts.csv",
             stringsAsFactors = F)
ll5 = ll4 %>% left_join(q, by = c("sample" = "individual"))
ll5$tot_coverage = round(ll5$tot_coverage / ll5$num_sites, 1)
ll5$num_sites = round(ll5$num_sites / 1e6, 1)
names(ll5) = c("individual", "nominal species",
               "revised taxon name", "latitude",
               "longitude", "# of loci", 
               "# of sites (Mb)", "avg. coverage")
ll5$SRA = "TBD"
ll5$`nominal species` = gsub("_", " ", ll5$`nominal species`)
ll5$`revised taxon name` = gsub("_", " ", ll5$`revised taxon name`)
write.csv(ll5, "~/Desktop/Table_S1.csv", quote = F, row.names = F)
