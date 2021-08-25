library(parallel)
library(data.table)
library(tidyverse)

if (exists("snakemake@threads[[1]]")) {
    threads = snakemake@threads[[1]]
} else {
    threads = 30
}
NCores = min(threads, detectCores())

#setwd("/home/chaodai/chao/heritability")

out.file.name = "output/GTEx_v8_sigif_eQTL_summ.txt"




SEL_COLS = c("phenotype_id", "variant_id", "slope", "maf")
eQTL.files = map_chr(dir("data/GTEx_Analysis_v8_eQTL_EUR", "signif_pairs.txt.gz"), ~paste("data/GTEx_Analysis_v8_eQTL_EUR/", .x, sep = ""))
eQTLs = mclapply(eQTL.files, fread, sep = "\t", header = T, select = SEL_COLS, mc.cores = NCores)
eQTLs = mclapply(eQTLs, function(df) mutate(df, gene_id = str_extract(phenotype_id, "ENSG[0-9]+")) %>% 
                            select(-phenotype_id) %>%
                            filter(maf > 0.05 & maf < 0.95))


### Compute total number of qualifying significant unique eQTLs across tissue samples ###

# for each dataframe in eQTL, select only gene_id and variant_id
# group by gene_id and nest to create a data column 
# convert each element of the data column to be a character vector
# then ungroup
# bind_rows for all processed df in eQTL
# use below defined function that takes in a list, with each element being a chr vector, report either a vector if desired or vector length
CombineChr = function(list_of_chr) {
    combined_chr = reduce(list_of_chr, union)
    return(combined_chr)
}

CombineChr_len = function(list_of_chr) {
    combined_chr = reduce(list_of_chr, union)
    return(length(combined_chr))
}


eQTL.variants = mclapply(eQTLs, function(df) {select(df, gene_id, variant_id) %>%
                             group_by(gene_id) %>% 
                             nest %>%
                             mutate_at("data", function(data.col) map(data.col, ~ unlist(.x, use.names = F)))
                                       },mc.cores = NCores)
eQTL.variants = bind_rows(eQTL.variants) %>%
                            group_by(gene_id) %>%
                            summarise(N_QTLs = CombineChr_len(data))


# for each QTL dataframe, do filtering and summarization
Sum_QTLs = function(df) {
    df = filter(df, maf > 0.05 & maf < 0.95) %>%
            group_by(gene_id) %>%
            summarise(max.abs.slope = max(abs(slope))) %>%
            ungroup
    return(df)
}
                                              
eQTLs = mclapply(eQTLs, Sum_QTLs, mc.cores = NCores) %>%
            bind_rows(.) %>% 
            group_by(gene_id) %>% 
            summarise(max.abs.slope = max(max.abs.slope)) %>% 
            ungroup
                                              
# now join max slope and N_eQTLs together
eQTLs = inner_join(eQTLs, eQTL.variants, by = "gene_id")
                                              
# export file
write_tsv(x = eQTLs, file = out.file.name)