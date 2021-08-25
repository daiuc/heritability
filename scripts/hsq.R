library(parallel)
library(data.table)
library(tidyverse)

if (exists("snakemake@threads[[1]]")) {
    threads = snakemake@threads[[1]]
} else {
    threads = 30
}
NCores = min(threads, detectCores())


#setwd("/project2/xuanyao/chao/heritability")

out.file.name = "output/hsq.txt"

# group folders
groups = dir("output/byGene", "group[0-9]{1,2}$", full.names = T)
#names(groups) = dir("output/byGene", "group[0-9]{1,2}$", full.names = F)

# file names for each group, name each element with gene name
files_in_groups = map(groups, ~dir(.x, ".+results$"))
all_result_files = map2(groups, files_in_groups, function(g,f) {map_chr(f, ~paste(g, "/", .x, sep = ""))} ) %>% unlist 
names(all_result_files) = str_split(all_result_files, "/", simplify = F) %>% 
                            map_chr(~.x[4]) %>% 
                            map_chr(~str_replace(.x, "_baseline.results", ""))


# function to read *.results files
read_result = function(fn) {
    gene = str_split(fn, "/", simplify = F) %>% 
                map_chr(~.x[4]) %>% 
                map_chr(~str_replace(.x, "_baseline.results", ""))
    rs = read_table(fn) %>% 
            select(Category, M_annot, M_tot, Coefficient) %>% 
            arrange(Category) %>%
            add_column("gene" = gene)
    return(rs)
}

# parallelize reading all *.result files, each list element is a df
res.ls = mclapply(all_result_files, read_result, mc.cores = NCores)


# compute tau_c_bar: verage per SNP heritability
res = bind_rows(res.ls)
tau_c_bar = group_by(res, Category) %>% 
                summarise(tau_c_bar = mean(Coefficient)) %>% 
                arrange(Category) %>% deframe

# [total heritability per gene] = [tau_c_bar] * [number of snps per category]
h2.tot = mclapply(res.ls, function(df) {mutate(df, h2 = M_annot * tau_c_bar) %>% pull(h2) %>% sum}, 
                  mc.cores = NCores)
h2.tot = enframe(unlist(h2.tot), name = "gene", value = "h2")

# get number of SNPs for each gene
m_tot = map_dbl(res.ls, ~.x$M_tot[1])

# add m_tot to dataframe
h2.tot = add_column(h2.tot, m_tot)
# output per gene heritability
write_tsv(x = h2.tot, file = out.file.name)