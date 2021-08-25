library(parallel)
library(data.table)
library(dplyr)
library(purrr)


NCores = min(22, detectCores())

# paralel read each chromosome
files = dir("data/cis-eQTL-sumStats/", "^cis.+\\.txt", full.names = T)
cis = mclapply(files, function(x) fread(x), mc.cores = NCores) 

# choose 5 genes from each chromosome
set.seed(22)
choose_genes = mclapply(cis, function(df) pull(df, GeneSymbol) %>% unique, mc.cores = NCores)

# parallel write to file, file named: GENE.stats.txt                        
mcmapply(function(x,y){
    walk(y, ~ filter(x, GeneSymbol == .x) %>% 
                    fwrite(file = paste("data/cis-eQTL-sumStats/byGene/", .x, ".stats.txt", sep=""), sep = "\t"))
    return("done")
}, x=cis, y=choose_genes, mc.cores = NCores)