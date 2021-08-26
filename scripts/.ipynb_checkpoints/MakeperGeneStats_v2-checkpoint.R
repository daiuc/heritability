library(furrr)
library(data.table)
library(dplyr)
library(purrr)
library(tictoc)



setwd('/project2/xuanyao/chao/heritability')
NCores = snakemake@threads[[1]]-2
print(paste("# Using", NCores, "cores."))
plan(multisession, workers = NCores)
future_globals_maxSize = 40*1024*1024^2 # this is effectively 30GB
options(future.globals.maxSize = future_globals_maxSize)

tic()
# paralel read each chromosome
files = dir("data/cis-eQTL-sumStats", "^cis.+\\.txt", full.names = T)
prefix = '/project2/xuanyao/chao/heritability/'
files = paste(prefix, files, sep="")

print("# Reading chromosome level sumstats")
cis = future_map(files, function(x) fread(x))         


# get all the gene names            
all.genes = future_map(cis, function(df) pull(df, GeneSymbol) %>% unique)
all.genes.v = unlist(all.genes)


# devide genes into 70 groups
genes.per.group = ceiling(length(all.genes.v) / 70)
groups = 1:length(all.genes.v) %/% genes.per.group + 1
groups = paste("group", groups, sep="")


print(paste("# Assigning", length(all.genes.v), "genes to", length(groups), "groups. On average", genes.per.group, "genes per group."))
# construct gene name, gene group lookup table dataframe
genes.groups.lookup = data.frame(gene = all.genes.v, gp = groups, stringsAsFactors = F)

# for each of the 22 chromosome devided dataset, create corresponding gene,genegroup lookup
#choose_genes = mclapply(cis, function(df) pull(df, GeneSymbol) %>% unique, mc.cores = NCores)
choose_genes = map(all.genes, ~ data.frame("gene" = .x, stringsAsFactors = F))
choose_genes = map(choose_genes, ~ left_join(.x, genes.groups.lookup, by = c("gene" = "gene")))

                        
# check write file path exists
for (g in unique(genes.groups.lookup$gp)) {
    p = paste("data/cis-eQTL-sumStats/byGene/", g, sep = "")
    if (!dir.exists(p)) {
        dir.create(path = p)
    }
}

# function to write file, expenct inputs are elements of choose_gene and cis list
write_to_file = function(df.cis, df.gene.lookup) {
    filename = map2_chr(df.gene.lookup$gene, df.gene.lookup$gp, 
                        ~ paste(prefix, "data/cis-eQTL-sumStats/byGene/", .y, "/", .x, ".stats.txt", sep=""))
    walk2(df.gene.lookup$gene, filename, 
          ~ filter(df.cis, GeneSymbol == .x) %>% 
             fwrite(file = .y, sep = "\t"))
}
                        
# write to file, file named: GENE.stats.txt  
# disabled parallel because RCC keeps sending kill signal.
#future_walk2(cis[1:2], choose_genes[1:2], write_to_file)
print("# Now writing sumstats for each gene to file. Using single core.")
walk2(cis, choose_genes, write_to_file)
print("# Success. Done!")
toc()