library(parallel)
library(data.table)
library(dplyr)
library(purrr)



NCores = min(22, detectCores())

# paralel read each chromosome
files = dir("data/cis-eQTL-sumStats/", "^cis.+\\.txt", full.names = T)
cis = mclapply(files, function(x) fread(x), mc.cores = NCores)
               

               
# get all the gene names

                        
all.genes = mclapply(cis, function(df) pull(df, GeneSymbol) %>% unique, mc.cores = NCores)
all.genes = unlist(all.genes)


# devide genes into 70 groups
genes.per.group = ceiling(length(all.genes) / 70)
groups = 1:length(all.genes) %/% genes.per.group + 1
groups = paste("group", groups, sep="")

# construct gene name, gene group lookup table dataframe
genes.groups.lookup = data.frame(gene = all.genes, gp = groups, stringsAsFactors = F)

# for each of the 22 chromosome devided dataset, create corresponding gene,genegroup lookup
choose_genes = mclapply(cis, function(df) pull(df, GeneSymbol) %>% unique, mc.cores = NCores)
choose_genes = map(choose_genes, ~ data.frame("gene" = .x, stringsAsFactors = F))
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
                        ~ paste("data/cis-eQTL-sumStats/byGene/", .y, "/", .x, ".stats.txt", sep=""))
    walk2(df.gene.lookup$gene, filename, 
          ~ filter(df.cis, GeneSymbol == .x) %>% 
             fwrite(file = .y, sep = "\t"))
}
                        
# parallel write to file, file named: GENE.stats.txt                        
mcmapply(write_to_file, df.cis=cis, df.gene.lookup=choose_genes, mc.cores = NCores)