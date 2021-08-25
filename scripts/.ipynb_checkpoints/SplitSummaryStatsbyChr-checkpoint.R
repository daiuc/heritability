library(purrr)
library(data.table)

debug_mode = snakemake@params[[1]]
chroms = 1:22
if (debug_mode == "run") {
    cis = fread(snakemake@input[[1]])
    file.names = paste(snakemake@params[[2]], chroms, ".txt", sep = "")
} else {
    cis = fread("data/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt", nrows = 100000)
    file.names = paste("data/cis-eQTL-sumStats/cis-eQTL-sumStats.", chroms, ".txt", sep = "")
}
dts = map(chroms, ~ cis[SNPChr == chroms[.x]][order(SNPPos)])
if (!dir.exists("data/cis-eQTL-sumStats")) dir.create("data/cis-eQTL-sumStats")
walk2(dts, file.names, ~ fwrite(.x, file = .y, sep = "\t", quote = F))