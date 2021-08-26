import glob
import re
#import pandas as pd

# GENES = glob.glob1("data/cis-eQTL-sumStats/byGene", "*.stats.txt")
# GENES = [re.sub("\.stats.txt", "", s) for s in GENES]

CHROMS = [i for i in range(1,23)]

# N_GENE_PER_GROUP=300
# Groups = ["group_" + str(i // N_GENE_PER_GROUP + 1) for i in range(len(GENES))]
# GroupLookup = dict(zip(GENES, Groups))

FOLDERS = glob.glob1("data/cis-eQTL-sumStats/byGene", "group*")


# def get_group_by_gene(wildcards):
#     '''Use gene name to lookup pre-assigned groupname, then pass to snakemake for resource allocation'''
#     return GroupLookup.get(wildcards.gene)
    

def get_group(wildcards):
    return wildcards.group

#####################################################################
#                  rules, rules, rules                              #
#####################################################################



localrules: all, getconda


rule all:
    input:
        expand("data/cis-eQTL-sumStats/byGene/{group}/{group}.done", group=FOLDERS),
        expand("output/byGene/{group}/{group}.done", group=FOLDERS)
#         expand("data/cis-eQTL-sumStats/byGene/{gene}.sumstats.gz", gene=GENES),
#         expand("output/byGene/{gene}_baseline.results", gene=GENES)


        
#-----------------------------------------------------------------------------
        
# run this rule once before submitting anything to clusters
# used to download and install required environments
rule getconda:
    input: "env/environment.yml"
    output: touch("getconda.done")
    conda:
        "env/environment.yml"
    shell:
        "echo done"


#-----------------------------------------------------------------------------



# split full cis-eQTL summary stats into 22 autosome subsets; helpful for parallelizatino
rule split_sumStats:
    input: "data/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
    output:
        done = touch("data/cis-eQTL-sumStats/cis-eQTL-sumStats.done")
    params:
        runmode = "run",
        out_prefix = "data/cis-eQTL-sumStats/cis-eQTL-sumStats."
    threads: 8
    resources:
        cpu=8,
        mem_mb=30000,
        time=120
    script:
        "scripts/SplitSummaryStatsbyChr.R"


# further split summar stats, each file contains all SNPs for the eQTLs for a particular gene
# this rule is created primarily for documentation purposes
rule make_perGene_sumStats:
    input:
        perChrStats = "data/cis-eQTL-sumStats/cis-eQTL-sumStats.1.txt" # the script doesn't need inputs from snakemake
    output:
        done = touch("data/cis-eQTL-sumStats/byGene/perGene.done")
    threads: 6
    script:
        "scripts/MakeperGeneStats_v2.R"

#-----------------------------------------------------------------------------


# Use ldsc's munge_sumstats to produce required summary statistics file
# NOTE, this is run on a per-gene bases
# for clusters that are slow to allocate resource, split all genes into a reasonable number of groups
rule munge_sumstats:
    input:
        snplist = "data/w_hm3.snplist"
    output: 
        done = touch("data/cis-eQTL-sumStats/byGene/{group}/{group}.done")
    group: get_group
    params:
        out_prefix_folder = "data/cis-eQTL-sumStats/byGene/{group}",
        munge = "../tools/ldsc/munge_sumstats.py"
    conda:
        "env/environment.yml"
    threads: 1
    resources:
        cpu=1,
        mem_mb=10000,
        time=2000 # min
    shell:
        '''
        if [[ -e {output.done} ]]; then
            truncate -s 0 {output.done}
        else
            echo
        fi 
        
        folder=data/cis-eQTL-sumStats/byGene/{wildcards.group}
        inputfiles=(${{folder}}/*.stats.txt)
        
        echo $(date)
        for input in ${{inputfiles[@]}}; do
            gene_name=$(basename -s .stats.txt $input)
            prefix="{params.out_prefix_folder}/$gene_name"
            echo python {params.munge} --sumstats $input \
                --merge-alleles {input.snplist} \
                --out $prefix \
                --a1-inc \
                --N-col NrSamples \
                --a1 AssessedAllele \
                --a2 OtherAllele \
                --signed-sumstats Zscore,0 >> {output.done}
            python {params.munge} --sumstats $input \
                --merge-alleles {input.snplist} \
                --out $prefix \
                --a1-inc \
                --N-col NrSamples \
                --a1 AssessedAllele \
                --a2 OtherAllele \
                --signed-sumstats Zscore,0
        done
        
        echo $(date)
        '''


# run ldsc partition heritability
# run on per gene basis
rule h2_baseline:
    input:
        sumstats = rules.munge_sumstats.output.done
    output: 
        done = "output/byGene/{group}/{group}.done"
    group: get_group
    params:
        baseline_prefix = "data/baseline/baseline.",
        weights_prefix = "data/weights_hm3_no_hla/weights.",
        freq_prefix = "data/1000G_frq/1000G.mac5eur.",
        out_prefix_folder = "output/byGene/{group}",
        ldsc = "../tools/ldsc/ldsc.py"
    conda:
        "env/environment.yml"
    threads: 1
    resources:
        cpu=1,
        mem_mb=10000,
        time=2000 # min
    shell:
        '''
        
        if [[ -e {output.done} ]]; then
            truncate -s 0 {output.done}
        else
            echo
        fi 
        
        folder=data/cis-eQTL-sumStats/byGene/{wildcards.group}
        inputfiles=(${{folder}}/*.sumstats.gz)
        
        
        echo $(date)
        
        for input in ${{inputfiles[@]}}; do
            gene_name=$(basename -s .sumstats.gz $input)
            prefix="{params.out_prefix_folder}/${{gene_name}}_baseline"
            echo python {params.ldsc} \
                --h2 $input \
                --w-ld-chr {params.weights_prefix} \
                --ref-ld-chr {params.baseline_prefix} \
                --overlap-annot \
                --frqfile-chr {params.freq_prefix} \
                --out $prefix \
                --print-coefficients \
                --no-intercept >> {output.done}
            python {params.ldsc} \
                --h2 $input \
                --w-ld-chr {params.weights_prefix} \
                --ref-ld-chr {params.baseline_prefix} \
                --overlap-annot \
                --frqfile-chr {params.freq_prefix} \
                --out $prefix \
                --print-coefficients \
                --no-intercept
        done
        echo $(date)

        '''

# compute total heritability per gene. This should be run manually in an interactive env, preferably with lots of cores
rule hsq:
    input: "output/byGene/group1/group1.done" # the script doesn't need inputs from snakemake
    output: touch("output/hsq.done")
    threads: 30
    script:
        "scripts/hsq.R"
        
# get summary of eQTLs from GTEx eQTL analyses
rule GTEx:
    input: "output/byGene/group1/group1.done" # the script doesn't need inputs from snakemake
    output: touch("output/GTEx_summary.done")
    threads: 30
    script:
        "scripts/sumGTExeQTLs.R"