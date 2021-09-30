#!/usr/bin/env Rscript
# This was heavily modified 
library(leafcutter)
library(dplyr)
setwd('/prj/MAGE/analysis/baltica')

args <- commandArgs(trailingOnly=TRUE)
cat(args, "\n")
arguments <- list(
  #counts_file='leafcutter/DCM-vs-CTRL/DCM-vs-CTRL_perind_numers.counts.gz', 
  counts_file=args[1],
  #groups_file='leafcutter/DCM-vs-CTRL/diff_introns.txt', 
  groups_file=args[2], 
  #output_prefix="leafcutter_HCM_age_discrete/",
  output_prefix=args[3],
  #age_datatype="discrete"
  age_datatype=args[4],
  max_cluster_size=Inf,
  min_samples_per_intron=5, 
  min_samples_per_group=3,
  min_coverage=20,
  timeout=30,
  num_threads=30,
  exon_file='leafcutter/exons.gtf.gz',
  init='smart',
  seed=12345)

counts_file=arguments$counts_file
groups_file=arguments$groups_file

cat("Loading counts from",counts_file,"\n")
if (!file.exists(counts_file)) stop("File ",counts_file," does not exist")
counts=read.table(counts_file, header=T, check.names = F)
# rownames(counts) <- counts$chrom
# counts$chrom <- NULL
  
cat("Loading metadata from",groups_file,"\n")
if (!file.exists(groups_file)) stop("File ",groups_file," does not exist")
meta=read.table(groups_file, header=F, stringsAsFactors = F)
colnames(meta)[1:2]=c("sample","group")

meta2 <- read.csv('MAGE_metadata.txt', stringsAsFactors = F)
meta2 <- meta2[c('Run', 'race', 'Age', 'sex')]
colnames(meta2)=c('sample', 'race', 'age', 'sex')
meta = left_join(meta, meta2, by='sample')

group_names=unique(meta$group) # keep order from groups_file unless numeric
if (is.numeric(meta$group)) group_names=sort(group_names)
meta$group=factor(meta$group, group_names)

stopifnot(length(group_names)==2)

cat("Encoding as",group_names[1],"=0,",group_names[2],"=1\n")
numeric_x=as.numeric(meta$group)-1
if(arguments$age_datatype == 'discrete') {
  cat('Transforming age to factor with 4 bins')
  meta$age = cut(meta$age, 4)
}  
cat(arguments$age_datatype)

confounders=NULL
if (ncol(meta)>2) {
    confounders=meta[,3:ncol(meta),drop=F]
    # scale continuous confounders
    for (i in seq_len(ncol(confounders)))
        if (is.numeric(confounders[,i]))
            confounders[,i]=scale(confounders[,i])
    # convert factors to one-of-K encoding
    confounders=model.matrix( ~., data=confounders )
    confounders=confounders[,2:ncol(confounders),drop=F] # remove intercept
}

minimum_group_size=min(sum(numeric_x==0),sum(numeric_x==1))
if (minimum_group_size < arguments$min_samples_per_intron)
  stop("The number of samples in the smallest group is less than min_samples_per_intron, which means no clusters are testable. You can reduce min_samples_per_intron using the -i argumentsion, but note that we have only carefully checked the calibration of leafcutter p-values down to n=4 samples per group.")
if (minimum_group_size < arguments$min_samples_per_group)
  stop("The number of samples in the smallest group is less than min_samples_per_group, which means no clusters are testable. You can reduce min_samples_per_intron using the -g argumentsion, but note that we have only carefully checked the calibration of leafcutter p-values down to n=4 samples per group.")

require(doMC)
registerDoMC(arguments$num_threads)

cat("Settings:\n")
print(arguments)

cat("Running differential splicing analysis...\n")
results <- differential_splicing(counts, numeric_x, confounders=confounders, max_cluster_size=arguments$max_cluster_size, min_samples_per_intron=arguments$min_samples_per_intron, min_samples_per_group=arguments$min_samples_per_group, min_coverage=arguments$min_coverage, timeout=arguments$timeout, init=arguments$init, seed=arguments$seed ) 

cat("Saving results...\n")

# Make cluster table
cluster_table          = cluster_results_table(results)
cluster_table$cluster  = add_chr(cluster_table$cluster)

# Add gene names to clusters if an exon file is available
if (!is.null(arguments$exon_file)) {
  cat("Loading exons from",arguments$exon_file,"\n")
  if (file.exists(arguments$exon_file)) {
     tryCatch( {
          exons_table     = read.table(arguments$exon_file, header=T, stringsAsFactors = F)
          intron_meta     = get_intron_meta(rownames(counts))
          exons_table$chr = add_chr(exons_table$chr)
          intron_meta$chr = add_chr(intron_meta$chr)
          clu_gene_map    = map_clusters_to_genes(intron_meta, exons_table)
          cluster_table   = merge(cluster_table, clu_gene_map, by.x="cluster", by.y="clu", all.x=TRUE)
     }, error=function(err) warning(as.character(err)) ) # ignore errors here
  } else warning("File ",arguments$exon_file," does not exist")
} else cat("No exon_file provided.\n")
write.table( cluster_table, paste0(arguments$output_prefix,"_cluster_significance.txt"), quote=F, sep="\t", row.names = F)

# Write effect size table
effect_size_table                = leaf_cutter_effect_sizes(results)
colnames(effect_size_table)[3:4] = group_names
effect_size_table$intron = add_chr(effect_size_table$intron)
write.table(effect_size_table, paste0(arguments$output_prefix,"_effect_sizes.txt"), quote=F, col.names = T, row.names = F, sep="\t")

# Save RData image
save.image(paste0(arguments$output_prefix,"_cluster_significance.RData"));
cat("All done, exiting\n")
