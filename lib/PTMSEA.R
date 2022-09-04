#!/usr/local/bioinf/bin/Rscript

#$ -S /usr/local/bioinf/bin/Rscript
#$ -pe smp 1
#$ -cwd
#$ -V
#$ -N PTMSEA
#$ -e ./$JOB_NAME.err
#$ -o ./$JOB_NAME.log
#$ -r y

# specify command line arguments
library("optparse")

option_list <- list(
  make_option( c("-r", "--script"), action='store', type='character',  dest='script', help='Path to R script'),
  make_option( c("-i", "--input"), action='store', type='character',  dest='input_ds', help='Path to input GCT file.'),
  make_option( c("-o", "--ouptut"), action='store', type='character',  dest='output_prefix', help='File prefix for output files.', default='out'),
  make_option( c("-d", "--db"), action='store', type='character',  dest='gene_set_databases', help='Path to gene set database (GMT format).'),
  make_option( c("-n", "--norm"), action='store', type='character',  dest='sample_norm_type', help='Sample normalization: "rank", "log", "log.rank" or "none".', default = 'rank'),
  make_option( c("-w", "--weight"), action='store', type='numeric',  dest='weight', help='When weight==0, all genes have the same weight; if weight>0 actual values matter and can change the resulting score.', default = 0.75),
  make_option( c("-c", "--correl"), action='store', type='character',  dest='correl_type', help='Correlation type: "rank", "z.score", "symm.rank".', default = 'z.score'),
  make_option( c("-t", "--test"), action='store', type='character',  dest='statistic', help='Test statistic: "area.under.RES", "Kolmogorov-Smirnov"', default = 'area.under.RES'),
  make_option( c("-s", "--score"), action='store', type='character',  dest='output_score_type', help='Score type: "ES" - enrichment score,  "NES" - normalized ES', default = 'NES'),
  make_option( c("-p", "--perm"), action='store', type='numeric',  dest='nperm', help='Number of permutations', default = 1000),
  make_option( c("-m", "--minoverlap"), action='store', type='numeric',  dest='min_overlap', help='Minimal overlap between signature and data set.', default = 10),
  make_option( c("-x", "--extendedoutput"), action='store', type='character',  dest='extended_output', help='If TRUE additional stats on signature coverage etc. will be included as row annotations in the GCT results files.', default = TRUE),
  make_option( c("-e", "--export"), action='store', type='character',  dest='export_signat_gct', help='For each signature export expression GCT files.', default = TRUE),
  make_option( c("-g", "--globalfdr"), action='store', type='character',  dest='global_fdr', help='If TRUE global FDR across all data columns is calculated.', default = FALSE),
  make_option( c("-l", "--lightspeed"), action='store', type='character',  dest='multi_core', help='If TRUE processing will be parallized across gene sets. (I ran out of single letters to define parameters...)', default = TRUE),
  make_option( c("-y", "--sparecores"), action='store', type='numeric',  dest='sparecores', help='Num. of cores to keep free', default = 1)
)

## #####################################

# parse command line parameters
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);
#opt <- parse_param(option_list) 

# source the actual script
script.dir <- opt$script
source(file.path(script.dir, "src", "ssGSEA2.0.R"))

# hard-coded parameters
log.file <- paste(opt$output.prefix, '_ssgsea.log.txt', sep='')


## ######################################################################################################
##
##                   run ssGSEA
##
## ######################################################################################################

set.seed(123)
res <- ssGSEA2(
  input.ds=opt$input_ds,
  output.prefix=opt$output_prefix,
  gene.set.databases=opt$gene_set_databases,
  sample.norm.type=opt$sample_norm_type,
  weight=opt$weight,
  statistic=opt$statistic,
  output.score.type=opt$output_score_type,
  nperm=opt$nperm,
  min.overlap=opt$min_overlap,
  correl.type=opt$correl_type,
  export.signat.gct=opt$export_signat_gct,
  extended.output=opt$extended_output,
  global.fdr=opt$global_fdr,
  par=opt$multi_core,
  spare.cores=opt$sparecores,
  log.file=log.file
)














