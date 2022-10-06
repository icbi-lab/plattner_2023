#!/usr/bin/env Rscript
#
# R script to generate Figure 2A from Plattner et al. 2022
#
# Author: Dietmar Rieder
# Date: 2022-08-25
#
#
library(ggplot2)
library(cowplot)
require(gridExtra)
library(reshape2)

library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(stringi)

library(conflicted)
conflict_prefer("get_legend", "cowplot")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")


##########
params = NULL
params$data_dir = "../../data/exomeseq"
params$result_dir = "../../results/exomeseq"

# create result dir
dir.create(params$result_dir, showWarnings = FALSE)

## Load data tables

# Exome coverage
coverage_data <- read_tsv(file.path(params$data_dir, "Exome_coverage_all.tsv"))
coverage_coding_data <- read_tsv(file.path(params$data_dir, "Exome_coverage_coding.tsv"))

sample_names <- unique(coverage_data$PatientID)

# Mutations
mutation_data <- read_tsv(file.path(params$data_dir, "SomaticVariant_summarytable_variantCounts.tsv"))
driver_data <- read_tsv(file.path(params$data_dir, "SomaticVariant_summarytable_codingVars_driver_genes.tsv"))


# immune gene data, load and make full data table
immuneGenes_data <- left_join(coverage_data,
                               read_tsv(file.path(params$data_dir, "ImmuneGenes_summarytable_codingVars_Immune_genes.tsv")),
                               by = "PatientID") |> 
  mutate(Gene = replace_na(Gene, "B2M")) |>
  select(-coverage)

# mismatch repair data, load and make full data table
mmr_data <- left_join(coverage_data,
                      read_tsv(file.path(params$data_dir, "MMR_summarytable_codingVars_MMR_genes.tsv")),
                      by = "PatientID") |> 
  mutate(Gene = replace_na(Gene, "MLH1")) |>
  select(-coverage)

# MSI
msi_data <- read_tsv(file.path(params$data_dir, "all_MSIscores.tsv"))

# NeoAntigens
neoantigen_data <- read_tsv(file.path(params$data_dir, "NeoAntigen_Counts.tsv"))

# CNVs
seqz_cnv_len_data <- read_tsv(file.path(params$data_dir, "sequenza_CNVs_len.tsv"))


## Gene lists

# load Intogen-DriverGene table with stats, first 27 (until CTNNB1) driver genes
intogen_driver_data <- read_tsv(file.path(params$data_dir, "IntOGen-DriverGenes_COREAD.tsv"), n_max = 27)

# load immune genes, filter to keep only those which we have data for
grasso_immuneGenes_data <- read_tsv(file.path(params$data_dir, "Grasso_et_al_2018_immuneGenes_order.tsv")) |>
  filter(SYMBOL %in% immuneGenes_data$Gene)


## Filter

# filter driver mutations and keep only genes from the first 27 Intogen CRC driver genes
driver_data <- driver_data |> filter(Gene %in% intogen_driver_data$Symbol)


## make variant data list

variant_data <- list(driver_data=driver_data, mmr_data=mmr_data, immuneGenes_data=immuneGenes_data)


## Clean up data

# clean up variants, keep only most severe consequences for transcripts
variant_data <- lapply(variant_data, function(vd) {
  vd$variant_type <- stri_replace_all_regex(vd$variant_type,
                                    pattern = c("splice_region_variant|\\&coding_sequence_variant",
                                                "\\&non_coding_transcript_variant",
                                                "\\&non_coding_transcript",
                                                "\\&NMD_transcript",
                                                "\\&",
                                                "inframe_insertion",
                                                "inframe_deletion",
                                                "_variant",
                                                "_"
                                                ),
                                    replacement = c('', 
                                                    '',
                                                    '',
                                                    '',
                                                    '',
                                                    'indel',
                                                    'indel',
                                                    '',
                                                    ' '),
                                    vectorize=FALSE)
  vd
})



## join tables and caclulate mutational load
mutational_load_data <- left_join(mutation_data, coverage_data, by="PatientID") |>
  rename(coverage_all = coverage) |>
  left_join(coverage_coding_data, by="PatientID") |>
  rename(coverage_coding = coverage) |>
  mutate(load_all_variants = all_variants / (coverage_all * 1e-06)) |>
  mutate(load_coding_variants = coding_variants / (coverage_coding * 1e-06)) |>
  left_join(msi_data, by="PatientID") |>
  mutate(microsatellite_status = ifelse(MSIscore > 5, "MSI", "MSS"))


# add mutational load to variant data 
variant_data <- lapply(variant_data, function(vd) {
  vd <- vd |>
    left_join(mutational_load_data |> 
                select(c("PatientID", "load_all_variants", "load_coding_variants")),
              by="PatientID")
  vd
})


# add mut load to neoantigen data
neoantigen_data <- neoantigen_data |>
  left_join(mutational_load_data |> 
              select(c("PatientID", "load_all_variants", "load_coding_variants")),
            by="PatientID")

# filter for genes that have mutations in our data
intogen_driver_data_filtered <- intogen_driver_data[intogen_driver_data$Symbol %in% unique(variant_data$driver_data$Gene),]


## make cnv data

# cnv data frame
cnv_len_data <- data.frame(PatientID=character(), cnv=double(), type=character(), stringsAsFactors=FALSE)

for (PatientID in sample_names) {
  if(PatientID %in% unique(seqz_cnv_len_data$PatientID)) {

    amp_len_only <- sum(seqz_cnv_len_data$length[seqz_cnv_len_data$PatientID == PatientID & seqz_cnv_len_data$nMajor > 1 & seqz_cnv_len_data$nMinor > 0 ]) / sum(seqz_cnv_len_data$length[seqz_cnv_len_data$PatientID == PatientID])
    del_len_only <- sum(seqz_cnv_len_data$length[seqz_cnv_len_data$PatientID == PatientID & seqz_cnv_len_data$nMinor < 1 & seqz_cnv_len_data$nMajor == 1]) / sum(seqz_cnv_len_data$length[seqz_cnv_len_data$PatientID == PatientID])
    amp_len_del <- sum(seqz_cnv_len_data$length[seqz_cnv_len_data$PatientID == PatientID & seqz_cnv_len_data$nMajor > 1 & seqz_cnv_len_data$nMinor == 0 ]) / sum(seqz_cnv_len_data$length[seqz_cnv_len_data$PatientID == PatientID])

    amp_len = amp_len_only + amp_len_del/2
    del_len = del_len_only + amp_len_del/2
    
    record <- list(PatientID, amp_len, "amp")
    cnv_len_data[nrow(cnv_len_data) + 1,] <- record  
    record <- list(PatientID, -del_len, "del")
    cnv_len_data[nrow(cnv_len_data) + 1,] <- record
  } else {

    record <- list(PatientID, NA, "amp")
    cnv_len_data[nrow(cnv_len_data) + 1,] <- record  
    record <- list(PatientID, NA, "del")
    cnv_len_data[nrow(cnv_len_data) + 1,] <- record
    
  }
  
}


# add mut load to CNV data
cnv_len_data <- cnv_len_data |>
  left_join(mutational_load_data |> 
              select(c("PatientID", "load_all_variants", "load_coding_variants")),
            by="PatientID")


### Start plot

# define MSI color palette
msi_Palette <- c("MSI" = "mediumpurple", "MSS"="chocolate1")

# define mutation color palette
mut_Palette <- c("frameshift" = "#56B4E9",
                 "indel" = "#F0E442",
                 "missense" = "#0072B2",
                 "splice acceptor" = "#D55E00",
                 "splice donor" = "#E69F00",
                 "splice donor exon" = "#009E73",
                 "stop gained" = "#FF423B",
                 na.value = "white")


## MIS status bars
msi_plot <- ggplot(mutational_load_data, aes(x=reorder(PatientID,-load_all_variants), y=1, fill=microsatellite_status, color=microsatellite_status, width = 1)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values=msi_Palette) +
  scale_colour_manual(values=msi_Palette) +
  scale_y_continuous(expand = c(0, 0, 0, 0.01),) +
  theme_void() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10))


## MMR gene mutation tiles
mmr_plot <- ggplot(variant_data$mmr_data, aes(x=reorder(PatientID,-load_all_variants), y=Gene, fill=variant_type, width=0.95, height=0.95)) +
  geom_tile() +
  scale_fill_manual(values=mut_Palette, drop = FALSE, na.value = 'white') +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ylab("MMR gene") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 9)) +
  labs(x = "")


## Mutational load data
ml_plot <- ggplot(mutational_load_data, aes(x=reorder(PatientID,-load_all_variants), y=load_all_variants)) + 
  geom_bar(stat = "identity", 
           color="black",
           fill="cyan3",
           position = "dodge") +
  scale_y_continuous(trans = 'log10', expand = c(0, 0, 0, 0.01),) +
  ylab("ML") +
  theme_bw() +
  theme(panel.spacing.x = unit(1, "mm"),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9))

## CNV data
cnv_plot <- ggplot(cnv_len_data, aes(x=reorder(PatientID,-load_all_variants), y=cnv, fill=type)) + 
  geom_bar(stat = "identity", color="black",
           position = position_stack()) +
  scale_fill_brewer(palette="Purples") +
  guides(fill = guide_legend(override.aes = list(colour = FALSE))) +
  ylab("CNV") +
  ylim(-1,1) +
  theme_bw() +
  theme(panel.spacing.x = unit(1, "mm"),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9))


## MHC Class I neoantigens: canonical
neoag_can_class_I_plot <- ggplot(neoantigen_data, aes(x=reorder(PatientID,-load_all_variants), y=(MHC_I_rank_filtered), fill="Class I")) + 
  geom_bar(stat = "identity",
           color="black",
           position = "dodge") +
  scale_fill_brewer(palette="Greens") +
  guides(fill = guide_legend(override.aes = list(colour = FALSE))) +
  scale_y_continuous(trans = 'log10', expand = c(0, 0, 0, 0.01),) +
  theme_bw() +
  theme(panel.spacing.x = unit(1, "mm"),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10))


## MHC Class I neoantigens: fusions
neoag_fus_class_I_plot <- ggplot(neoantigen_data, aes(x=reorder(PatientID,-load_all_variants), y=(MHC_I_neofuse_filtered), fill="Class I gene fusions")) + 
  geom_bar(stat = "identity",
           color="black",
           position = "dodge") +
  scale_fill_brewer(palette="Reds") +
  guides(fill = guide_legend(override.aes = list(colour = FALSE))) +
  scale_y_continuous(trans = 'log10', expand = c(0, 0, 0, 0.01),) +
  theme_bw() +
  theme(panel.spacing.x = unit(1, "mm"),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10))

## MHC Class II neoantigens: canonical
neoag_can_class_II_plot <- ggplot(neoantigen_data, aes(x=reorder(PatientID,-load_all_variants), y=(MHC_II_rank_filtered), fill="Class II")) + 
  geom_bar(stat = "identity",
           color="black",
           position = "dodge") +
  scale_fill_brewer(palette="Blues") +
  guides(fill = guide_legend(override.aes = list(colour = FALSE))) +
  scale_y_continuous(trans = 'log10', expand = c(0, 0, 0, 0.01),) +
  theme_bw() +
  theme(panel.spacing.x = unit(1, "mm"),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title=element_blank(),
        legend.text = element_text(size = 10))


## MHC Class II neoantigens: fusions
neoag_fus_class_II_plot <- ggplot(neoantigen_data, aes(x=reorder(PatientID,-load_all_variants), y=(MHC_II_neofuse_filtered), fill="Class II gene fusions")) + 
  geom_bar(stat = "identity",
           color="black",
           position = "dodge") +
  scale_fill_brewer(palette="Oranges") +
  guides(fill = guide_legend(override.aes = list(colour = FALSE))) +
  scale_y_continuous(trans = 'log10', expand = c(0, 0, 0, 0.01),) +
  theme_bw() +
  theme(panel.spacing.x = unit(1, "mm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 10, angle=45, vjust=0.1),
        legend.text = element_text(size = 10))

## CRC driver genes
driver_plot <- ggplot(variant_data$driver_data, aes(x=reorder(PatientID,-load_all_variants), y=Gene, fill=variant_type, width=0.95, height=0.95)) +
  geom_tile() +
  scale_fill_discrete(na.value = 'white') +
  scale_fill_manual(values=mut_Palette, guide = guide_legend(ncol = 3), limits = names(mut_Palette[1:7]), drop = FALSE) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(limits=rev(intogen_driver_data$Symbol[order(intogen_driver_data$Samples, decreasing = TRUE)]), expand = c(0, 0)) +
  ylab("CRC driver gene") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 9)) +
  labs(x = "")

## immune gens
immune_plot <- ggplot(variant_data$immuneGenes_data, aes(x=reorder(PatientID,-load_all_variants), y=Gene, fill=variant_type, width=0.95, height=0.95))

immune_plot <- immune_plot + geom_tile() +
  scale_fill_discrete(na.translate = F) +
  scale_fill_manual(values=mut_Palette,na.translate = F) +
  scale_y_discrete(limits=rev(grasso_immuneGenes_data$SYMBOL[order(grasso_immuneGenes_data$NR, decreasing = FALSE)]), expand = c(0, 0)) +
  ylab("immune related genes") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 9)) +
  labs(x = "")



## get legend from sub plots

legend <- plot_grid(
  plot_grid(
    get_legend(msi_plot), get_legend(cnv_plot), ncol = 2, align = 'vh', rel_heights=c(0.25, 0.25), rel_widths = c(0.25, 0.25)
  ),
  plot_grid(get_legend(driver_plot)),
  plot_grid(
    get_legend(neoag_can_class_I_plot), get_legend(neoag_fus_class_I_plot),
    get_legend(neoag_can_class_II_plot), get_legend(neoag_fus_class_II_plot),
    ncol = 2, nrow = 2, align = 'hv', axis = "b", rel_heights=c(0.125, 0.125, 0.125, 0.125), rel_widths = c(0.1, 0.3)
  ),
  align = 'vh', axis = 'br', ncol = 3, rel_heights = c(0.25, 0.25, 0.15), rel_widths = c(0.2, 0.55, 0.25)
)


msi_plot <- msi_plot + theme(legend.position='none')
mmr_plot <- mmr_plot + theme(legend.position='none')
ml_plot <- ml_plot + theme(legend.position='none')
cnv_plot <- cnv_plot + theme(legend.position='none')
neoag_can_class_I_plot <- neoag_can_class_I_plot + theme(legend.position='none')
neoag_fus_class_I_plot <- neoag_fus_class_I_plot + theme(legend.position='none')
neoag_can_class_II_plot <- neoag_can_class_II_plot + theme(legend.position='none')
neoag_fus_class_II_plot <- neoag_fus_class_II_plot + theme(legend.position='none')
driver_plot <- driver_plot + theme(legend.position='none')
immune_plot <- immune_plot + theme(legend.position='none')


# make final figure
final_figure <- ggdraw(
  plot_grid(
    plot_grid(
      msi_plot, 
      mmr_plot, 
      ml_plot, 
      cnv_plot, 
      driver_plot, 
      immune_plot, 
      neoag_can_class_I_plot, 
      neoag_fus_class_I_plot, 
      neoag_can_class_II_plot, 
      neoag_fus_class_II_plot, 
      ncol=1, align='v',
      rel_heights=c(0.02, 0.1, 0.075, 0.1, 0.35, 0.25, 0.05, 0.05, 0.05, 0.1)
      ),
    plot_grid(
      legend,
      ncol=1),
    rel_heights=c(0.40, 0.035),
    nrow = 2,
    axis = 'br', align = 'hv') +
    draw_label("neoantigen load", x=-0.002, y=0.2, vjust=1.9, angle=90, size = 10) + 
    theme(
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

# Save figure as png, svg and pdf
ggsave(
  file.path(params$result_dir, "Figure_2A.png"),
  plot = final_figure,
  width = 12.32,
  height = 13.27,
  dpi = 300
)
ggsave(
  file.path(params$result_dir, "Figure_2A.svg"),
  plot = final_figure,
  width = 8.85,
  height = 13.27,
  dpi = 300
)
ggsave(
  file.path(params$result_dir, "Figure_2A.pdf"),
  plot = final_figure,
  width = 8.85,
  height = 13.27,
  dpi = 300
)
