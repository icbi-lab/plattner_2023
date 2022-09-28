


### PHOSPHOPROTEOMICS ----------------------------------------------------------


### mapDIA ----


makeDIAcontrasts <- function(design, formula, ref = NULL){
  
  # Make contrast matrix as used by mapDIA
  # e.g. ~ organoid|treatment to compare treatment levels in each organoid subgroup
  
  form_terms <- strsplit(x = labels(terms(formula)), split = "|", fixed = TRUE)[[1]]
  form_terms <- trimws(form_terms)
  
  if (length(form_terms) > 1) split_by <- form_terms[1] else split_by <- NULL
  if (length(form_terms) > 1) compare_class <- form_terms[2] else compare_class <- form_terms[1]
  
  df <- data.frame(design[,c(split_by, compare_class), drop = FALSE])
  if (!is.null(split_by)) df$group <- paste0(df[[compare_class]], "_", df[[split_by]]) else df$group <- df[[compare_class]]
  df$size <- sapply(df$group, function(tmp) sum(tmp == df$group))
  
  
  orig_check <- df$group
  
  
  cdf <- df[!duplicated(df$group),]
  rownames(cdf) <- cdf$group
  groups <- cdf$group
  contr <- data.frame(matrix(data = 0, nrow = length(groups), ncol = length(groups), dimnames = list(groups, groups)))
  diag(contr) <- "-"
  
  if (is.null(ref)) ref <- levels(df[[compare_class]])[1]
  if (is.null(ref)) ref <- df[[compare_class]][1]
  
  controls <- unique(df$group[df[[compare_class]] == ref])
  
  if (!is.null(split_by)){
    for (con in controls){
      tmpclass <- unique(df[[split_by]][df$group == con])
      comp_samples <- unique(df$group[df[[split_by]] == tmpclass & df$group != con])
      contr[comp_samples,con] <- 1
    }
  } else {
    for (con in controls){
      comp_samples <- unique(df$group[df$group != con])
      contr[comp_samples,con] <- 1
    }
  }
  
  if( !all(df$group == rep(cdf$group, cdf$size)) ) stop("Error: All replicates must be grouped!")
  
  return(list(samples = rownames(df), design = cdf, contrasts = contr))
}




setmapDIAparams <- function(design, formula, input_data, ref = "control",
                            params_file = "mapDIA_params.txt", normalization = "tis",
                            SDF= 2, MIN_CORREL= 0.25, MIN_OBS = 0, MIN_FRAG_PER_PEP = 3, MAX_FRAG_PER_PEP = 5, MIN_PEP_PER_PROT = 1,
                            MIN_DE = 0.01, MAX_DE = 0.99, MAX_PEP_PER_PROT = 10){
  
  
  ### DESIGN AND CONTRASTS
  
  # check sample names
  header <- read.delim2(input_data, dec = ".", nrows = 5)
  anno_cols <- c("ProteinName", "PeptideSequence", "FragmentIon", "RT")
  samples <- colnames(header)[!colnames(header) %in% anno_cols]
  
  if ("RT" %in% colnames(header)){ stopifnot(which(colnames(header) == "RT") == ncol(header)) } # RT must be last
  
  stopifnot(all(rownames(design) %in% samples))
  design <- design[samples,]
  
  mapdia <- makeDIAcontrasts(design = design, formula = formula, ref = ref)
  
  stopifnot(all(mapdia$samples == samples))
  stopifnot(sum(mapdia$design$size) == length(samples))
  
  
  ### PARAMS FILE
  
  # Input file
  cat("###input file", "\n", file = params_file)
  cat(paste0("FILE=", basename(input_data)), "\n\n", file = params_file, append = TRUE)
  
  # Module
  cat("### MODULE data through MRF model", "\n", file = params_file, append = TRUE)
  cat("MODULE = none", "\n\n", file = params_file, append = TRUE)
  
  # Experimental design
  cat("### Experimental design", "\n", file = params_file, append = TRUE)
  cat("EXPERIMENTAL_DESIGN=IndependentDesign", "\n\n", file = params_file, append = TRUE)
  
  # Normalization
  cat("### Normalization", "\n", file = params_file, append = TRUE)
  cat(paste0("NORMALIZATION= ", normalization), "\n\n", file = params_file, append = TRUE)
  
  # Filter
  cat("### Filter", "\n", file = params_file, append = TRUE)
  cat(paste0("SDF= ", SDF), "\n", file = params_file, append = TRUE)
  cat(paste0("MIN_CORREL= ", MIN_CORREL), "\n", file = params_file, append = TRUE)
  cat(paste0("MIN_OBS= ", paste(rep(MIN_OBS, ncol(mapdia$contrasts)), collapse = " ")), "\n", file = params_file, append = TRUE)
  cat(paste0("MIN_FRAG_PER_PEP= ", MIN_FRAG_PER_PEP), "\n", file = params_file, append = TRUE)
  cat(paste0("MAX_FRAG_PER_PEP= ", MAX_FRAG_PER_PEP), "\n", file = params_file, append = TRUE)
  cat(paste0("MIN_PEP_PER_PROT= ", MIN_PEP_PER_PROT), "\n\n", file = params_file, append = TRUE)
  
  # Sample information
  cat("### Sample information", "\n", file = params_file, append = TRUE)
  cat(paste0("LABELS= ", paste(mapdia$design$group, collapse = " ")), "\n", file = params_file, append = TRUE)
  cat(paste0("SIZE= ", paste(mapdia$design$size, collapse = " ")), "\n\n", file = params_file, append = TRUE)
  
  # min. max. DE
  cat("### min. max. DE", "\n", file = params_file, append = TRUE)
  cat(paste0("MIN_DE=", MIN_DE), "\n", file = params_file, append = TRUE)
  cat(paste0("MAX_DE= ", MAX_DE), "\n\n", file = params_file, append = TRUE)
  
  
  # Contrast matrix
  cat("### Contrast matrix for group comparison", "\n", file = params_file, append = TRUE)
  cat("CONTRAST= ", "\n", file = params_file, append = TRUE)
  for (i in 1:nrow(mapdia$contrasts)){
    cat(paste(mapdia$contrasts[i,], collapse = " "), "\n", file = params_file, append = TRUE)
  }
  cat("\n", file = params_file, append = TRUE)
  
  
  cat("### protein_level.txt", "\n", file = params_file, append = TRUE)
  cat(paste0("MAX_PEP_PER_PROT= ", MAX_PEP_PER_PROT), "\n\n", file = params_file, append = TRUE)
  
  message(paste0("MapDIA parameters written to '", params_file, "'"))
}




getmapdiaContrasts <- function(file){
  
  pf <- readLines(file)
  
  groups <- pf[grep("LABELS", pf)]
  groups <- strsplit(trimws(gsub("LABELS.", "", groups)), " ")[[1]]
  
  pf1 <- pf[grep("CONTRAST", pf)+1]
  ngroups <- nchar(gsub(" ", "", pf1))
  
  cm <- pf[grep("CONTRAST", pf)+1:ngroups]
  cmn <- t(sapply(strsplit(cm, " "), as.numeric))
  dimnames(cmn) <- list(groups, groups)
  
  return(cmn)
}










### Heatmaps ----


heatmap <- function(heatdata,
                    coltitle = NULL, rowtitle = NULL,
                    scale_rows = TRUE, scale_cols = FALSE, # Z-scaling
                    cluster_rows = FALSE, plot_row_clusters = FALSE, cluster_cols = TRUE, plot_col_clusters = TRUE, coldend_size = 0.1, rowdend_size = 0.1, # Clustering
                    top_df = NULL, top_df_colors = NULL, top_df_height = 1, plot_top_legend = FALSE, # Top annotation
                    left_df = NULL, left_df_colors = NULL, left_df_width = 1, plot_left_legend = FALSE, # Left annotation
                    bottom_df = NULL, bottom_df_colors = NULL, bottom_df_width = 1, plot_bottom_legend = FALSE, # Bottom annotation
                    right_df = NULL, right_df_colors = NULL, right_df_width = 1, plot_right_legend = FALSE, # Right annotation
                    matchnames = TRUE, 
                    field_border_color = NA, field_border_width = 0,
                    legend_separate_file = FALSE, anno_leg_size = 1, # Legend parameters
                    colorscale = NULL, colorscale_legend_title = "", col_leg_height = 1000, col_leg_width = 100, # Color scale (heatmap values)
                    plotfontsize = 10, rowlabsize = 10, collabsize = 10, annofontsize = 10, leg_fontsize = 10, use_bold = TRUE, # Font sizes
                    bg = NA, NAcolor = rgb(0.90, 0.90, 0.90), annogridcolor = rgb(0.25,0.25,0.25),
                    colnames_height = 0.2, rownames_width = 0.2,
                    scale_margins = c(1,1,1,1),
                    pval_mat = NULL, pval_vec = NULL, pval_tres = 0.05, res = 300,
                    width = 4000, height = 4000, pdffile = NULL, tifffile = NULL, pngfile = NULL, svgfile=NULL,
                    return_results = FALSE, ...){
  
  
  # Documentation: https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html
  require("circlize", quietly = TRUE)
  require("ComplexHeatmap", quietly = TRUE)
  require("dendsort", quietly = TRUE)
  
  ##### FUNCTIONS
  
  getAnnoOrder <- function(anno_df){
    # Workaround for:
    # To match colors, need to use characters
    # To order the levels in the legend, need to use factor levels (otherwise as returned by 'unique')
    
    anno_df <- droplevels(anno_df)
    anno_df_order <- sapply(anno_df, levels)
    tmp <- lapply(anno_df[,sapply(anno_df_order, is.null),drop=FALSE], unique)
    anno_df_order[names(tmp)] <- tmp
    
    return(anno_df_order)
  }
  
  ##### PARAMS
  
  # Use 4000 * 4000 as base size
  plotcex <- sqrt(width * height) / sqrt(4000 * 4000) # scalefactor
  wcex <- width / 4000
  hcex <- height / 4000
  
  leg_title_gap <- unit(0.5*hcex, "cm")
  field_borders <- gpar(col = field_border_color, lwd = field_border_width)
  margins = unit(c(0.5*hcex*scale_margins[1], 0.5*wcex*scale_margins[2], 0.5*hcex*scale_margins[3], 0.5*wcex*scale_margins[4]), "cm") #bottom, left, top, right
  
  
  
  ##### ANNOTATIONS
  
  top_annotation <- left_annotation <- bottom_annotation <- right_annotation <- NULL
  topLegend <- leftLegend <- bottomLegend <- rightLegend <- NULL
  plotlegend <- c(plot_top_legend, plot_left_legend, plot_bottom_legend, plot_right_legend)
  
  
  
  ### Top annotation ----
  if (!is.null(top_df)){
    
    top_df_order <- getAnnoOrder(top_df)
    top_df[] <- lapply(top_df, as.character)
    
    
    if(matchnames == TRUE){ top_df <- top_df[colnames(heatdata),] }
    
    if (is.null(top_df_colors)){
      
      top_df_colors <- lapply(col2list(top_df), function(tmp){
        cols <- rand_color(n = length(unique(tmp)), hue = "red")
        names(cols) <- unique(tmp)
        return(cols)
      })
    }
    
    top_annotation <- HeatmapAnnotation(df = top_df,
                                        annotation_height = unit(100*top_df_height*hcex*(72/300), "points"),
                                        gap = unit(0, "cm"),
                                        border = FALSE,
                                        gp = gpar(col = annogridcolor),
                                        annotation_name_gp = gpar(fontsize = annofontsize*plotcex, fontface = ifelse(use_bold, "bold", "plain")),
                                        simple_anno_size_adjust = TRUE,
                                        annotation_name_side = "left",
                                        show_legend = FALSE,
                                        col = top_df_colors)
    
    topLegend <- lapply(names(top_df_order), function(tmp){
      
      Legend(labels = top_df_order[[tmp]],
             grid_height = unit(0.5*anno_leg_size*plotcex, "cm"),
             grid_width = unit(0.5*anno_leg_size*plotcex, "cm"),
             title = tmp,
             gap = unit(1*plotcex, "cm"),
             title_gap = leg_title_gap,
             title_gp = gpar(fontsize = leg_fontsize*plotcex, fontface = ifelse(use_bold, "bold", "plain")),
             legend_gp = gpar(fill = top_df_colors[[tmp]][top_df_order[[tmp]]]),
             labels_gp = gpar(fontsize = leg_fontsize*plotcex*0.85))
    })
    
  }
  
  
  
  
  ### Left annotation ----
  if (!is.null(left_df)){
    
    left_df_order <- getAnnoOrder(left_df)
    left_df[] <- lapply(left_df, as.character)
    
    if(matchnames == TRUE){ left_df <- left_df[rownames(heatdata),] }
    
    if (is.null(left_df_colors)){
      left_df_colors <- lapply(col2list(left_df), function(tmp){
        cols <- rand_color(n = length(unique(tmp)), hue = "green")
        names(cols) <- unique(tmp)
        return(cols)
      })
    }
    
    left_annotation <- HeatmapAnnotation(df = left_df,
                                         annotation_width = left_df_width,
                                         which = "row",
                                         gap = unit(0, "cm"),
                                         border = FALSE,
                                         gp = gpar(col = annogridcolor),
                                         annotation_name_gp = gpar(fontsize = annofontsize*plotcex, fontface = ifelse(use_bold, "bold", "plain")),
                                         simple_anno_size_adjust = TRUE,
                                         annotation_name_side = "bottom",
                                         show_legend = FALSE,
                                         col = left_df_colors)
    
    leftLegend <- lapply(1:ncol(left_df), function(tmp){
      
      Legend(labels = left_df_order[[tmp]],
             grid_height = unit(0.5*anno_leg_size*plotcex, "cm"),
             grid_width = unit(0.5*anno_leg_size*plotcex, "cm"),
             title = tmp,
             gap = unit(1*plotcex, "cm"),
             title_gap = leg_title_gap,
             title_gp = gpar(fontsize = leg_fontsize*plotcex, fontface = ifelse(use_bold, "bold", "plain")),
             legend_gp = gpar(fill = left_df_colors[[tmp]][left_df_order[[tmp]]]),
             labels_gp = gpar(fontsize = leg_fontsize*plotcex*0.85))
    })
    
    
  }
  
  
  ### Bottom annotation ----
  if (!is.null(bottom_df)){
    
    bottom_df_order <- getAnnoOrder(bottom_df)
    bottom_df[] <- lapply(bottom_df, as.character)
    
    if(matchnames == TRUE){ bottom_df <- bottom_df[colnames(heatdata),] }
    
    if (is.null(bottom_df_colors)){
      bottom_df_colors <- lapply(col2list(bottom_df), function(tmp){
        cols <- rand_color(n = length(unique(tmp)), hue = "blue")
        names(cols) <- unique(tmp)
        return(cols)
      })
    }
    
    bottom_annotation <- HeatmapAnnotation(df = bottom_df,
                                           annotation_height = bottom_df_height,
                                           gap = unit(0, "cm"),
                                           border = FALSE,
                                           gp = gpar(col = annogridcolor),
                                           annotation_name_gp = gpar(fontsize = annofontsize*plotcex, fontface = ifelse(use_bold, "bold", "plain")),
                                           simple_anno_size_adjust = TRUE,
                                           annotation_name_side = "left",
                                           show_legend = FALSE,
                                           col = bottom_df_colors)
    
    bottomLegend <- lapply(1:ncol(bottom_df), function(tmp){
      
      Legend(labels = bottom_df_order[[tmp]],
             grid_height = unit(0.5*anno_leg_size*plotcex, "cm"),
             grid_width = unit(0.5*anno_leg_size*plotcex, "cm"),
             title = tmp,
             gap = unit(1*plotcex, "cm"),
             title_gap = leg_title_gap,
             title_gp = gpar(fontsize = leg_fontsize*plotcex, fontface = ifelse(use_bold, "bold", "plain")),
             legend_gp = gpar(fill = bottom_df_colors[[tmp]][bottom_df_order[[tmp]]]),
             labels_gp = gpar(fontsize = leg_fontsize*plotcex*0.85))
    })
    
  }
  
  
  ### Right annotation ----
  if (!is.null(right_df)){
    
    right_df_order <- getAnnoOrder(right_df)
    right_df[] <- lapply(right_df, as.character)
    
    if(matchnames == TRUE){ right_df <- right_df[rownames(heatdata),] }
    
    if (is.null(right_df_colors)){
      right_df_colors <- lapply(col2list(right_df), function(tmp){
        cols <- rand_color(n = length(unique(tmp)), hue = "purple")
        names(cols) <- unique(tmp)
        return(cols)
      })
    }
    
    right_annotation <- HeatmapAnnotation(df = right_df,
                                          annotation_width = right_df_width,
                                          which = "row",
                                          gap = unit(0, "cm"),
                                          border = FALSE,
                                          gp = gpar(col = annogridcolor),
                                          annotation_name_gp = gpar(fontsize = annofontsize*plotcex, fontface = ifelse(use_bold, "bold", "plain")),
                                          simple_anno_size_adjust = TRUE,
                                          annotation_name_side = "bottom",
                                          show_legend = FALSE,
                                          col = right_df_colors)
    
    rightLegend <- lapply(1:ncol(right_df), function(tmp){
      
      Legend(labels = right_df_order[[tmp]],
             grid_height = unit(0.5*anno_leg_size*plotcex, "cm"),
             grid_width = unit(0.5*anno_leg_size*plotcex, "cm"),
             title = tmp,
             gap = unit(1*plotcex, "cm"),
             title_gap = leg_title_gap,
             title_gp = gpar(fontsize = leg_fontsize*plotcex, fontface = ifelse(use_bold, "bold", "plain")),
             legend_gp = gpar(fill = right_df_colors[[tmp]][right_df_order[[tmp]]]),
             labels_gp = gpar(fontsize = leg_fontsize*plotcex*0.85))
      
    })
    
  }
  
  legendlist <- unlist(list(topLegend, leftLegend, bottomLegend, rightLegend)[plotlegend], recursive = FALSE)
  
  
  
  
  ### Scaling ----
  
  matScale <- function(data, scale_rows = scale_rows, scale_cols = scale_cols){
    
    names.org <- dimnames(data)
    
    if (scale_rows == TRUE){
      data <- t(apply(data, 1, function(tmp) as.numeric(scale(tmp))))
    }
    
    if (scale_cols == TRUE){
      data <- apply(data, 2, function(tmp) as.numeric(scale(tmp)))
    }
    
    dimnames(data) <- names.org
    return(data)
  }
  
  
  heatdata <- matScale(heatdata, scale_rows = scale_rows, scale_cols = scale_cols)
  
  ### Colors ----
  
  if (is.null(colorscale)){
    scale_range <- 10^round(-log10(diff(range(heatdata, na.rm = TRUE))/10))
    if (scale_range < 1){ scale_range <- 1 }
    
    if (any(heatdata < 0, na.rm = TRUE)){
      heatpal = colorRamp2( c(floor(min(heatdata, na.rm = TRUE)*scale_range)/scale_range,0,ceiling(max(heatdata, na.rm = TRUE)*scale_range)/scale_range), c("blue", "white", "red"))
    } else {
      heatpal <-  colorRamp2( c(0,ceiling(max(heatdata, na.rm = TRUE)*scale_range)/scale_range), c("white", "red"))
    }
    
  } else {
    heatpal <- colorscale
  }
  
  
  ### Legends ----
  legend.colorscale <- Legend(col_fun = heatpal,
                              direction = "vertical",
                              title_gap = leg_title_gap,
                              border = rgb(1,1,1),
                              at = attr(heatpal, "breaks"),
                              title = colorscale_legend_title,
                              labels_gp = gpar(fontsize = leg_fontsize*plotcex*0.85),
                              legend_height = unit(col_leg_height*1000*hcex*(72/300), "points"),
                              title_gp = gpar(fontsize = leg_fontsize*plotcex, fontface = "bold"),
                              grid_width = unit(col_leg_width*80*wcex*(72/300), "points"),
                              title_position = "topleft")
  
  
  if (!is.null(legendlist)){
    legends.packed <- packLegend(list = c(legend.colorscale, legendlist),
                                 row_gap = unit(10*hcex, "mm"),
                                 column_gap = unit(10*wcex, "mm"))
  } else {
    legends.packed <- packLegend(legend.colorscale, direction = "vertical")
  }
  
  
  
  ### Clustering ----
  row.clusters <- FALSE
  col.clusters <- FALSE
  
  # rows
  if (is.logical(cluster_rows)){
    if (cluster_rows == TRUE){
      message("Cluster rows...")
      row.clusters <- dendsort(hclust(dist(heatdata)))
    }
  } else if (class(cluster_rows) == "hclust"){
    row.clusters <- cluster_rows
  }
  
  # columns
  if (is.logical(cluster_cols)){
    if (cluster_cols == TRUE){
      message("Cluster columns...")
      col.clusters <- dendsort(hclust(dist(t(heatdata))))
    }
  } else if (class(cluster_cols) == "hclust"){
    col.clusters <- cluster_cols
  }
  
  
  
  
  ### P-values ----
  
  if (!is.null(pval_mat)){
    # Caution: This function is called later (when using 'draw', thus pval_mat can only be defined once in a heatmap list)
    signiFUN <- function(j, i, x, y, width, height, fill) {if (pval_mat[i,j] <= pval_tres){grid.points(x, y, pch = 8, size = unit(4*plotcex, "pt"))}}
  } else {
    signiFUN <- NULL
  }
  
  if (!is.null(pval_vec)){
    rownames(heatdata) <- paste0(rownames(heatdata), " (p=", as.character(round(pval_vec, 2)), ")")
    rownames(heatdata)[pval_vec <= pval_tres] <- paste0(rownames(heatdata)[pval_vec <= pval_tres], "*")
  }
  
  
  ### Heatmap ----
  
  hm <- Heatmap(matrix = heatdata,
                row_title = rowtitle,
                column_title = coltitle,
                rect_gp = field_borders,
                column_dend_height = unit(coldend_size, "npc"),
                row_dend_width = unit(rowdend_size, "npc"),
                show_column_dend = plot_col_clusters,
                show_row_dend = plot_row_clusters,
                column_names_max_height = unit(colnames_height, "npc"),
                row_names_max_width = unit(rownames_width, "npc"),
                na_col = NAcolor,
                show_heatmap_legend = FALSE,
                cell_fun = signiFUN,
                row_names_gp = gpar(fontsize = rowlabsize*plotcex),
                column_names_gp = gpar(fontsize = collabsize*plotcex),
                col = heatpal,
                cluster_rows = row.clusters,
                cluster_columns = col.clusters,
                top_annotation = top_annotation,
                left_annotation = left_annotation,
                bottom_annotation = bottom_annotation,
                right_annotation = right_annotation, ...)
  
  
  
  ### Print to file(s) ----
  
  if (legend_separate_file != TRUE){
    legends_to_draw <- legends.packed
  } else {
    legends_to_draw <- NULL
  }
  
  message("Drawing heatmap...")
  
  if (!is.null(pdffile)){
    pdf(file = paste0(pdffile, ".pdf"), width = width/res, height = height/res, pointsize = plotfontsize)
    ComplexHeatmap::draw(hm, annotation_legend_list = legends_to_draw, padding = margins)
    dev.off()
  }
  
  if (!is.null(tifffile)){
    tiff(filename = paste0(tifffile, ".tiff"), width = width, height = height, pointsize = plotfontsize, res = 300)
    ComplexHeatmap::draw(hm, annotation_legend_list = legends_to_draw, padding = margins)
    dev.off()
  }
  
  if (!is.null(pngfile)){
    png(filename = paste0(pngfile, ".png"), width = width, height = height, units = "px", pointsize = plotfontsize, res = 300, bg=bg, type = "cairo")
    ComplexHeatmap::draw(hm, annotation_legend_list = legends_to_draw, padding = margins)
    dev.off()
  }
  
  if (!is.null(svgfile)){
    svg(filename = paste0(pngfile, ".svg"), width = width*0.6, height = height, pointsize = plotfontsize, bg = bg)
    ComplexHeatmap::draw(hm, annotation_legend_list = legends_to_draw, padding = margins)
    dev.off()
  }
  
  
  if (legend_separate_file == TRUE){
    
    message("Drawing legend...")
    
    if (!is.null(pdffile)){
      pdf(file = paste0(pdffile, "_legend.pdf"), width = width*0.6/res, height = height/res, pointsize = plotfontsize)
      draw(legends.packed)
      dev.off()
    }
    
    if (!is.null(tifffile)){
      tiff(filename = paste0(tifffile, "_legend.tiff"), width = width*0.6, height = height, pointsize = plotfontsize, res = 300)
      draw(legends.packed)
      dev.off()
    }
    
    if (!is.null(pngfile)){
      png(filename = paste0(pngfile, "_legend.png"), width = width*0.6, height = height, units = "px", pointsize = plotfontsize, res = 300, bg = bg, type = "cairo")
      draw(legends.packed)
      dev.off()
    }
    
    if (!is.null(svgfile)){
      svg(filename = paste0(pngfile, "_legend.svg"), width = width*0.6/res, height = height/res, pointsize = plotfontsize, bg = bg)
      draw(legends.packed)
      dev.off()
    }
    
  }
  
  message("Finished.")
  
  if ( (is.null(pdffile) & is.null(tifffile) & is.null(pngfile)) | return_results == TRUE ){
    return(list(hm = hm, ha = legends.packed, clusters = list(row = row.clusters, col = col.clusters)))
  }
  
  
}







PTMSEAheatmap <- function(heatdata, heatdesign, file, heatpal, trts, heatpval, collabs = TRUE, top_df_colors, htgap = 0.8, legcol = 1, legdir = "vertical",  legend_separate_file = FALSE, lh = 4, lsq = 1, width = 4000, height = 4000, fontsize.small = 14, fontsize.large = 18, fontsize.cols = 9, fontsize.rows = 9, pval_tres = 0.05, sigsymbol = 19, sigsize = 2.7, sigcol = rgb(0.45, 0.45, 0.45)){
  
  ### 1. Compute row clusters for the overall data
  
  clustdata <- heatdata
  clustdata[is.na(clustdata)] <- 0
  row.clusters <- dendsort(hclust(dist(clustdata)))
  
  
  ### 2. Generate a single legend (coloramp) and top annotation
  
  legend.colorscale <- Legend(col_fun = heatpal,
                              title_gap = unit(0.5, "cm"),
                              direction = "vertical",
                              at = attr(heatpal, "breaks"),
                              title = "PTM-SEA NES",
                              title_gp = gpar(fontsize = fontsize.small, fontface = "bold"),
                              labels_gp = gpar(fontsize = fontsize.small),
                              legend_height = unit(lh, "cm"),
                              grid_width = unit(1, "cm"),
                              title_position = "topleft") #"leftcenter-rot"
  
  topLegend <- Legend(labels = unique(heatdesign$organoid),
                      ncol = legcol, 
                      grid_height = unit(lsq, "cm"),
                      grid_width = unit(lsq, "cm"),
                      title = "organoid",
                      gap = unit(1, "cm"),
                      title_gap = unit(0.5, "cm"),
                      title_gp = gpar(fontsize = fontsize.small, fontface = "bold"),
                      legend_gp = gpar(fill = top_df_colors$organoid[unique(heatdesign$organoid)]),
                      labels_gp = gpar(fontsize = fontsize.small))
  
  
  
  ### 3. Make the list of heatmaps
  
  hml <- lapply(trts, function(trt){
    
    tmpdata <- heatdata[ ,heatdesign$treatment == trt ]
    tmpdesign <- heatdesign[ colnames(tmpdata), ]
    pval_mat <- heatpval[rownames(heatdata) ,colnames(tmpdata) ]
    
    pval_mat[is.na(pval_mat)] <- 1
    assign(paste0("pval_mat_", trt), value = pval_mat)
    
    # this functions is called in the end (when using 'draw') and then uses the last saved value of pval_mat, so it has to be neamed dynamically
    signiFUN <- function(j, i, x, y, width, height, fill) {if ( eval(parse(text = paste0("pval_mat_", trt)))[i,j] <= pval_tres){grid.points(x, y, gp = gpar(col = sigcol), pch = sigsymbol, size = unit(sigsize, "pt"))}}
    
    top_df <- heatdesign[colnames(tmpdata),"organoid",drop=FALSE]
    
    top_HA <- HeatmapAnnotation(df = top_df,
                                show_annotation_name = FALSE,
                                gap = unit(1.5, "cm"),
                                border = FALSE,
                                gp = gpar(col = rgb(0,0,0)),
                                annotation_name_gp = gpar(fontsize = fontsize.small, fontface = "bold"),
                                simple_anno_size_adjust = TRUE,
                                annotation_name_side = "left",
                                show_legend = FALSE,
                                col = top_df_colors)
    
    if (collabs == TRUE){
      collabs <- colnames(tmpdata)
    } else {
      collabs <- rep("", ncol(tmpdata))
    }
    
    hm <- Heatmap(matrix = tmpdata,
                  column_labels = collabs, 
                  row_names_max_width = unit(11, "cm"),
                  row_dend_width = unit(2.5, "cm"),
                  column_title = trt,
                  column_title_gp = gpar(fontsize = fontsize.large, fontface = "bold"),
                  show_column_names = TRUE,
                  rect_gp = gpar(col = rgb(1,1,1), lwd = unit(2, "pt")),
                  na_col = rgb(1, 1, 1),
                  show_heatmap_legend = FALSE,
                  row_names_gp = gpar(fontsize = fontsize.rows),
                  column_names_gp = gpar(fontsize = fontsize.cols),
                  col = heatpal,
                  cluster_rows = row.clusters,
                  cluster_columns = FALSE,
                  top_annotation = top_HA,
                  cell_fun = signiFUN)
    
    return(hm)
  })
  
  
  ### 4. Combine heatmaps into heatmap list and add the legend tolegendlist
  # Properties like row order are taken from the first (main) heatmap
  
  hm_combined <- Reduce("+", hml)
  legendlist <- packLegend(list = c(legend.colorscale, topLegend), direction = legdir, gap = unit(2, "cm"))
  
  
  ### 5. Draw heatmap and legend list
  
  if (legend_separate_file == TRUE){
    
    png(filename = paste0(gsub(pattern = ".png", "", file), "__legend", ".png"), type = "cairo", units = "px", width = height, height = height, res = 300)
    ComplexHeatmap::draw(legendlist)
    dev.off()
    
    png(filename = file, type = "cairo", units = "px", width = width, height = height, res = 300)
    ComplexHeatmap::draw(hm_combined,
                         ht_gap = unit(0.8, "cm"),
                         row_dend_side = "left",
                         heatmap_legend_side = "left",
                         annotation_legend_side = "left")
    dev.off()
    
  } else {
    
    png(filename = file, type = "cairo", units = "px", width = width, height = height, res = 300)
    ComplexHeatmap::draw(hm_combined,
                         annotation_legend_list = legendlist,
                         ht_gap = unit(htgap, "cm"),
                         row_dend_side = "left",
                         heatmap_legend_side = "left",
                         annotation_legend_side = "left")
    dev.off()
    
  }
  
}






`%row<%` <- function(mat,tres){
  res <- rowSums(naf(mat <= tres)) > 0
  return(res)
}


matchNodes <- function(IG, phosphodata){
  
  ignames <- V(IG)$name
  
  gix <- which(colnames(phosphodata$annotation) == "gene")
  syix <- which(colnames(phosphodata$annotation) == "synonyms")
  tmpann <- phosphodata$annotation[!duplicated(phosphodata$annotation$gene),]
  
  signorname <- sapply(as.list(data.frame(t(tmpann))), function(tmp){
    if (!tmp[gix] %in% ignames){
      tmpsyn <- strsplit(tmp[syix], ",")[[1]]
      return(paste(tmpsyn[tmpsyn %in% ignames], collapse = ","))
    } else{
      return(tmp[gix])
    }
  })
  names(signorname) <- tmpann$gene
  #signorname[nchar(signorname)==0] <- tmpann$gene[nchar(signorname)==0]
  
  signorname <- signorname[!is.na(signorname)]
  # data.frame( SIGNOR = signorname[which(signorname != phosphodata$annotation$gene)],  ANNO = phosphodata$annotation$gene[which(signorname != phosphodata$annotation$gene)])
  
  return(signorname)
}


impute_groups <- function(data, design, group = "CRC", subgroup = "trt", subgroup_ctrl = "control", min_value = NULL, max_na_ctrl = 0, max_na_trt = 0, min_trts_for_ctrl_imp = 1, max_na_keep = 0.5){
  
  ### Function to perform vectorized replacement of missing values and removal of non-replicated values
  
  # For each CRC in design:
  # If both replicates of control samples are measured, impute missing values in treatment groups
  # If both replicates in one treatment group are measured, impute missing values in the control group
  # If there are too many missing values in any group, set all values to NA
  
  # max_na_ctrl: max. allowed fraction of missing values to consider a ctrl group as a "true" signal
  # max_na_trt: max. allowed fraction of missing values to consider a trt group as a "true" signal
  # min_trts_for_ctrl_imp: how many trts groups must be "true" to do imputation of the ctrl group?
  # max_na_keep: Maximum allowed fraction of missing values in a group
  
  rn <- rownames(data)
  data <- data.matrix(data)
  
  # Imputation function
  impFUN <- function(data, scale = 0.2){
    imp_vals <- matrixStats::rowMins(data, na.rm = TRUE) * scale
    imp_vals[!is.finite(imp_vals)] <- NA
    imp <- matrix(rep(imp_vals, ncol(data)), ncol = ncol(data))
    return(imp)
  }
  
  ### 1. Remove very low values
  if (!is.null(min_value)){
    message("Setting ", signif(sum(data < min_value, na.rm = TRUE)/prod(dim(data)) * 100, 3), "% of values to NA because they are below ", min_value)
    message("These values may be imputed later...")
    data[data < min_value] <- NA
  }
  
  ### 2. Pre-calculate substitute values using impFUN
  imp <- impFUN(data)
  
  ### 3. Determine which values to replace with imp
  do_imp <- set_na <- matrix(data = FALSE, nrow = nrow(data), ncol = ncol(data), dimnames = dimnames(data))
  
  # For each organoid:
  for (tmpgroup in unique(design[,group])){
    
    message("Processing ", tmpgroup)
    ix <- design[,group] == tmpgroup
    subgroups <- design[ix,subgroup]
    
    # Calculate fraction of missing values per group
    na_fraction <- sapply(unique(subgroups), function(tmptrt) rowMeans(is.na(data[,ix][,subgroups %in% tmptrt,drop=FALSE])) )
    ctrl_measured <- na_fraction[,subgroup_ctrl] == max_na_ctrl # determine if enough control samples are measured
    trt_measured <- subset(na_fraction, select = -control) == max_na_trt # determine if enough treatment samples are measured
    
    # If enough treatment replicates are measured, impute the controls
    impute_ctrl <- !ctrl_measured & rowSums(trt_measured) >= min_trts_for_ctrl_imp # row-index of controls (impute/not impute)
    do_imp[,ix][,subgroups == subgroup_ctrl] <- impute_ctrl # replace in all columns
    
    # If all control replicates are measured, impute treatments
    impute_trts <- ctrl_measured & !trt_measured # row-indexes of all trts (impute/not impute)
    for (trt in subgroups[subgroups != subgroup_ctrl]){
      do_imp[,ix][,trt == subgroups] <- impute_trts[,trt] # replace in all columns
    }
    
    # If NA fraction per group is too high, set all values of that group to NA - after imputation
    na_fraction_imp <- sapply(unique(subgroups), function(tmp) rowMeans(!do_imp[,ix][,subgroups == tmp,drop=FALSE] & is.na(data[,ix][,subgroups == tmp,drop=FALSE] )) )
    for (trt in unique(subgroups)){
      set_na[,ix][,trt == subgroups] <- na_fraction_imp[,trt] > max_na_keep
    }
    
  }
  
  ### 4. Replace the missing values with imputed values and remove filtered values
  do_imp <- do_imp & is.na(data) # only replace values that are NA
  data[do_imp] <- imp[do_imp]
  data[set_na] <- NA
  rownames(data) <- rn
  
  message(signif(sum(do_imp)/prod(dim(data))*100, 2), "% of values were imputed")
  message(signif(sum(set_na)/prod(dim(data))*100, 2), "% of values were set to NA")
  
  return(list("data" = data, "impute" = do_imp, "remove" = set_na))
}





uniprot2gene <- function(ids, sep = ",", targetid = "GENES"){
  
  # Function to convert uniprot ids (possibly separated by ";") to gene symbols
  
  require(UniProt.ws, quietly = TRUE)
  
  if (!exists("uniprot.hs")){
    uniprot.hs <- UniProt.ws(taxId=9606)
  }
  
  ids.uni <- unique(ids)
  
  genelib <- UniProt.ws::select(uniprot.hs, ids.uni, targetid)
  colnames(genelib) <- c("UNIPROTKB", "GENES")
  
  genelib$gene <- sapply(strsplit(split = " ", genelib$GENES), function(tmp) tmp[1])
  genelib$synonyms <- sapply(strsplit(split = " ", genelib$GENES), function(tmp) paste(tmp, collapse = sep))
  
  multiprot <- strsplit(split = ";", ids.uni[ grep(pattern = ";", ids.uni) ])
  names(multiprot) <- ids.uni[ grep(pattern = ";", ids.uni) ]
  
  multigenes <- lapply(multiprot, function(tmpids){
    tryCatch(
      expr = {
        return(select(uniprot.hs, tmpids, targetid))
      },
      error = function(e){
        return(NA)
      }
    )
  })
  
  multigenes.gene <- sapply(multigenes, function(tmpgenes) {
    if(!all(is.na(tmpgenes))){
      return(paste(sapply(strsplit(split = " ", tmpgenes$GENES[!is.na(tmpgenes$GENES)]), function(tmp) tmp[1] ), collapse = ";") )
    } else {
      return(NA)
    }
  })
  
  multigenes.synonyms <- sapply(multigenes, function(tmpgenes) {
    if(!all(is.na(tmpgenes))){
      return(paste(gsub(pattern = " ", replacement = sep, tmpgenes$GENES[!is.na(tmpgenes$GENES)]), collapse = ";"))
    } else {
      return(NA)
    }})
  
  genelib[ match(names(multigenes.gene), genelib$UNIPROTKB), ]$gene <- multigenes.gene
  genelib[ match(names(multigenes.gene), genelib$UNIPROTKB), ]$synonyms <- multigenes.synonyms
  
  rownames(genelib) <- genelib$UNIPROTKB
  res <- genelib[,c("UNIPROTKB", "gene","synonyms")]
  colnames(res) <- c("uniprot", "gene","synonyms")
  
  return(res)
}





identifyPsites <- function(phosphoproteins, phosphopeptides, filter = "phospho", modpattern = "\\[.*?\\]", flankup = 7, flankdown = 7, sep = "|", uniprot.obj = NULL){
  
  if (!"UniProt.ws" %in% installed.packages()){
    BiocManager::install("UniProt.ws")
  }
  
  
  # ---- Functions ----
  
  getUniProt <- function(ids, uniprot.obj){
    # Download protein sequence from UniProt
    
    if (!"UniProt.ws" %in% installed.packages()){
      BiocManager::install("UniProt.ws") # in case it needs to be installed on each node
    }
    require(UniProt.ws, quietly = TRUE)
    
    # if (!exists("uniprot.obj")){
    #   uniprot.obj <- UniProt.ws(taxId = taxname2taxid("HUMAN"))
    #   message("'uniprot.obj' added to environment")
    # }
    
    k <- 0
    repeat({
      Sys.sleep(0.5)
      k <- k+1
      print(k)
      tryCatch({
        res <- select(uniprot.obj, keys = ids, columns = "SEQUENCE", keytype = "UNIPROTKB")
        break()
      }, error = function(err){
      })
      
      if (k > 10){
        res <- NULL
        stop("Error in select(uniprot.obj, ...)")
        break()
      }
    })
    
    res <- res[!is.na(res$UNIPROTKB),]
    
    tryCatch({
      rownames(res) <- res$UNIPROTKB
    }, error=function(err){
      message("Could not set rownames.")
    })
    
    return(res)
  }
  
  
  getPsites <- function(peptide_seq, protein_seq){
    
    ### Function that returns phosphosite positions
    
    if (is.na(peptide_seq) | is.na(protein_seq)){
      return(NA)
    }
    
    tmp.sequence <- gsub(modpattern, "", peptide_seq)
    
    # get modification names
    tmp.mod <- regmatches( peptide_seq, gregexpr(modpattern, peptide_seq))[[1]]
    
    # get modified amino acids
    tmp.split <- unlist(strsplit(peptide_seq, split = modpattern))
    if (length(tmp.mod) < length(tmp.split)){
      tmp.res <- tmp.split[-length(tmp.split)] # if mod is not at the last position
    } else {
      tmp.res <- tmp.split # if mod is at the last position
    }
    tmp.aa <- substr(tmp.res, nchar(tmp.res), nchar(tmp.res))
    
    # get positions
    tmp.pos <- cumsum(nchar(tmp.res))
    peptide_PTM_pos <- tmp.pos-1
    startpos_peptide <- unlist(gregexpr(pattern = tmp.sequence, protein_seq )) # get PTM position from string comparison between protein and peptide
    protein_PTM_pos <- sapply(startpos_peptide, function(tmp) tmp + peptide_PTM_pos)
    names(protein_PTM_pos) <- NULL
    
    aa <- sapply(tmp.pos, function(pos){ substr(tmp.sequence, start = pos, stop = pos) } ) # get PTM aa
    protein_PTM_sites <- paste(aa, protein_PTM_pos, sep = "")
    
    # filter for phosphosites
    if (!is.null(filter)){
      protein_PTM_sites <- protein_PTM_sites[grep(filter, tmp.mod, ignore.case = TRUE)]
    }
    
    return(paste(unique(protein_PTM_sites), collapse = sep))
  }
  
  
  
  getFlankingSeqs <- function(peptide_seq, protein_seq, flankup = 7, flankdown = 7){
    
    ### Function that returns flanking sequences
    
    if (is.na(peptide_seq) | is.na(protein_seq)){
      return(NA)
    }
    
    tmp.sequence <- gsub(modpattern, "", peptide_seq)
    
    # get modification names
    tmp.mod <- regmatches( peptide_seq, gregexpr(modpattern, peptide_seq))[[1]]
    
    # get modified amino acids
    tmp.split <- unlist(strsplit(peptide_seq, split = modpattern))
    if (length(tmp.mod) < length(tmp.split)){
      tmp.res <- tmp.split[-length(tmp.split)] # if mod is not at the last position
    } else {
      tmp.res <- tmp.split # if mod is at the last position
    }
    tmp.aa <- substr(tmp.res, nchar(tmp.res), nchar(tmp.res))
    
    # get positions
    tmp.pos <- cumsum(nchar(tmp.res))
    peptide_PTM_pos <- tmp.pos-1
    startpos_peptide <- unlist(gregexpr(pattern = tmp.sequence, protein_seq )) # get PTM position from string comparison between protein and peptide
    protein_PTM_pos <- sapply(startpos_peptide, function(tmp) tmp + peptide_PTM_pos)
    names(protein_PTM_pos) <- NULL
    
    # get flanking sequence(s)
    ext_right <- nchar(protein_seq) - (protein_PTM_pos + flankup)
    ext_right[ext_right>0] <- 0
    ext_left <- protein_PTM_pos - flankdown + 1
    ext_left[ext_left>0] <- 0
    flankseqs <- sapply(protein_PTM_pos, function(ptmpos) substr(protein_seq, start = ptmpos - flankdown, stop = ptmpos + flankup))
    flankseqs <- sapply(seq_along(flankseqs), function(i) paste0(paste(rep("_", -ext_left[i]), collapse = ""), flankseqs[i], paste(rep("_", -ext_right[i]), collapse = "")) )
    
    # filter for phosphosites
    if (!is.null(filter)){
      flankseqs <- flankseqs[grep(filter, tmp.mod, ignore.case = TRUE)]
    }
    
    return(paste(unique(flankseqs), collapse = sep))
    
  }
  
  
  
  
  # ---- Main ----
  
  if (is.null(uniprot.obj)){
    uniprot.obj <- UniProt.ws(taxId = taxname2taxid("HUMAN"))
  }
  
  ids.ppep <- paste(phosphoproteins, phosphopeptides, sep = "_")
  names(phosphoproteins) <- ids.ppep
  names(phosphopeptides) <- ids.ppep
  
  ids <- unique(phosphoproteins)
  
  message("Going parallel...")
  require("parallel")
  
  cmax <- detectCores()
  cmax <- cmax - ceiling(0.1*cmax)
  
  minchunk <- 100
  nclust <- min(cmax, ceiling(length(ids)/minchunk))
  ids.parts <- split(ids, ceiling(seq_along(ids)/ceiling(length(ids)/nclust)))
  
  
  # execute parallel
  cluster <- makeCluster(nclust)
  
  message("Downloading sequences...")
  seq.part <- parLapply(cluster, ids.parts, getUniProt, uniprot.obj = uniprot.obj)
  
  stopCluster(cluster)
  
  uniprotdata <- c()
  for (i in 1:length(seq.part)){
    uniprotdata <- rbind(uniprotdata, seq.part[[i]])
  }
  
  alldata <- uniprotdata[phosphoproteins,]
  protein_sequences <- alldata$SEQUENCE
  names(protein_sequences) <- ids.ppep
  
  # get PTM positions
  message("Computing PTM positions...")
  psites <- sapply(ids.ppep, function(id){ getPsites(peptide_seq = phosphopeptides[id], protein_seq = protein_sequences[id])} )
  psites <- psites[ids.ppep]
  
  # get flanking sequences
  message("Computing flanking sequences...")
  flankseqs <- sapply(ids.ppep, function(id){ getFlankingSeqs(peptide_seq = phosphopeptides[id], protein_seq = protein_sequences[id], flankup = flankup, flankdown = flankdown)} )
  flankseqs <- flankseqs[ids.ppep]
  
  psites[ psites == "numeric(0)" ] <- NA
  psites[ grep("NA", psites) ] <- NA
  psites[nchar(psites) == 0] <- NA
  
  message("Finished.")
  return(list("psites" = psites, "flankseqs" = flankseqs))
  
}



expandPsites <- function(sitenames){
  
  splitprot <- strsplit(sitenames, split = ";", fixed = TRUE)
  baseids <- lapply(splitprot, function(tmp) strsplit(tmp, split = "_", fixed = TRUE) )
  newids <- lapply(baseids, function(splitids){
    sapply(splitids, function(tmp){
      mod <- strsplit(tmp[2], split = "|", fixed = TRUE)[[1]]
      paste(paste(tmp[1], mod, sep = "_"), collapse = "|")
    })
  })
  res <- sapply(newids, function(tmp) paste(tmp, collapse = ";"))
  
  return(res)
}


selectMultiples <- function(v, vnames = NULL, FUN = max, ...){
  
  # v: vector of values based on which multiples are selected
  # FUN: function that selects which value to use
  # sel: logical index of whether a row is selected
  
  if (!is.null(vnames)){
    names(v) <- vnames
  }
  
  sel <- !logical(length(v))
  
  for (i in 1:length(v)){
    tmpdata <- v[names(v)[i] == names(v)] # get all duplicates/multiples
    tmpdata <- tmpdata[!is.na(tmpdata)]
    if (length(tmpdata) > 1){
      if (sum(tmpdata) == 0){ tmpdata <- Inf }
      sel[i] <- (v[i] == FUN(tmpdata, ...))[1] # select
    }
    
  }
  
  return(sel)
}


replaceInf <- function(data){
  data <- data.matrix(data)
  maxval <- ceiling(max(data[is.finite(data)], na.rm = TRUE))
  minval <- floor(min(data[is.finite(data)], na.rm = TRUE))
  data[naf(data > maxval)] <- maxval
  data[naf(data < minval)] <- minval
  return(data)
}


naf <- function(data){
  data[is.na(data)] <- FALSE
  return(data)
}


rowDemultiplex <- function(data, meandata = NULL, sep = "|"){
  
  # Split row identifiers and add as additional rows
  
  ids <- rownames(data)
  
  dmp <- unlist(sapply(seq_along(ids), function(tmp){
    idsep <- unique(strsplit(ids[tmp], split = sep, fixed = TRUE)[[1]])
    #idsep <- paste0(gsub("-p", "", idsep), "-p")
    if (length(idsep) > 1){
      idsep <- idsep[!idsep %in% ids] # remove all that are already on a unique peptide
    }
    res <- rep(tmp, length(idsep))
    names(res) <- idsep
    return(res)
  }))
  
  if (!is.null(meandata)){
    sel <- selectMultiples(v = meandata[dmp], vnames = names(dmp), FUN = max)
    dmp <- dmp[sel]
  }
  
  
  data.dmp <- data[dmp,]
  rownames(data.dmp) <- names(dmp)
  
  return(data.dmp)
}



fisherGS <- function(genes, universe, genesets, padjust = TRUE){
  
  # Test function
  getpval <- function(GSgenes, genes, universe){
    
    # Contingency table
    contable <- matrix(nrow=2, ncol=2)
    rownames(contable) <- c("genes", "universe")
    colnames(contable) <- c("GS", "notGS")
    contable["genes", "GS"] <- length(intersect(genes,GSgenes))
    contable["genes", "notGS"] <- length(setdiff(genes, GSgenes))
    contable["universe", "GS"] <- length(intersect(setdiff(universe,genes), GSgenes))
    contable["universe", "notGS"] <- length(setdiff(setdiff(universe,genes), GSgenes))
    
    # Fisher test
    pval <- fisher.test(contable, alternative="greater")$p.value
    return(pval)
  }
  
  pval <- sapply(genesets, getpval, genes = genes, universe = universe)
  
  if (padjust == TRUE){
    pval <- p.adjust(pval)
  }
  
  return(pval)
}



transformPhosphoGS <- function(GS, IG, activating_only = FALSE){
  
  phoshoGS <- lapply(GS, function(genes){
    
    psites <- unique(unlist(lapply( intersect(genes, V(IG)$name), function(tmpgene){
      res <- E(IG)$PSITE_NAME[ E(IG)$TARGET == tmpgene & E(IG)$EFFECT == "activation" ]
      res2 <- E(IG)$PSITE_NAME[ E(IG)$SOURCE == tmpgene & E(IG)$EFFECT == "activation" ]
      unique(c(res[!is.na(res)], res2[!is.na(res2)]))
    })))
    
    return(psites)
  })
  
  return(phoshoGS)
}


cutNames <- function(v, maxchar = 35){
  ix <- nchar(v) > maxchar
  v <- substr(v, 1, maxchar)
  v[ix] <- paste0(v[ix], "...")
  return(v)
}


matScale <- function(data, scale_rows = TRUE, scale_cols = FALSE, min = NULL, max = NULL){
  
  names.org <- dimnames(data)
  
  if (scale_rows == TRUE){
    data <- t(apply(data, 1, function(tmp) as.numeric(scale(tmp))))
  }
  
  if (scale_cols == TRUE){
    data <- apply(data, 2, function(tmp) as.numeric(scale(tmp)))
  }
  
  dimnames(data) <- names.org
  
  if (!is.null(max)) heatdata[heatdata > max] <- max
  if (!is.null(min)) heatdata[heatdata < min] <- min
  
  
  return(data)
}


ggvolcano <- function(data, x = NULL, y = NULL, color = NULL, label = NULL, shape = NULL,
                      nlabels = NULL, lab_size = 12, repel = 1.5, attract = NULL, box.padding = 0.5, max_overlaps = Inf, seed = 123,
                      ptres = 0.05, clip = FALSE, symlim = TRUE, expand = c(0,0), nbreaks_x = 7, nbreaks_y = 7,
                      color_up = "#eb9d0e", color_down = "#146bc7", color_nonsig = "#4d4d4d",
                      title = NULL, title_size = NULL, point_size = 2, scale_size = FALSE, axis_size = NULL, leg_size = NULL,
                      lwd = 0.8, at_zero = FALSE, ...){
  
  ### Function to plot volcano plots using ggplot.
  
  # SETUP
  require("ggplot2", quietly = TRUE)
  require("ggrepel", quietly = TRUE)
  require("dplyr", quietly = TRUE)
  require("rlang", quietly = TRUE)
  require("magrittr", quietly = TRUE)
  require("scales", quietly = TRUE)
  
  getLimits <- function(x, clip = TRUE, expand = xe){
    
    x <- x[!is.na(x)]
    x <- x + x*expand
    
    if (clip == TRUE){
      h <- hist(x, plot = FALSE, breaks = 30)
      xd <- h$counts > 3
      xmin <- h$breaks[which(xd)[1]]
      
      
      xmax <- rev(h$breaks)[which(rev(xd))[1]]
    } else {
      xmin <- NA
      xmax <- NA
    }
    
    if (is.na(xmin)){ xmin <- floor(min(x)*10^-floor(-log10(abs(min(x)))))/10^-floor(-log10(abs(min(x)))) }
    if (is.na(xmin)){xmin <- 0}
    if (xmin > 0){xmin <- 0}
    if (is.na(xmax)){ xmax <- ceiling(max(x)*10^-floor(-log10(abs(max(x)))))/10^-floor(-log10(abs(max(x)))) }
    
    return(c("min" = xmin, "max" = xmax))
  }
  
  
  
  # PARSE INPUTS
  
  data <- as.data.frame(data)
  
  x <- enquo(x)
  y <- enquo(y)
  
  if (quo_is_null(x)){ x <- sym(grep("lfc|log2FoldChange|logFC|log2FC|nes", names(data), value = TRUE, ignore.case = TRUE)[1]) }
  if (quo_is_null(y)){ y <- sym(grep("padj|fdr", names(data), value = TRUE, ignore.case = TRUE)[1]) }
  
  data$x <- data[[as_name(x)]]
  data$y <- -log10(data[[as_name(y)]])
  data <- data[!is.na(data$x) & !is.na(data$y),]
  data$score <- abs(as.numeric(scale(data$x, center = FALSE))) * abs(as.numeric(scale(data$y, center = FALSE)))
  
  data$class <- "not signif."
  data$class[data[[as_name(y)]] <= ptres & data$x > 0] <- "up"
  data$class[data[[as_name(y)]] <= ptres & data$x < 0] <- "down"
  
  
  if (is.null(title_size)) title_size <- lab_size
  if (is.null(axis_size)) axis_size <- lab_size
  if (is.null(leg_size)) leg_size <- lab_size
  
  
  shape <- enquo(shape)
  
  
  # LABELS
  
  label <- enquo(label)
  
  if (quo_is_null(label)){
    data[["label"]] <- rownames(data)
  } else {
    data[["label"]] <- data[[as_name(label)]]
  }
  
  data <- data[order(data$score, decreasing = TRUE),]
  
  if (is.null(nlabels)){ nlabels <- min(20, ceiling(nrow(data)/10)) }
  if (is.infinite(nlabels)){ nlabels <- nrow(data) }
  data$do_label <- FALSE
  nlabels_left <- ceiling(nlabels/2 * max(subset(data, x < 0)$score, na.rm = TRUE) / max(data$score, na.rm = TRUE))
  nlabels_right <- ceiling(nlabels/2 * max(subset(data, x > 0)$score, na.rm = TRUE) / max(data$score, na.rm = TRUE))
  data$do_label[data$x < 0][1:nlabels_left] <- TRUE
  data$do_label[data$x > 0][1:nlabels_right] <- TRUE
  
  
  if (sum(data$do_label) < nlabels){ data$do_label[!data$do_label][1:(nlabels-sum(data$do_label))] <- TRUE }
  data$label[!data$do_label] <- ""
  data$do_label[data$label == ""] <- FALSE
  
  # COLORS
  
  color <- enquo(color)
  
  if (quo_is_null(color)){
    color <- sym("class")
    colorvals <-  c("up" = color_up, "down" = color_down, "not signif." = color_nonsig)
  } else {
    col_levels <- unique(data[[as_name(color)]])
    colorvals <-  setNames(scales::muted(rainbow(length(col_levels))), col_levels)
  }
  
  
  # LIMITS
  
  xylimits <- list(xlim = getLimits(data$x, clip = clip, expand = expand[1]), ylim = getLimits(data$y, clip = clip, expand = expand[2]))
  if (symlim == TRUE){ xylimits$xlim <- c("min" = -max(abs(xylimits$xlim)), "max" = max(abs(xylimits$xlim))) }
  data$xorg <- data$x
  data$yorg <- data$y
  
  
  # set limits and axis tick labels
  # if clip==TRUE and any points are cut off, and recalculate limits based on breaks
  # if clip==TRUE and no points are cut off, use breaks as is
  # if clip==FALSE and no points are cut off, use breaks as is
  
  # X breaks
  xclip_min <- any(data$xorg < xylimits$xlim["min"])
  xclip_max <- any(data$xorg > xylimits$xlim["max"])
  
  xbreaks <- scales::pretty_breaks(n = nbreaks_x)(xylimits$xlim, n = nbreaks_x)
  if (clip & xclip_min){
    xylimits$xlim["min"] <- min(xbreaks)
    xclip_min <- any(data$xorg < xylimits$xlim["min"])
  }
  if (clip & xclip_max){
    xylimits$xlim["max"] <- max(xbreaks)
    xclip_max <- any(data$xorg > xylimits$xlim["max"])
  }
  
  
  xylimits$xlim <- xylimits$xlim + c(-diff(xylimits$xlim), diff(xylimits$xlim)) * c(!xclip_min, !xclip_max)*0.02
  
  names(xbreaks) <- as.character(xbreaks)
  if (xclip_min){ names(xbreaks)[1] <- paste0("<", xbreaks[1]) }
  if (xclip_max){ names(xbreaks)[length(xbreaks)] <- paste0(">", xbreaks[length(xbreaks)]) }
  
  
  
  # Y breaks
  yclip_min <- any(data$yorg < xylimits$ylim["min"])
  yclip_max <- any(data$yorg > xylimits$ylim["max"])
  
  ybreaks <- scales::pretty_breaks(n = nbreaks_y)(xylimits$ylim, n = nbreaks_y)
  
  if (clip & yclip_min){
    xylimits$ylim["min"] <- min(ybreaks)
    yclip_min <- any(data$xorg < xylimits$ylim["min"])
  }
  if (clip & yclip_max){
    xylimits$ylim["max"] <- max(ybreaks)
    yclip_max <- any(data$xorg > xylimits$ylim["max"])
  }
  
  
  xylimits$ylim <- xylimits$ylim + c(-diff(xylimits$ylim), diff(xylimits$ylim)) * c(!yclip_min & !at_zero, !yclip_max)*c(0.01, 0.05)
  
  names(ybreaks) <- as.character(ybreaks)
  if (yclip_min){ names(ybreaks)[1] <- paste0("<", ybreaks[1]) }
  if (yclip_max){ names(ybreaks)[length(ybreaks)] <- paste0(">", ybreaks[length(ybreaks)]) }
  
  
  data$x[ data$x < xylimits$xlim["min"] ] <- xylimits$xlim["min"]
  data$x[ data$x > xylimits$xlim["max"] ] <- xylimits$xlim["max"]
  data$y[ data$y < xylimits$ylim["min"] ] <- xylimits$ylim["min"]
  data$y[ data$y > xylimits$ylim["max"] ] <- xylimits$ylim["max"]
  
  
  
  ### GGPLOT ###
  
  data <- data[order(data$score, decreasing = FALSE),]
  
  gg <- data %>% ggplot(aes(x = x, y = y, label = label, color = !!color, shape = !!shape, ...))
  
  gg %<>% + theme_bw(base_size = 20)
  gg %<>% + theme(text = element_text(color = "black", size = lab_size),
                  rect = element_rect(color = "black", size = lwd),
                  line = element_line(size = lwd),
                  legend.text = element_text(color = "black", size = leg_size),
                  legend.title = element_text(color = "black", size = leg_size),
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_line(size = lwd, color = rgb(0.9,0.9,0.9)),
                  panel.border = element_rect(colour = "black", fill = NA, size = lwd),
                  strip.background = element_blank(),
                  strip.text = element_text(color = "black", size = title_size),
                  axis.ticks = element_line(color = "black", size = lwd),
                  axis.line = element_blank(),
                  plot.margin = unit(c(1,1,1,1), "cm"),
                  plot.title = element_text(size = title_size, hjust = 0.5, lineheight = 1.5),
                  axis.title = element_text(size = axis_size, face = "bold"),
                  axis.text = element_text(size = axis_size, color = "black"))
  
  if (!is.null(ptres)){ gg %<>% + geom_hline(yintercept = -log10(ptres), linetype = "dashed", color = rgb(0.3,0.3,0.3)) }
  
  # points
  if (scale_size == FALSE){
    gg %<>% + geom_point(size = point_size, alpha = 0.8)
  } else {
    gg %<>% + geom_point(aes(size = score), alpha = 0.8)
    gg %<>% + scale_size_continuous(range = c(point_size/5, point_size*2), guide = FALSE)
  }
  
  
  
  gg %<>% + scale_colour_manual(values = colorvals)
  gg %<>% + labs(title = title, y = paste0("-log10 ", as_name(y)), x = as_name(x), size = "none")
  
  
  gg %<>% + scale_x_continuous(expand = expansion(mult = c(0,0)),
                               limits = xylimits$xlim,
                               breaks = xbreaks,
                               labels = names(xbreaks))
  
  gg %<>% + scale_y_continuous(expand = expansion(mult = c(0,0)),
                               limits = xylimits$ylim,
                               breaks = ybreaks,
                               labels = names(ybreaks))
  
  # point labels
  if (is.null(attract)) attract <- sqrt(repel)
  gg %<>% + geom_text_repel(data = subset(data, do_label == TRUE),
                            size = lab_size/ggplot2:::.pt,
                            seed = seed,
                            xlim = xylimits$xlim - c(-diff(xylimits$xlim), diff(xylimits$xlim))*0.18,
                            ylim = xylimits$ylim - c(-diff(xylimits$ylim)*0.3, diff(xylimits$ylim)*0.02),
                            force = repel, force_pull = attract,  max.overlaps = max_overlaps,
                            point.padding = 0.35, box.padding = box.padding, max.time = 30, max.iter = 10^6, min.segment.length = 0, vjust = 0, color = rgb(0.0,0.0,0.0), segment.alpha = 0.6)
  
  gg %<>% + coord_cartesian(clip = "off")
  
  return(gg)
}




### NETWORKS -------------------------------------------------------------------



### DATA HANDLING AND PROCESSING ------------------------------------------------------

signor2net <- . %>%
  mutate(RESIDUE = sub("Ser", "S", RESIDUE)) %>% 
  mutate(RESIDUE = sub("Thr", "T", RESIDUE)) %>%
  mutate(RESIDUE = sub("Tyr", "Y", RESIDUE)) %>%
  mutate(from = ENTITYA, to = ENTITYB, PSITE = RESIDUE) %>%
  mutate(PSITE_NAME = paste0(ENTITYB, "_", PSITE)) %>%
  mutate(ID = paste(ENTITYA, ENTITYB, RESIDUE, sep = "_")) %>% mutate(ID = gsub('_$', '', ID)) %>%
  mutate(EFFECT = gsub("down.*", "inhibition", EFFECT)) %>%
  mutate(EFFECT = gsub("up.*", "activation", EFFECT)) %>%
  mutate(EFFECT = gsub("form complex", "binding", EFFECT)) %>%
  mutate(PSITE = "[<-"(PSITE, !MECHANISM %in% c("phosphorylation", "dephosphorylation"), "")) %>%
  mutate(MECHANISM = "[<-"(MECHANISM, !MECHANISM %in% names(ltys), "other")) %>%
  arrange(desc(SCORE)) %>%
  arrange(factor(EFFECT, ordered = TRUE, levels = c("activation", "inhibition", "binding", "unknown"))) %>%
  arrange(desc(nchar(PSITE))) %>%
  subset(!duplicated(ID)) %>%
  relocate(from, to)


importMutationdata <- function(mutfile, coding = TRUE){
  
  mutdata <- read.delim(mutfile)
  
  CRCs <- unique(mutdata$PatientID)
  names(CRCs) <- CRCs
  mutres <- lapply(CRCs, function(CRC){
    
    mutraw <- mutdata[ mutdata$PatientID == CRC, ]
    
    mdf <- data.frame(GENE = mutraw$Gene, TYPE = mutraw$variant_type, stringsAsFactors = FALSE)
    mgenes <- unique(mdf$GENE)
    
    mgenesnames <- sapply(mgenes, function(mgene) paste( unique(mdf[naf(mdf$GENE == mgene),]$TYPE), collapse = ","))
    names(mgenes) <- mgenesnames
    return(mgenes)
  })
  
  
  if (coding == TRUE){
    
    filter_vars <- c("missense_variant",
                     "frameshift_variant",
                     "splice_acceptor_variant",
                     "inframe_deletion",
                     "stop_gained",
                     "splice_donor_variant",
                     "start_lost",
                     "inframe_insertion",
                     "stop_lost",#
                     "coding_sequence_variant")
    
    mutres <- lapply(mutres, function(M){
      x <- strsplit(names(M), split = ",|&")
      y <- sapply(x, function(xx){ unique(xx[xx %in% filter_vars]) })
      z <- sapply(y, paste, collapse = ",")
      names(M) <- z
      M[nchar(z) > 0]
    })
  }
  
  
  return(mutres)
}



### NETWORK PROCESSING AND ANALYSIS ------------------------------------------------------


setGeneric(name = "select", def = function(.data, ...){ standardGeneric("select") } )

setMethod(f = "select", signature = "ANY", definition = function(.data, ...){
  dplyr::select(.data, ...)
})

setMethod(f = "select", signature = "tbl_graph", definition = function(.data, ..., env = parent.frame()){
  
  graph <- .data
  args <- enquos(...)
  
  args <- args[[length(args)]]
  
  ix <- rlang::eval_tidy(args, data = graph$design, env = env)
  if (!is.null(rownames(graph$design)) & !any(duplicated(rownames(graph$design)))) ix <- rownames(graph$design)[ix]
  
  if (all(ix == FALSE)) return(NULL)
  graph <- subset_data(graph, ix)
  if (!is.null(graph$design)) graph$design <- graph$design[ix,,drop = FALSE]
  return(graph)
  
})


setMethod(f = "select", signature = "igraph", definition = function(.data, ..., env = parent.frame()){
  graph <- .data
  args <- enquos(...)
  
  args <- args[[length(args)]]
  
  ix <- rlang::eval_tidy(args, data = graph$design, env = env)
  if (!is.null(rownames(graph$design)) & !any(duplicated(rownames(graph$design)))) ix <- rownames(graph$design)[ix]
  
  if (all(ix == FALSE)) return(NULL)
  graph <- subset_data(graph, ix)
  if (!is.null(graph$design)) graph$design <- graph$design[ix,,drop = FALSE]
  return(graph)
  
})


setMethod(f = "melt", signature = "igraph", definition = function(data, ...){
  
  # Melt edge data, return new graph with an edge for each column
  args <- list(...)
  graph <- data
  if (vcount(graph) <= 1 & ecount(graph) == 0) return(as_tbl_graph(graph))
  
  edges_long <- getdata(graph, type = "edges", args[[1]]) %>% melt(value.name = args[[1]])
  colnames(edges_long) <- c("ID", "sample", args[[1]])
  edges_long$ID <- as.character(edges_long$ID)
  edges_long[[args[[2]]]] <- getdata(graph, type = "edges", args[[2]]) %>% melt(value.name = args[[2]]) %>% subset(select = args[[2]]) %>% unlist()
  
  edf <- igraph::as_data_frame(graph)
  edges_long <- data.frame(edges_long, edf[match(edges_long$ID, edf$ID),])
  edges_long %<>% relocate(c(from, to))
  
  graph <- graph_from_data_frame(edges_long, vertices = igraph::as_data_frame(graph, what = "vertices")) %>% as_tbl_graph()
  return(graph)
})


subset_data <- function(graph, ix){
  
  # nodes
  n <- gsub("DATA__", "", grep("DATA__", vertex_attr_names(graph), value = TRUE))
  if (length(n) > 0){
    ndata <- lapply(setNames(n, n), getdata, graph = graph)
    ndata <- lapply(ndata, function(tmpdata) tmpdata[, ix, drop = FALSE] )
    for (i in 1:length(ndata)){
      
      graph <- do.call(setdata, c(list(graph = graph, fill = FALSE), ndata[i]))
      if (ncol(ndata[[i]]) == 1){
        vertex_attr(graph = graph, name = names(ndata)[i], index = V(graph)) <- setNames(as.vector(ndata[[i]]), V(graph)$name)
      }
    }
  }
  
  # edges
  e <- gsub("DATA__", "", grep("DATA__", edge_attr_names(graph), value = TRUE))
  if (length(e) > 0){
    edata <- lapply(setNames(e, e), getdata, graph = graph, type = "edges")
    edata <- lapply(edata, function(tmpdata) tmpdata[, ix, drop = FALSE] )
    for (i in 1:length(edata)){
      graph <- do.call(setdata, c(list(graph = graph, type = "edges", fill = FALSE), edata[i]))
      if (ncol(edata[[i]]) == 1){
        edge_attr(graph = graph, name = names(edata)[i], index = E(graph)) <- setNames(as.vector(edata[[i]]), E(graph)$ID)
      }
    }
  }
  
  return(graph)
}



tolist <- function(data){
  itemlist <- as.list(as.data.frame(t(data)))
  itemlist <- lapply(itemlist, setNames, colnames(data))
  return(itemlist)
}



fromlist <- function(list){
  cols <- unique(unlist(lapply(list, names)))
  cols <- cols[!is.na(cols)]
  list[sapply(list, is.null)] <- NA
  if (!is.null(cols)) list <- lapply(list, function(tmp) tmp[match(cols, names(tmp))] )
  data <- t(as.data.frame(list))
  rownames(data) <- NULL
  return(data)
}



setdata <- function(graph, type = "nodes", fill = TRUE, ...){
  
  # Set dataframe as node/edge data
  
  args <- list(...)
  stopifnot(length(args) == 1)
  names(args) <- paste0("DATA__", names(args))
  
  data <- args[[1]]
  if (!is.null(graph$design) & fill){
    s <- rownames(graph$design)
    nots <- s[!s %in% colnames(data)]
    nadf <- data.frame(matrix(NA, nrow = nrow(data), ncol = length(nots), dimnames = list(rownames(data), nots)))
    data <- cbind(data, nadf)[,rownames(graph$design)]
  }
  
  if (type == "nodes"){
    vertex_attr(graph = graph, name = names(args), index = V(graph)) <- tolist(data)[V(graph)$name]
  }
  if (type == "edges"){
    stopifnot(!is.null(E(graph)$ID))
    edge_attr(graph = graph, name = names(args), index = E(graph)) <- tolist(data)[E(graph)$ID]
  }
  
  return(graph)
  
}


getdata <- function(graph, type = "nodes", ...){
  
  # Get dataframe from node/edge data
  
  args <- list(...)
  if (length(args) == 0){ args[[1]] <- type; type <- "nodes" }
  stopifnot(length(args) == 1)
  args <- paste0("DATA__", args)
  
  if (type == "nodes"){
    data <- vertex_attr(graph = graph, name = args, index = V(graph))
    data <- fromlist(data)
    rownames(data) <- V(graph)$name
  }
  if (type == "edges"){
    stopifnot(!is.null(E(graph)$ID))
    data <- edge_attr(graph = graph, name = args, index = E(graph))
    data <- fromlist(data)
    rownames(data) <- E(graph)$ID
  }
  
  return(data)
  
}






lcc <- function(IG){
  
  # Function that returns the largest connected component of a network
  
  subgraphs <- decompose.graph(IG) # decompose graph into disconnected subgraphs/components
  subgraphs_nodes <- sapply(subgraphs, vcount) # get the number of nodes for each subgraph
  
  if (any(subgraphs_nodes)){
    
    ixlargest <- (1:length(subgraphs_nodes))[subgraphs_nodes == max(subgraphs_nodes)]
    subGlargest <- subgraphs[[ixlargest[1]]] # select largest component
    
    # messages
    Vkept <- length(as_ids(V(subGlargest)))
    Vorig <- length(as_ids(V(IG)))
    
    message(paste("Number of disconnected components is ", length(subgraphs), ".", sep = ""))
    message(paste("The largest component has ", Vkept, " nodes.", sep = ""))
    message(paste(Vorig - Vkept, " nodes were removed.", sep = ""))
    
    return(subGlargest)
    
  } else {
    print("No nodes found in graph.")
    return(graph.empty())
  }
}


rgraph <- function(n = 20){
  
  # Generate random graph with n nodes
  rnodes <- paste0(LETTERS, rep(1:3, each = length(LETTERS)))[1:n]
  graph <- play_erdos_renyi(n = n, p = 0.1, loops = TRUE)
  edf <- igraph::as_data_frame(graph)
  edf <- rbind(edf, edf[sample(1:nrow(edf), round(nrow(edf)*0.3)),])
  graph <- graph_from_data_frame(edf) %>% as_tbl_graph()
  graph %<>% activate(edges) %>% uniqueEID() %>% mutate(ID = paste0("edge_", ID))
  graph %<>% activate(nodes) %>% mutate(name = rnodes[1:vcount(graph)])
  
  return(graph)
}







subnet <- function(graph, nodes = NULL, edges = NULL, weights = NULL, targets = NULL, max_edges = NULL, direction = "all", edge_direction = "all", edge_order = 1){
  
  ## Function for subgraph extraction based on shortest paths between selected nodes/edges
  # nodes: node attribute (logical vector) for selection of nodes
  # edges: edge attribute (logical vector) for selection of edges
  # weights: edge attribute used as weight
  # targets: optional target nodes to/from which to calculate shortest paths
  # max_edges: max. search radius for finding shortest paths between nodes
  # direction: whether to use undirected ("all"), incoming ("in") or outgoing ("out") edges for shortest paths between nodes
  
  # edges must have 'ID' attribute
  
  # 1. Input
  nodes <- rlang::enquo(nodes)
  edges <- rlang::enquo(edges)
  weights <- rlang::enquo(weights)
  
  ndf <- igraph::as_data_frame(graph, what = "vertices")
  edf <- igraph::as_data_frame(graph, what = "edges")
  
  if (!is.null(targets)){
    targets <- targets[targets %in% V(graph)$name]
    not_in_nodes <- targets[!targets %in% V(graph)$name] # use only targets that are also in graph
    if (length(not_in_nodes) > 0) message(paste0("Warning: ", paste(not_in_nodes, collapse = ", "), " not found in graph"))
  }
  
  
  if (!rlang::quo_is_null(nodes)) node_select <- dplyr::pull(.data = ndf, !!nodes) else node_select <- rep(NA, vcount(graph)) # get selected nodes
  if (!rlang::quo_is_null(edges)) edge_select <- dplyr::pull(.data = edf, !!edges) else edge_select <- rep(NA, ecount(graph)) # get selected edges
  node_select[is.na(node_select)] <- FALSE
  edge_select[is.na(edge_select)] <- FALSE
  nodes <- as_ids(V(graph)[node_select])
  stopifnot(!is.null(E(graph)$ID))
  edges_ids <- E(graph)[edge_select]$ID
  
  weight_vec <- NULL
  if (!rlang::quo_is_null(weights)) weight_vec <- dplyr::pull(.data = edf, !!weights) # get weights
  if (!is.null(weight_vec)) names(weight_vec) <- E(graph)$ID
  
  stopifnot(!any(duplicated(nodes)))
  stopifnot(!any(duplicated(edges_ids)))
  
  message(paste0("Subgraph extraction using ", length(nodes), " nodes and ", length(edges_ids), " edges..."))
  
  # 2. Get edges that are on shortest paths between the input nodes
  if (length(nodes) > 1){
    # get shortest paths as sequences of nodes:
    vseq <- shortestPaths(graph = graph, from = nodes, to = targets, direction = direction, max_edges = max_edges, weights = weight_vec)
    # get all edge IDs for node the sequences:
    spedges <- vseqEdges(vseq = vseq, graph = graph, use_weights = TRUE)
  } else {
    spedges <- NULL
  }
  
  
  # 3. Get selected edges that are adjacent to the subgraph:
  seledges <- NULL
  if (!is.infinite(edge_order) & edge_order != 0){
    tmpgraph <- subgraph.edges(graph, eids = match(spedges, E(graph)$ID))
    tmpnodes <- unique(unlist(sapply(neighborhood(graph = graph, order = edge_order, mode = edge_direction, nodes = V(tmpgraph)$name), as_ids)))
    adj_graph <- induced_subgraph(graph = graph, v = tmpnodes) # might add many edges between the new nodes...
    adj_edf <- igraph::as_data_frame(adj_graph, what = "edges")
    if (!rlang::quo_is_null(edges)) edge_select <- dplyr::pull(.data = adj_edf, !!edges) else edge_select <- rep(FALSE, ecount(adj_graph))
    seledges <- E(adj_graph)$ID[edge_select]
  } else if (edge_order != 0){
    seledges <- E(graph)$ID[edge_select]
  }
  
  # 4. Make subgraph from all selected edges
  all_edges <- unique(c(spedges, seledges))
  subgraph <- subgraph.edges(graph, eids = match(all_edges, E(graph)$ID))
  message(paste0("Subgraph contains ", vcount(subgraph), " nodes and ", ecount(subgraph), " edges."))
  
  return(subgraph)
}


simplifySubnet <- function(graph, edge_select = "DE", node_select = "DE", rm_loops = TRUE, rm_leafs = TRUE, rm_multiples = TRUE){
  
  ### Function to simplify graphs (remove loops and multiple edges between the same pairs of nodes) based on predefined rules.
  
  if (ecount(graph) < 2) return(graph)
  
  edf <- igraph::as_data_frame(graph)
  stopifnot(!is.null(edf$ID))
  if (any(duplicated(edf$ID))) edf$ID2 <- E(uniqueEID(graph))$ID else edf$ID2 <- edf$ID
  
  edf$uniqueSELECT <- TRUE
  if (!is.null(edf$PSITE_NAME)) edf$uniqueSELECT <- !edf$PSITE_NAME %in% edf$PSITE_NAME[duplicated(edf$PSITE_NAME)] & edf[[edge_select]]
  edf$from_to <- with(edf, paste0(from, "_", to))
  edges <- setNames(rep(TRUE, nrow(edf)), edf$ID2)
  
  # remove loops
  if (rm_loops == TRUE) edges <- edges  &  !( edf$from == edf$to & !edf[[edge_select]] )
  
  
  # remove multiples
  if (rm_multiples == TRUE){
    for (i in 1:length(edges)){
      tmpdf <- edf[edf$from_to[i] == edf$from_to,,drop = FALSE]
      if (nrow(tmpdf) > 1 & any(tmpdf[[edge_select]], na.rm = TRUE)) tmpdf <- tmpdf[tmpdf[[edge_select]],,drop = FALSE]
      if (nrow(tmpdf) > 1 & any(nchar(tmpdf$PSITE) > 0, na.rm = TRUE)) tmpdf <- tmpdf[nchar(tmpdf$PSITE) > 0,,drop = FALSE]
      if (nrow(tmpdf) > 1 & !is.null(tmpdf$EFFECT)) tmpdf <- tmpdf[!tmpdf$EFFECT %in% c("unknown", "other"),,drop = FALSE]
      if (nrow(tmpdf) > 1 & length(grep("DATA_", colnames(edf))) > 1){
        measured <- rowSums(!is.na(edf[,grep("DATA_", colnames(edf))]))
        tmpdf <- tmpdf[measured == max(measured),,drop = FALSE]
      }
      edges[i] <- edges[i]  &  ( edf$ID2[i] %in% tmpdf$ID2 )
    }
  }
  
  sg <- subgraph.edges(graph = graph, eids = E(graph)[edges])
  
  if (rm_multiples == TRUE) sg <- de_edges_2(sg, edge_select = edge_select)
  
  
  # remove non-selected leafs with duplicate p-sites
  
  edf <- igraph::as_data_frame(sg)
  edf$uniqueSELECT <- TRUE
  if (!is.null(edf$PSITE_NAME)) edf$uniqueSELECT <- !edf$PSITE_NAME %in% edf$PSITE_NAME[duplicated(edf$PSITE_NAME)] & edf[[edge_select]]
  edf$from_to <- with(edf, paste0(from, "_", to))
  edges <- setNames(rep(TRUE, nrow(edf)), edf$ID2)
  
  if (rm_leafs == TRUE){
    nn <- sapply(setNames(neighborhood(sg, order = 1), V(sg)$name), length)-1
    leafs <- V(sg)$name[ V(sg)$name %in% names(nn)[nn == 1] & !vertex_attr(sg, node_select) ]
    edges <- edges  &  !( (edf$from %in% leafs | edf$to %in% leafs) & !edf$uniqueSELECT)
  }
  
  sg <- subgraph.edges(graph = sg, eids = E(sg)[edges])
  
  return(sg)
}


getCRCunion <- function(IG, subgraphs){
  
  ### Function to merge subgraphs by CRC
  
  # Group by CRC
  CRCs <- as.character(IG$design$CRC)
  IGLcrc <- lapply(unique(CRCs), function(CRCtmp) return((subgraphs[CRCs == CRCtmp]) ))
  names(IGLcrc) <- unique(CRCs)
  
  # Make graph for each treatment
  V(IG)$ID <-  V(IG)$name
  
  res <- lapply(IGLcrc, function(CRCgraphs){
    
    UGcrc <- lcc(do.call(igraph::union, args = c(CRCgraphs, byname = TRUE)))
    V(UGcrc)$ID <- V(UGcrc)$name
    
    UGcrc <-  transfer_attr(IGfrom = IG, IGto = UGcrc, delete = TRUE)
    
    CRCres <- lapply(CRCgraphs, function(TRTgraph){
      
      tmpIG <- UGcrc
      
      eattr <- edge_attr_names(TRTgraph)[ grep(paste(c("DE","PVAL", "PADJ", "LFC", "NES"), collapse = "|"), edge_attr_names(TRTgraph), ignore.case = TRUE) ]
      vattr <- vertex_attr_names(TRTgraph)[ grep(paste(c("DE","PVAL", "PADJ", "LFC", "MUT", "NES"), collapse = "|"), vertex_attr_names(TRTgraph), ignore.case = TRUE) ]
      
      tmpIG <- transfer_attr(IGfrom = TRTgraph, IGto = tmpIG, sel_eattr = eattr, sel_vattr = vattr)
      
      V(tmpIG)$orig <- V(tmpIG)$ID %in% V(TRTgraph)$ID
      E(tmpIG)$orig <- E(tmpIG)$ID %in% E(TRTgraph)$ID
      
      attr(tmpIG, "CRC") <- attr(TRTgraph, "CRC")
      attr(tmpIG, "TRT") <- attr(TRTgraph, "TRT")
      
      return(tmpIG)
    })
    
    return(list("crc" = UGcrc, "trts" = CRCres))
  })
  
  return(res)
  
}


transfer_attr <- function(IGfrom, IGto, sel_eattr = NULL, sel_vattr = NULL, delete = TRUE){
  
  ### Function to transfer attributes from one graph to another
  
  if (is.null(sel_eattr)){ sel_eattr <- edge_attr_names(IGfrom)[!edge_attr_names(IGfrom) %in% c("name", "ID")] }
  if (is.null(sel_vattr)){ sel_vattr <- vertex_attr_names(IGfrom)[!vertex_attr_names(IGfrom) %in% c("name", "ID")] }
  
  if (delete == TRUE){ # delete old attributes from IGto
    for (attr in edge_attr_names(IGto)[!edge_attr_names(IGto) %in% c("name", "ID")]) {IGto <- delete_edge_attr(IGto, attr)}
    for (attr in vertex_attr_names(IGto)[!vertex_attr_names(IGto) %in% c("name", "ID")]) {IGto <- delete_vertex_attr(IGto, attr)}
  } else {
    for (attr in intersect(sel_eattr, edge_attr_names(IGto))) {IGto <- delete_edge_attr(IGto, attr)}
    for (attr in intersect(sel_vattr, vertex_attr_names(IGto))) {IGto <- delete_vertex_attr(IGto, attr)}
  }
  
  if (is.null(V(IGfrom)$ID) | is.null(V(IGto)$ID)){
    message("Warning: Cannot use vertex ids.")
    Vids.from <- as_ids(V(IGfrom))
    Vids.to <- as_ids(V(IGto))
  } else {
    Vids.from <- V(IGfrom)$ID
    Vids.to <- V(IGto)$ID
  }
  
  if (is.null(E(IGfrom)$ID) | is.null(E(IGto)$ID)){
    message("Warning: Cannot use edge ids.")
    Eids.from <- as_ids(E(IGfrom))
    Eids.to <- as_ids(E(IGto))
  } else {
    Eids.from <- E(IGfrom)$ID
    Eids.to <- E(IGto)$ID
  }
  
  
  # edge attributes
  eattr.from <- edge_attr(IGfrom)[sel_eattr]
  eattr.to <- lapply(eattr.from, function(tmpattr){
    names(tmpattr) <- Eids.from
    tmpattr.to <- tmpattr[Eids.to]
    names(tmpattr.to) <- Eids.to
    return(tmpattr.to)
  })
  edge_attr(IGto) <- c(edge_attr(IGto), eattr.to)
  
  # vertex attributes
  vattr.from <- vertex_attr(IGfrom)[sel_vattr]
  vattr.to <- lapply(vattr.from, function(tmpattr){
    names(tmpattr) <- Vids.from
    tmpattr.to <- tmpattr[Vids.to]
    names(tmpattr.to) <- Vids.to
    return(tmpattr.to)
  })
  vertex_attr(IGto) <- c(vertex_attr(IGto), vattr.to)
  
  return(IGto)
}


shortestPaths <- function(graph, from, to = NULL, max_edges = NULL, direction = "all", weights = NULL){
  
  ### Function returning all shortest path from nodes 'from' to nodes 'to' using edge IDs (E(graph)$ID).
  
  stopifnot(!is.null(E(graph)$ID) | any(duplicated(E(graph)$ID)))
  gin <- graph
  if (is.null(to)) to <- from
  
  
  vseq <- lapply(setNames(from, from), function(tmpnode){
    if (!is.null(max_edges)) gin <- make_ego_graph(graph, order = max_edges, mode = direction, nodes = tmpnode)[[1]]
    if (is.null(weights)) w <- NA else w <- weights[E(gin)$ID]
    
    sp <- all_shortest_paths(gin, from = tmpnode, to = to[to != tmpnode & to %in% V(gin)$name],
                             mode = direction,
                             weights = w)
    
    res <- lapply(sp$res, as_ids)
    if (!is.null(res)) names(res) <- sapply(res, function(tmp) paste(tmp[1], tmp[length(tmp)], sep = "_") )
    return(res)
  })
  vseq <- unique(unlist(vseq, recursive = FALSE))
  names(vseq) <- sapply(vseq, function(tmp) paste(tmp[1], tmp[length(tmp)], sep = "_") )
  
  return(vseq)
}


uniqueEID <- function(graph){
  
  edf <- igraph::as_data_frame(graph)
  edf$dupIDs <- paste0(edf$from, "_", edf$to)
  
  edf$duplicated <- edf$dupIDs %in% unique(edf$dupIDs[duplicated(edf$dupIDs)])
  edf$running <- 0
  for (dupID in unique(edf$dupIDs)){
    edf$running[edf$dupIDs == dupID] <- cumsum(subset(edf, dupIDs == dupID)$duplicated)
  }
  
  edf$uniIDs <- paste0(edf$dupIDs, "_", edf$running)
  edf$uniIDs <- sub("_0", "", edf$uniIDs)
  
  stopifnot(!any(duplicated(edf$uniIDs)))
  
  E(graph)$ID <- edf$uniIDs
  return(graph)
}


nodeMatrix <- function(graphlist, attr = NULL){
  
  allnodes <- sort(unique(unlist(lapply(graphlist, function(tmp) V(tmp)$name))))
  res <- matrix(data = NA, nrow = length(allnodes), ncol = length(graphlist), dimnames = list(allnodes, names(graphlist)))
  
  if (!is.null(attr)){
    for (i in 1:length(graphlist)){
      tmp <- graphlist[[i]]
      res[,i] <- vertex_attr(tmp, attr)[match(allnodes, V(tmp)$name)]
    }
  } else {
    for (i in 1:length(graphlist)){
      res[,i] <- allnodes %in% V(graphlist[[i]])$name
    }
  }
  
  return(res)
}



de_edges_2 <- function(graph, edge_select = "DE"){
  
  edf <- igraph::as_data_frame(graph)
  edf$fromto <- paste0(edf$from, "_", edf$to)
  
  edf_notde <- edf[!edf[[edge_select]],,drop=FALSE]
  edf_de <- edf[edf[[edge_select]],,drop=FALSE]
  
  # edges that are not DE
  edf_notde_uni <- edf_notde[!duplicated(edf_notde$fromto) & !edf_notde$fromto %in% edf_de$fromto,]
  
  if (nrow(edf_notde) == 0) return(graph)
  
  if (!is.null(edf_notde_uni$label)) edf_notde_uni$label <- ""
  
  # merge effects/types
  
  if (!is.null(edf_notde_uni$PSITE)){
    edf_notde_uni$PSITE <- sapply(edf_notde_uni$fromto, function(tmpe){
      psites <- unique(edf_notde$PSITE[edf_notde$fromto == tmpe])
      if (length(psites) > 1) return(paste0(psites, collapse = "/")) else return(psites)
    })
    
    edf_notde_uni$PSITE_NAME <- paste0(edf_notde_uni$to, "_", edf_notde_uni$PSITE)
    
  }
  
  if (!is.null(edf_notde_uni$EFFECT)){
    edf_notde_uni$EFFECT <- sapply(edf_notde_uni$fromto, function(tmpe){
      effs <- unique(edf_notde$EFFECT[edf_notde$fromto == tmpe])
      if (length(effs) > 1) return("unknown") else return(effs)
    })
  }
  
  
  
  if (!is.null(edf_notde_uni$MECHANISM)){
    edf_notde_uni$MECHANISM <- sapply(edf_notde_uni$fromto, function(tmpe){
      effs <- unique(edf_notde$MECHANISM[edf_notde$fromto == tmpe])
      if (length(effs) > 1) return("other") else return(effs)
    })
    
  }
  
  
  
  # edges that are DE
  edf_new <- dplyr::full_join(edf_notde_uni, edf_de)
  graph <- graph_from_data_frame(edf_new, vertices = igraph::as_data_frame(graph, what = "vertices"))
  
  return(graph)
}




de_edges <- function(IG, edgepadj){
  
  edf <- igraph::as_data_frame(IG)
  edf$fromto <- paste0(edf$from, "_", edf$to)
  
  # edges that are not DE
  edf$DE <- edf$fromto %in% edf$fromto[naf(edf$padj <= 0.05)]
  edf_notde <- subset(edf, !DE)
  edf_notde_uni <- edf_notde[!duplicated(edf_notde$fromto),]
  edf_notde_uni$label <- ""
  # merge effects/types
  
  edf_notde_uni$EFFECT <- sapply(edf_notde_uni$fromto, function(tmpe){
    effs <- unique(edf_notde$EFFECT[edf_notde$fromto == tmpe])
    if (length(effs) > 1) return("unknown") else return(effs)
  })
  
  edf_notde_uni$MECHANISM <- sapply(edf_notde_uni$fromto, function(tmpe){
    effs <- unique(edf_notde$MECHANISM[edf_notde$fromto == tmpe])
    if (length(effs) > 1) return("other") else return(effs)
  })
  
  # edges that are DE
  edf_de <- subset(edf, DE)
  edf_de <- edf_de[naf(edf_de$padj <= 0.05),]
  IG <- graph_from_data_frame(rbind(edf_notde_uni, edf_de), vertices = igraph::as_data_frame(IG, what = "vertices"))
  
  return(IG)
}



unionNet <- function(subgraphs, graph, by = NULL){
  
  ### Function to generate subgraph union networks using a parent reference graph.
  
  if (rlang::is_formula(by)) by <- labels(terms(by))
  stopifnot(length(by) == 1)
  
  if (!is.null(rownames(graph$design)) & !is.null(names(subgraphs))) graph$design <- graph$design[names(subgraphs),]
  groups <- graph$design[[by]]
  stopifnot(length(groups) == length(subgraphs))
  
  # get subgraph node and edge ids
  ids <- lapply(setNames(unique(groups), unique(groups)), function(tmpgroup){
    vids <- lapply(subgraphs[groups == tmpgroup], function(tmpgraph) V(tmpgraph)$name )
    eids <- lapply(subgraphs[groups == tmpgroup], function(tmpgraph) E(tmpgraph)$ID )
    list("vids" = vids, "eids" = eids)
  })
  
  # get subgraph unions containing the node and edge ids
  uniong <- lapply(ids, function(tmp){
    
    edf <- igraph::as_data_frame(graph, what = "edges")
    ndf <- igraph::as_data_frame(graph, what = "vertices")
    
    edf_sub <- subset(edf, ID %in% unique(unlist(tmp$eids)))
    ndf_sub <- subset(ndf, name %in% unique(unlist(tmp$vids)))
    
    edf_orig <- as.data.frame(row.names = edf_sub$ID, lapply(tmp$eids, function(tmpe) edf_sub$ID %in% tmpe ))
    ndf_orig <- as.data.frame(row.names = ndf_sub$name, lapply(tmp$vids, function(tmpv) ndf_sub$name %in% tmpv ))
    
    tmpg <- igraph::graph_from_data_frame(edf_sub, vertices = ndf_sub)
    tmpg <- setdata(tmpg, SUBGRAPHS = edf_orig, type = "edges")
    tmpg <- setdata(tmpg, SUBGRAPHS = ndf_orig)
    
    tmpg
  })
  
  # subset data to the respective groups
  uniong2 <- lapply(setNames(names(uniong), names(uniong)), function(tmp){
    tmpg <- subset_data(graph = uniong[[tmp]], names(subgraphs)[tmp == groups])
    tmpg$design <- graph$design[names(subgraphs)[tmp == groups],]
    tmpg
  })
  
  
  uniong2
}






### NETWORK VISUALIZATION ------------------------------------------------------


ggnet <- function(graph, layout = NULL, node_mapping = aes(), label_mapping = aes(label = name), edge_mapping = aes(), arrow_mapping = aes(), edge_geom = ggraph::geom_edge_arc,
                   node_color = NULL, node_size = NULL, node_labsize = NULL, arrow_length = NULL,pie_nodesize=2,
                   edge_colour = NULL, edge_labsize = NULL, edge_width = NULL, edge_endgap = NULL, edge_startgap = NULL, edge_labdist = NULL, ...){
  
  ### Function to plot igraph/tidygraph objects with different line or arrow types or with pie charts as nodes.
  
  ## Functions ----
  update_aes <- function(mapping, ...) {
    ggplot2:::rename_aes(modifyList(mapping, ...))  
  }
  
  get_edges2 <- function (format = "short", collapse = "none", ix, ...) 
  {
    if (!collapse %in% c("none", "all", "direction")) {
      stop("Collapse must be either \"none\", \"all\" or \"direction\"")
    }
    
    function(layout) {
      
      edges <- ggraph:::collect_edges.layout_tbl_graph(layout)
      
      edges <- switch(collapse, none = edges, all = ggraph:::collapse_all_edges(edges), 
                      direction = ggraph:::collapse_dir_edges(edges))
      edges <- switch(format, short = ggraph:::format_short_edges(edges, 
                                                                  layout), long = ggraph:::format_long_edges(edges, layout), 
                      stop("Unknown format. Use either \"short\" or \"long\""))
      edges <- do.call(cbind, c(list(edges), lapply(list(...), 
                                                    rep, length.out = nrow(edges)), list(stringsAsFactors = FALSE)))
      attr(edges, "type_ggraph") <- "edge_ggraph"
      subset(edges, ix)
      
    }
  }
  
  
  ## Input data  ----
  edf <- igraph::as_data_frame(graph)
  E(graph)$fromto <- paste0(edf$from, "_", edf$to)
  graph <- as_tbl_graph(graph)
  graph %<>% activate(edges) %>% mutate(multi = E(graph)$fromto %in% E(graph)$fromto[duplicated(E(graph)$fromto)])
  args <- list(...)
  
  
  # set sizes and colors
  if (!is.null(node_mapping[["size"]])){
    node_size <- node_mapping[["size"]]
  }
  if (is.null(node_size)) node_size <- 40/sqrt(vcount(graph)+5)
  if (is.null(node_labsize)) node_labsize <- 5/log(sqrt(vcount(graph)+5))
  if (is.null(edge_startgap)) edge_startgap <- 0
  if (is.null(edge_endgap)) edge_endgap <- 2 + 25/(vcount(graph)) #3 + (node_size)/5
  if (is.null(edge_labsize)) edge_labsize <- node_labsize * 0.4
  if (is.null(edge_labdist)) edge_labdist <- unit(0.1/sqrt(ecount(graph)+2), 'npc')
  if (is.null(arrow_length)) arrow_length <- 1 + 10/(vcount(graph)) # 10/sqrt(ecount(graph)+2)
  if (is.null(edge_width)) edge_width <- 1/log(sqrt(ecount(graph)+2))
  if (is.null(edge_colour) & is.null(edge_mapping[["colour"]])) edge_colour <- rgb(0.4, 0.4, 0.4)
  if (is.null(node_color) & is.null(node_mapping[["colour"]])) node_color <- rgb(0.6, 0.6, 0.6)
  
  
  
  # layout
  if (is.null(layout)) layout <- cbind(x = V(graph)$x, y = V(graph)$y)
  if (is.null(layout)){
    layout <- layout_with_stress(graph)
  }
  if (!is.null(layout)) dimnames(layout) <- list(V(graph)$name, c("x", "y"))
  
  if ("x" %in% vertex_attr_names(graph)) graph <- graph %N>% dplyr::select(-c(x, y))
  
  
  ## Ggraph ----
  gg <- ggraph(graph, layout = layout) + 
    theme_graph(base_family = "sans") + 
    theme(text = element_text(family = "sans"), plot.title = element_text(hjust=0.5))
  
  
  ## Edges ----
  
  # edges
  edge_params <- args[["edges"]]
  if (is.null(edge_params)){
    
    edge_params <- list(start_cap = circle(edge_startgap/100, 'npc'),
                        end_cap = circle(edge_endgap/100, 'npc'),
                        angle_calc = "along",
                        label_dodge = edge_labdist,
                        label_size = edge_labsize,
                        label_colour = rgb(0.3, 0.3, 0.3))
    edge_params$colour <- edge_colour
  }
  
  edge_params2 <- edge_params
  edge_params2$end_cap <- circle(edge_endgap/100 * 1.5, 'npc')
  
  
  
  # arrows
  length_attr <- arrow_length
  type_attr <- "open"
  ends_attr <- "last"
  angle_attr <- 60
  
  if (!is.null(arrow_mapping[["length"]])) length_attr <- as.numeric(rlang::as_name(arrow_mapping[["length"]]))
  if (!is.null(arrow_mapping[["type"]])) type_attr <- rlang::as_name(arrow_mapping[["type"]])
  if (!is.null(arrow_mapping[["ends"]])) ends_attr <- rlang::as_name(arrow_mapping[["ends"]])
  if (!is.null(arrow_mapping[["angle"]])) angle_attr <- arrow_mapping[["angle"]]
  
  
  if (!is.null(arrow_mapping[["angle"]])){
    
    angle_params <- args[["angle"]]
    if (!rlang::as_name(angle_attr) %in% edge_attr_names(graph)) stop("Please add ", rlang::as_name(angle_attr), " as edge attribute to graph!")
    all_angle_attrs <- graph %E>% pull(!!angle_attr) %>% unique()
    
    if (is.null(angle_params) | !all(all_angle_attrs %in% names(angle_params))) angle_params <- all_angle_attrs %>% setNames((1:length(.)) * 90/length(.), .)
    
    arrows <- lapply(setNames(names(angle_params), names(angle_params)), function(tmp){
      arrow(length = unit(length_attr/100, 'npc'),
            angle = angle_params[tmp],
            ends = ends_attr,
            type = type_attr)
    })
  } else {
    all_angle_attrs <- "default"
    arrows <- list("default" = arrow(length = unit(length_attr/100, 'npc'),
                                     angle = angle_attr,
                                     ends = ends_attr,
                                     type = type_attr))
  }
  
  args <- args[!names(args) %in% c("edges", "angle")]
  
  # draw edges
  for (aattr in all_angle_attrs){
    
    # draw solid arrows
    if (length(all_angle_attrs) == 1) aes_tmp <- update_aes(edge_mapping, aes(linetype = NULL))
    if (length(all_angle_attrs) > 1) aes_tmp <- update_aes(edge_mapping, aes(filter = !!angle_attr == !!aattr, linetype = NULL))
    gg <- gg + do.call(edge_geom, args = c("data" = get_edges2(ix = !E(graph)$multi), "mapping" = list(aes_tmp), "arrow" = list(arrows[[aattr]]), width = edge_width, linetype = "solid", args, edge_params))
    gg <- gg + do.call(ggraph::geom_edge_fan, args = c("data" = get_edges2(ix = E(graph)$multi), "mapping" = list(aes_tmp), "arrow" = list(arrows[[aattr]]), strength = 1, width = edge_width, linetype = "solid", args[names(args) != "strength"], edge_params)) #
    gg <- gg + do.call(ggraph::geom_edge_loop, args = c("mapping" = list(aes_tmp), "arrow" = list(arrows[[aattr]]), width = edge_width, linetype = "solid", edge_params))
    
    # overlay solid lines with (shorter) white line
    if (length(all_angle_attrs) == 1) aes_tmp <- update_aes(edge_mapping, aes(linetype = NULL))
    if (length(all_angle_attrs) > 1) aes_tmp <- update_aes(edge_mapping, aes(filter = !!angle_attr == !!aattr, linetype = NULL))
    gg <- gg + do.call(edge_geom, args = c("data" = get_edges2(ix = !E(graph)$multi), "mapping" = list(aes_tmp), linetype = "solid", width = edge_width * 1.09, color = "white", args, edge_params2)) 
    gg <- gg + do.call(ggraph::geom_edge_fan, args = c("data" = get_edges2(ix = E(graph)$multi), "mapping" = list(aes_tmp), strength = 1, linetype = "solid", width = edge_width * 1.09, color = "white", args[names(args) != "strength"], edge_params2)) 
    gg <- gg + do.call(ggraph::geom_edge_loop, args = c("mapping" = list(aes_tmp), linetype = "solid", width = edge_width * 1.09, color = "white", edge_params2))
    
    # draw edges of line-type line without arrow
    if (length(all_angle_attrs) == 1) aes_tmp <- edge_mapping
    if (length(all_angle_attrs) > 1) aes_tmp <- update_aes(edge_mapping, aes(filter = !!angle_attr == !!aattr))
    gg <- gg + do.call(edge_geom, args = c("data" = get_edges2(ix = !E(graph)$multi), "mapping" = list(aes_tmp), width = edge_width, args, edge_params2))
    gg <- gg + do.call(ggraph::geom_edge_fan, args = c("data" = get_edges2(ix = E(graph)$multi), "mapping" = list(aes_tmp), width = edge_width, strength = 1, args[names(args) != "strength"], edge_params2))
    gg <- gg + do.call(ggraph::geom_edge_loop, args = c("mapping" = list(aes_tmp), width = edge_width, edge_params2))
  }
  
  
  ## Nodes ----
  
  if ("pie" %in% names(node_mapping)){
    
    
    gg <- gg + geom_node_point(color =  "lightgray",size=pie_nodesize)
    
    # if (!is.null(node_mapping[["size"]])) node_size <- node_mapping[["size"]] else node_size <- 0.15
    node_size <- node_size/60
    node_mapping <- update_aes(node_mapping, aes(x = x, y = y))
    pie_attr <- node_mapping[["pie"]]
    col_attr <- node_mapping[["fill"]]
    node_mapping <- update_aes(node_mapping, aes(pie = NULL))
    
    piegraph <- to_subgraph(graph, V(graph)$pie, subset_by = "nodes")$subgraph
    piedata <- piegraph %N>% pull(!!pie_attr)
    piedata <- t(sapply(piedata, function(tmp) tmp ))
    rownames(piedata) <- piegraph %N>% pull(name)
    
    pies <- reshape2::melt(data.matrix(piedata), varnames = c("node", "column"),  value.name = rlang::as_name(col_attr))
    pies$node <- as.character(pies$node)
    pies$column <- as.character(pies$column)
    pies$value <- 1
    pies <- cbind(pies, layout[match(pies$node, rownames(layout)),]) # adds x, y
    
    gg <- gg + geom_scatterpie(mapping = node_mapping,
                               color = NA,
                               cols = "column",
                               legend_name = "pie",
                               data = pies,
                               pie_scale = ncol(piedata),
                               long_format = T)
    
    # edit the last layer to replace the default scatterpie scale
    node_mapping <- update_aes(node_mapping, aes(x = NULL, y = NULL, x0 = x, y0 = y, r0 = 0, r = node_size, amount = 1))
    gg$layers[[length(gg$layers)]]$mapping <- node_mapping
    
    gg <- gg + coord_equal(clip = "off")
    
    
    
    
  } else {
    
    node_args <- list("mapping" = node_mapping, "size" = node_size)
    node_args$colour <- node_color
    gg <- gg + do.call(ggraph::geom_node_point, args = node_args)
    
    gg <- gg + coord_cartesian(clip = "off")
  }
  
  
  gg <- gg + ggraph::geom_node_text(mapping = label_mapping, size = node_labsize)
  
  
  
  return(gg)
}









ggcolorbar <- function(..., title = NULL, min = -1, max = 1, nticks = 5, FUN = scale_fill_gradient2, tickcol = "white", textcol = "black", res = 500, extend = 1.02, fontsize = 12){
  
  df <- data.frame(x = seq(min, max, 1/res), y = 0)
  dflab <- data.frame(x = scales::breaks_extended(n = nticks)(df$x))
  df$x <- df$x * extend
  
  gg <- ggplot(df, aes(x = x, y = y, fill = x)) +
    theme_void(base_size = fontsize, base_family = "sans") +
    geom_raster(width = 1.01/res, show.legend = FALSE) + ylim(c(-1,1)) +
    do.call(FUN, list(...)) + 
    geom_segment(data = dflab, aes(group = x, x = x, xend = x, yend = -0.5, y = -0.4),  color = tickcol) +
    geom_segment(data = dflab, aes(group = x, x = x, xend = x, yend = 0.5, y = 0.4),  color = tickcol) +
    geom_text(data = dflab, aes(x = x, y = -0.55, label = x), color = textcol, vjust = 1, size = fontsize / ggplot2::.pt)
  
  if (!is.null(title)){
    gg <- gg + annotate("text", x = 0, y = 0.55, label = title, vjust = 0, size = fontsize / ggplot2::.pt)
  }
  
  return(gg)
}


ggraph_legend <- function(gg, leg.fontsize = 8, labels = NULL, heights = NULL, ...){
  
  require("ggforce", quietly = TRUE)
  require("cowplot", quietly = TRUE)
  
  leg.edge.params <- gg$plot_env$edge_params
  leg.edge.arrow <- gg$plot_env$arrows
  leg.edge.arrow <- lapply(leg.edge.arrow, function(tmp){
    tmp$length <- unit(2/100, 'npc') 
    tmp
  })
  
  plot_pie <- !is.null(gg$plot_env$pies)
  
  theme_leg <- theme_void(base_size = leg.fontsize) + 
    theme(plot.title = element_text(size = leg.fontsize)) +
    theme(plot.margin = margin(0.04, 0.01, 0.12, 0, unit = "npc"))
  
  
  # Node and edge colorscales
  
  sc_aes <- sapply(gg$scales$scales, function(sc) sc$aesthetics)
  
  leg_edgecol <- NULL
  if ("edge_colour" %in% sc_aes){
    br <- gg$scales$scales[sc_aes == "edge_colour"][[1]]$break_info()
    brcol <- gg$scales$scales[sc_aes == "edge_colour"][[1]]$break_positions()
    names(brcol) <- br$major_source
    leg_edgecol <- ggcolorbar(low = brcol[1], mid = brcol[ceiling(length(brcol)/2)], high = brcol[length(brcol)],
                              midpoint = as.numeric(names(brcol)[ceiling(length(brcol)/2)]), min = br$range[1], max = br$range[2], fontsize = leg.fontsize) + theme_leg
  }
  
  leg_nodecol <- NULL
  if ("colour" %in% sc_aes){
    br <- gg$scales$scales[sc_aes == "colour"][[1]]$break_info()
    brcol <- gg$scales$scales[sc_aes == "colour"][[1]]$break_positions()
    names(brcol) <- br$major_source
    leg_nodecol <- ggcolorbar(low = brcol[1], mid = brcol[ceiling(length(brcol)/2)], high = brcol[length(brcol)],
                              midpoint = as.numeric(names(brcol)[ceiling(length(brcol)/2)]), min = br$range[1], max = br$range[2], fontsize = leg.fontsize) + theme_leg
  }
  
  # Pie charts
  if (plot_pie == TRUE){
    piedata <- gg$plot_env$piedata
    # df <- data.frame(fraction = sub("CRC.*_", "", colnames(piedata))) %>% mutate(value = 1)
    df <- data.frame(fraction = colnames(piedata)) %>% mutate(value = 1)
    df$start <- ((0:nrow(df))*2*pi/(nrow(df)))[1:nrow(df)]
    df$end <- ((0:nrow(df))*2*pi/(nrow(df)))[-1]
    df$mid <- rowMeans(df[,c("start", "end")])
    df$offset <- 0
    df$offset[1:round(nrow(df)/2)] <- 0.1
    df$offset[(round(nrow(df)/2)+1):nrow(df)] <- -0.1
    #df$offset[1] <- 0.1
    df$offset[nrow(df)] <- -0.3
    
    leg_pie <- ggplot(data = df) + 
      theme_leg +
      geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 0.9, start = start, end = end), fill = "grey", color = "white", show.legend = FALSE) + #fill = fraction
      geom_text(aes(x = 1.4*sin(mid), y = 1.4*cos(mid), label = fraction), color = "black", nudge_x = df$offset, size = leg.fontsize / ggplot2::.pt, show.legend = FALSE) + #color = fraction
      coord_fixed(ratio = 1, xlim = c(-1.5, 3), ylim = c(-1.5, 1.5), clip = "off")+
      ggtitle("")
    
    if ("fill" %in% sc_aes){
      br <- gg$scales$scales[sc_aes == "fill"][[1]]$break_info()
      brcol <- gg$scales$scales[sc_aes == "fill"][[1]]$break_positions()
      names(brcol) <- br$major_source
      leg_nodecol <- ggcolorbar(low = brcol[1], mid = brcol[ceiling(length(brcol)/2)], high = brcol[length(brcol)],
                                midpoint = as.numeric(names(brcol)[ceiling(length(brcol)/2)]), min = br$range[1], max = br$range[2], fontsize = leg.fontsize) + theme_leg
    }
    
  }
  
  
  
  # Linetypes and arrows
  elt <- gg$scales$scales[sc_aes == "edge_linetype"][[1]]$palette(0)
  df1 <- data.frame(name = names(elt), par = elt)
  df1$name <- factor(df1$name, ordered = TRUE, levels = rev(df1$name))
  leg_lty <- ggplot(df1, aes(x = 0, y = name)) +
    theme_leg + xlim(0,3) +
    geom_segment(aes(xend = rep(1, nrow(df1)), yend = name, linetype = par)) +
    scale_linetype_identity() +
    scale_y_discrete(position = "right") + 
    geom_text(aes(label = name, x = 1.05), hjust = 0, size = leg.fontsize / ggplot2::.pt)+
    ggtitle("")
  
  df2 <- data.frame(name = names(leg.edge.arrow), par = sapply(leg.edge.arrow, "[[", "angle"), par2 = sapply(leg.edge.arrow, "[[", "length"))
  df2$par2 <- unit(df2$par2, "npc") * 10
  df2$name <- factor(df2$name, ordered = TRUE, levels = rev(df2$name))
  leg_eff <- ggplot(df2, aes(x = 0, y = name)) +
    theme_leg + xlim(0,3) +
    geom_segment(aes(xend = rep(1, nrow(df2)), yend = name), arrow = arrow(angle = df2$par, length = df2$par2)) +
    scale_y_discrete(position = "right") + 
    geom_text(aes(label = name, x = 1.05), hjust = 0, size = leg.fontsize / ggplot2::.pt)+
    ggtitle("")
  
  
  # Combine
  
  if (plot_pie == TRUE){
    if (is.null(heights)) heights <- c(0.8, 0.6, 0.6, 1, 1/12, 1/12, 0.8)
    if (is.null(labels)) labels <- c("phosphosite LFC", "kinase NES", "treatment", "mechanism", "effect")
    heights[5] <- heights[5]*nrow(df1) + 0.25
    heights[6] <- heights[6]*nrow(df2) + 0.25
    leg_comb <- plot_grid(NULL, leg_edgecol, leg_nodecol, leg_pie, leg_lty, leg_eff, NULL,
                          labels = c("", labels, ""),
                          hjust = 0, vjust = 1, label_fontfamily = "sans", label_size = leg.fontsize*1, label_fontface = "bold",
                          ncol = 1, rel_heights = heights, axis = "l")
  } else {
    if (is.null(heights)) heights <- c(0.8, 0.6, 0.6, 1/12, 1/12, 1.4)
    if (is.null(labels)) labels <- c("phosphosite LFC", "kinase NES", "mechanism", "effect")
    heights[4] <- heights[4]*nrow(df1) + 0.25
    heights[5] <- heights[5]*nrow(df2) + 0.25
    leg_comb <- plot_grid(NULL, leg_edgecol, leg_nodecol, leg_lty, leg_eff, NULL,
                          labels = c("", labels, ""),
                          hjust = 0, vjust = 1, label_fontfamily = "sans", label_size = leg.fontsize*1, label_fontface = "bold",
                          ncol = 1, rel_heights = heights, axis = "l")
  }
  
  
  return(leg_comb)
}








addLayout <- function(IG, layout){
  
  if (is.null(dim(layout))) return(layout)
  
  if (is.null(colnames(layout))){
    colnames(layout) <- c("x", "y")
  }
  
  if(!is.null(rownames(layout))){
    layout <- layout[V(IG)$name,]
  }
  
  if (nrow(layout) == length(V(IG))){
    V(IG)$x <- layout[,"x"]
    V(IG)$y <- layout[,"y"]
  } else {
    message("Error: Layout has wrong number of nodes!")
  }
  
  return(IG)
}


getLayout <- function(IG){
  
  layout <- cbind(V(IG)$x, V(IG)$y)
  if (is.null(layout)){
    layout <- matrix(nrow = vcount(IG), ncol = 2)
  }
  colnames(layout) <- c("x", "y")
  rownames(layout) <- V(IG)$name
  
  return(layout)
}


layout_remove_overlaps <- function(layout, label_width = NULL, label_heigth = NULL, maxiter = 30, doplot = FALSE){
  
  # Function to adjust network layout so that labels do not overlap
  # Set label width and height as fraction of overall plot width and height
  
  
  #### Functions
  
  getOverlap <- function(a, b){
    
    # x overlap
    startpos <- max(a["x1"], b["x1"])
    endpos <- min(a["x2"], b["x2"])
    dx <- endpos - startpos
    dx[dx<0] <- 0 # no overlap
    
    # y overlap
    startpos <- max(a["y1"], b["y1"])
    endpos <- min(a["y2"], b["y2"])
    dy <- endpos - startpos
    dy[dy<0] <- 0 # no overlap
    
    overlap_area <- dx * dy
    return(overlap_area)
  }
  
  
  #### Main
  
  if (doplot == TRUE){
    require("circlize", quietly = TRUE)
    cols <- rand_color(nrow(layout))
    print(plot(layout, pch = 20, axes = FALSE, col = cols, xlab = "", ylab = ""))
  }
  
  
  if (is.null(label_heigth) | is.null(label_width)){
    dim.layout <- c("x" = diff(range(layout[,1])), "y" = diff(range(layout[,2]))) # dimensions of layout
    area_per_node <- prod(dim.layout) * 0.04 / nrow(layout)
    aspect_ratio <- 2
    x <- sqrt(area_per_node/aspect_ratio)
    label_heigth <- x
    label_width <- aspect_ratio*x
  }
  
  
  
  orig.layout <- layout
  count <- 0
  overlaps <- Inf
  
  repeat {
    
    count <- count+1
    
    # get base areas
    dim.layout <- c("x" = diff(range(layout[,1])), "y" = diff(range(layout[,2]))) # dimensions of layout
    dim.lab <- c( dim.layout["x"]*label_width, dim.layout["y"]*label_heigth) # label width and height fractions converted to unit layout
    node.areas <- cbind(x1 = layout[,1] - dim.lab[1]/2, x2 = layout[,1] + dim.lab[1]/2, y1 = layout[,2] - dim.lab[2]/2, y2 = layout[,2] + dim.lab[2]/2) # area (x- and y-range) occupied by each vertex
    
    all.areas <- sum(abs(node.areas[,1] - node.areas[,2]) * abs(node.areas[,3] - node.areas[,4]))
    if (all.areas > prod(dim.layout)){
      message("Warning: Not enough space!")
    }
    
    
    # calculate overlaps for all node pairs (using the previously calculated label positions)
    overlaps <- apply(node.areas, 1, function(tmppos1){
      apply(node.areas, 1, function(tmppos2){getOverlap(tmppos1, tmppos2)})
    })
    diag(overlaps) <- 0
    
    if ( naf(sum(overlaps != 0) == 0) | naf(count > maxiter)){
      message("Finished.")
      break }
    
    message(paste0("Iteration ", count, "..."))
    
    # get all overlapping node labels
    id <- t(sapply(1:nrow(overlaps), function(tmp1) {
      sapply(1:ncol(overlaps), function(tmp2){
        c(tmp1, tmp2)
      }, simplify = F)
    }))
    overlapping_pairs <- id[ overlaps != 0 ]
    
    # loop over all overlapping pairs
    max_overlap <- prod(dim.lab)
    for (i in 1:length(overlapping_pairs)){
      
      ## for each pair, move the overlapping pairs apart from each other and update their positions:
      
      tmppair <- overlapping_pairs[[i]]
      
      dev_orig_nodeA <- layout[tmppair[1],] - orig.layout[tmppair[1],]
      dev_orig_nodeB <- layout[tmppair[2],] - orig.layout[tmppair[2],]
      
      force <- 1 + (overlaps[tmppair[1],tmppair[2]] / max_overlap)^2 # make force dependent on overlap size
      backforce_nodeA <- 1/(1 + (dim.layout*0.05 + dev_orig_nodeA/dim.layout)^-4) # in [0,1]
      backforce_nodeB <- 1/(1 + (dim.layout*0.05 + dev_orig_nodeB/dim.layout)^-4) # in [0,1]
      
      oldpos_nodeA <- layout[tmppair[1],]
      oldpos_nodeB <- layout[tmppair[2],]
      
      x0 <- mean(c(oldpos_nodeA[1], oldpos_nodeB[1])) # x-center
      y0 <- mean(c(oldpos_nodeA[2], oldpos_nodeB[2])) # y-center
      
      # apply force
      newpos_nodeA <- c(x0 + (oldpos_nodeA[1] - x0)*force, y0 + (oldpos_nodeA[2] - y0)*force)
      newpos_nodeB <- c(x0 + (oldpos_nodeB[1] - x0)*force, y0 + (oldpos_nodeB[2] - y0)*force)
      
      # apply backforce as linear combination of new and original positions
      newpos_nodeA <- rowSums(cbind("new" = newpos_nodeA, "orig" = orig.layout[tmppair[1],]) * cbind("new" = 1-backforce_nodeA, "old" = backforce_nodeA))
      newpos_nodeB <- rowSums(cbind("new" = newpos_nodeB, "orig" = orig.layout[tmppair[2],]) * cbind("new" = 1-backforce_nodeB, "old" = backforce_nodeB))
      
      layout[tmppair[1],] <- newpos_nodeA
      layout[tmppair[2],] <- newpos_nodeB
      
    }
    
    if (doplot == TRUE){
      print(plot(layout, pch = 20, axes = FALSE, col = cols, xlab = "", ylab = ""))
    }
    
    
  }
  
  return(layout)
}


vseqEdges <- function(vseq, graph, use_weights = TRUE, use_score = TRUE){
  
  ### Function to select and convert sequences of nodes into edge IDs
  
  if (is.null(E(graph)$ID)) stop("Error: Please provide unique edge IDs!")
  if (is.null(E(graph)$weight)) E(graph)$weight <- 1
  if (is.null(E(graph)$SCORE)) E(graph)$SCORE <- 1
  edf <- igraph::as_data_frame(graph)
  
  edf$SCORE[is.na(edf$SCORE)] <- 0
  edf$weight[is.na(edf$weight)] <- max(edf$weight, na.rm = TRUE)
  
  # update loops
  vseq_loops <- vseq[sapply(vseq, length) == 1]
  vseq[sapply(vseq, length) == 1] <- lapply(vseq_loops, function(tmp) c(tmp, tmp))
  
  # get edges
  
  # For each consecutive pair of nodes, get all edges between them
  eseq <- lapply(vseq, function(tmpseq){
    v_steps <- setNames(1:(length(tmpseq)-1), paste0(tmpseq[-length(tmpseq)], "_", tmpseq[-1]))
    e_steps <- lapply(v_steps, function(i){
      subset(edf, (from == tmpseq[i] & to == tmpseq[i+1]) | (from == tmpseq[i+1] & to == tmpseq[i]) )$ID
    })
    return(e_steps)
  })
  
  # eseq:from_to:step1:edges
  
  allpaths <- unique(names(eseq))
  eseq_filtered <- lapply(setNames(allpaths, allpaths), function(tmp){
    
    tmppaths <- eseq[tmp == names(eseq)]
    
    if (length(tmppaths) == 1) return(tmppaths)
    
    names(tmppaths) <- paste0("path", seq(tmppaths))
    keep_paths <- setNames(rep(FALSE, length(tmppaths)), names(tmppaths))
    
    edf_list <- lapply(setNames(names(tmppaths), names(tmppaths)), function(tmppath_name){
      tmppath <- tmppaths[[which(tmppath_name == names(tmppaths))]]
      lapply(setNames(names(tmppath), names(tmppath)), function(tmpstep_name){
        tmpstep <- tmppath[[which(tmpstep_name == names(tmppath))]]
        data.frame(PATH = tmppath_name, STEP = tmpstep_name, subset(edf, ID %in% tmpstep))
      })
      
    })
    
    df <- Reduce(rbind, unlist(edf_list, recursive = FALSE))
    
    
    # keep all with any DE nodes or edges
    keep_paths <- keep_paths | sapply(names(tmppaths), function(tmppath) any(subset(df, PATH == tmppath)$DE) )
    
    # select path with lowest weight
    if (!any(keep_paths) & use_weights == TRUE){
      sum_weight <- sapply(names(tmppaths), function(tmppath) sum( aggregate(weight ~ STEP, subset(df, PATH == tmppath), FUN = min)[,"weight"] ) )
      keep_paths <- keep_paths | sum_weight == min(sum_weight)
    }
    
    # select path with highest average annotation score
    if (!any(keep_paths) & use_score == TRUE){
      p_scores <- sapply(names(tmppaths), function(tmppath) mean( aggregate(SCORE ~ STEP, subset(df, PATH == tmppath), FUN = max)[,"SCORE"] ) )
      keep_paths <- keep_paths | p_scores == max(p_scores)
    }
    
    
    # if none selected, return all
    if (!any(keep_paths)) keep_paths <- keep_paths | TRUE
    
  })
  
  
  
  ### add directions????
  
  
  edges <- unique(unlist(eseq_filtered))
  return(edges)
  
}


matchdf <- function(df, ids, by = NULL, ...){
  if (is.null(by)) x <- rownames(df) else x <- df[[by]]
  df <- df[match(ids, x, ...),]
  if (any(duplicated(ids))){ df <- data.matrix(df); message("Warning: Duplicate rownames!!!") }
  rownames(df) <- ids
  return(df)
}


collapseSets <- function(setlist, i){
  
  which_dups <- sapply(setlist, function(tmp1) sapply(setlist, function(tmp2) {setequal(tmp1[[i]], tmp2[[i]]) } ))
  diag(which_dups) <- FALSE
  dup_names <- sapply(colnames(which_dups), function(tmp1) rownames(which_dups)[which_dups[,tmp1]] )
  dedup_names <- sapply(names(dup_names), function(tmp){
    if (length(dup_names[[tmp]]) > 1){
      paste(sort(c(tmp, dup_names[[tmp]])), collapse = "/")
    } else {
      tmp
    }
  })
  
  names(setlist) <- dedup_names
  return(setlist[unique(dedup_names)])
}




writeTables <- function(dflist, file, rownames = TRUE){
  
  if (is.null(names(rownames)) & length(rownames) > 1) names(rownames) <- names(dflist)
  
  wb <- openxlsx::createWorkbook()
  invisible(  lapply(names(dflist), function(tmpname){
    openxlsx::addWorksheet(wb, tmpname)
    openxlsx::writeData(wb, sheet = tmpname, x = as.data.frame(dflist[[tmpname]]), rowNames = rownames)
  }))
  openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)
}


ptmsigdb2kinase <- function(PTMSEAresults, ids, match = "KINASE|SIGNOR"){
  
  stopifnot(all.equal(rownames(PTMSEAresults$NES), rownames(PTMSEAresults$pval)))
  stopifnot(all.equal(rownames(PTMSEAresults$NES), rownames(PTMSEAresults$FDR)))
  stopifnot(all.equal(rownames(PTMSEAresults$NES), PTMSEAresults$signatures$signatures))
  
  all <- rownames(PTMSEAresults$NES)
  PTMSEAresults$kinases <- data.frame(orig = all)
  
  all_split <- setNames(strsplit(all, split = ")/", fixed = TRUE), all)
  all_is_kinase <- lapply(all_split, grepl, pattern = match)
  all_kinases <- mapply(all_split, all_is_kinase, FUN = function(x,y) x[y] )
  all_kinases <- lapply(all_kinases, function(tmp) tmp %>%
                          gsub(")", "", .) %>%
                          gsub(" (SIGNOR", "_(SIGNOR)", ., fixed = TRUE) %>%
                          gsub(" (PTMSigDB", "_(PTMSigDB)", ., fixed = TRUE) %>%
                          gsub("KINASE-PSP_", "", .))
  
  
  id_kinases <- lapply(all_kinases, function(tmp){
    if (length(tmp) == 0) return(NA)
    tmp <- tmp %>% gsub("_(.*)", "", .) %>%
      strsplit(., split = "/", fixed = TRUE) %>% unlist() %>%
      strsplit(., split = "|", fixed = TRUE) %>% unlist()
    ids[tolower(ids) %in% tolower(tmp)]
  })
  
  db_kinases <- sapply(all_kinases, function(tmp){
    if (length(tmp) == 0) return(NA)
    tmpr <- sapply(tmp, function(tmptmp) paste0("SIGNOR"[grepl("SIGNOR", tmptmp)], "PTMSigDB"[grepl("PTMSigDB", tmptmp)]) )
    paste0(tmpr, collapse = "|")
  })
  
  PTMSEAresults$kinases$all <- setNames(sapply(all_kinases, paste0, collapse = "|"), NULL)
  PTMSEAresults$kinases$ID <- setNames(sapply(id_kinases, paste0, collapse = "|"), NULL)
  PTMSEAresults$kinases$ID[is.na(PTMSEAresults$kinases$ID) | PTMSEAresults$kinases$ID == "NA"] <- ""
  PTMSEAresults$kinases$DB <- db_kinases
  
  PTMSEAresults$kinases$n <- strsplit(PTMSEAresults$signatures$sites, split = "|", fixed = T) %>% sapply(., length)
  # PTMSEAresults$kinases <- PTMSEAresults$kinases[order(rowSums(is.na(PTMSEAresults$NES))),]
  PTMSEAresults$kinases <- PTMSEAresults$kinases[order(PTMSEAresults$kinases$n, decreasing = TRUE),]
  rownames(PTMSEAresults$kinases) <- PTMSEAresults$kinases$orig
  PTMSEAresults <- lapply(PTMSEAresults, function(tmp) tmp[PTMSEAresults$kinases$orig,] )
  
  PTMSEAresults$kinases$ID[duplicated(PTMSEAresults$kinases$ID)] <- ""
  kinnames <- PTMSEAresults$kinases$ID
  PTMSEAresults <- lapply(PTMSEAresults, function(tmp) tmp[nchar(kinnames) > 0,] )
  
  stopifnot(all.equal(PTMSEAresults$kinases$orig, rownames(PTMSEAresults$NES)))
  rownames(PTMSEAresults$NES) <- PTMSEAresults$kinases$ID
  rownames(PTMSEAresults$FDR) <- PTMSEAresults$kinases$ID
  rownames(PTMSEAresults$pval) <- PTMSEAresults$kinases$ID
  rownames(PTMSEAresults$signatures) <- PTMSEAresults$kinases$ID
  
  
  
  return(PTMSEAresults)
}



### UTILS ----------------------------------------------------------------------

### pipes ----

`%L>%` <- function(lhs, rhs){
  lhs <- substitute(lhs)
  rhs <- substitute(rhs)
  pe <- parent.frame()
  res <- listpipe(lhs = lhs, rhs = rhs, env = pe)
  return(res)
}

`%<L>%` <- function(lhs, rhs){
  lhs <- substitute(lhs)
  rhs <- substitute(rhs)
  pe <- parent.frame()
  res <- invisible(listpipe(lhs = lhs, rhs = rhs, env = pe))
  assign(deparse(lhs), res, envir = pe)
}


`%S>%` <- function(lhs, rhs){
  lhs <- substitute(lhs)
  rhs <- substitute(rhs)
  pe <- parent.frame()
  res <- listpipe(lhs = lhs, rhs = rhs, env = pe, simplify = TRUE)
  return(res)
}

`%<S>%` <- function(lhs, rhs){
  lhs <- substitute(lhs)
  rhs <- substitute(rhs)
  pe <- parent.frame()
  res <- invisible(listpipe(lhs = lhs, rhs = rhs, env = pe, simplify = TRUE))
  assign(deparse(lhs), res, envir = pe)
}


`%P>%` <- function(lhs, rhs){
  lhs <- substitute(lhs)
  rhs <- substitute(rhs)
  pe <- parent.frame()
  res <- listpipe(lhs = lhs, rhs = rhs, env = pe, parallel = TRUE)
  return(res)
}

`%<P>%` <- function(lhs, rhs){
  lhs <- substitute(lhs)
  rhs <- substitute(rhs)
  pe <- parent.frame()
  res <- invisible(listpipe(lhs = lhs, rhs = rhs, env = pe, parallel = TRUE))
  assign(deparse(lhs), res, envir = pe)
}


### list pipe function ----


listpipe <- function(lhs, rhs, env, parallel = FALSE, simplify = FALSE, qsub = FALSE){
  
  ### LHS (list) ----
  
  lhs <- eval(lhs, envir = env)
  if (is.null(names(lhs)) & (is.character(lhs) | is.integer(lhs))) names(lhs) <- as.character(lhs)
  
  ### RHS (function) ----
  
  FUN <- rlang::call_standardise(rhs)
  args <- rlang::call_args(FUN)
  
  # args are different for an anonymous function
  anon_f <- any(unlist(sapply(args, grepl, pattern = "function")))
  if (anon_f){
    FUN <- args[[length(args)]]
    FUN <- as.call(parse(text = as.character(FUN)))
    args <- unlist(args[-(length(args)-0:1)], recursive = FALSE)
  }
  
  # check if dot is specified in function args
  if (length(args) > 0){
    
    names(args)[names(args)==""] <- NA
    dix <- which(sapply(args, function(tmp) all(tmp == ".") ))
    
    # if dot is named, use name(s) later
    dot_args <- names(args)[dix]
    dot_args <- dot_args[!is.na(dot_args)]
    if (length(dot_args) == 0) dot_args <- NULL
    
    # if dot is unnamed, remove
    drm <- dix[is.na(names(args)[dix]) | is.null(names(args)[dix])]
    if (any(drm)) args <- args[-drm]
    if (length(args) == 0) args <- NULL
    
  } else dot_args <- NULL
  
  # remove args from function call (unless they are unnamed and passed to ...)
  args <- args[!is.na(names(args))]
  rm_args <- lapply(setNames(names(args), names(args)), function(tmp) rlang::zap() )
  FUN <- rlang::call_modify(FUN, ... = rlang::zap())
  FUN <- rlang::call_modify(FUN, !!!rm_args)
  
  # remove dots from call
  ftxt <- as.character(FUN)
  ftxt <- ftxt[ftxt != "."]
  FUN <- as.call(parse(text = ftxt))
  
  # remove dot from args
  if (length(dot_args) != 0) args[dot_args] <- NULL
  
  
  ### Iterate function over list ----
  
  itfun <- function(tmp, dot_args, args, FUN, env){
    # add tmp as dot or unnamed arg
    tmp_args <- list(tmp)
    if (!is.null(dot_args)) tmp_args <- rep(tmp_args, length(dot_args))
    all_args <- c(setNames(tmp_args, dot_args), args) # update current arguments
    FUN <- rlang::call_modify(FUN, !!!all_args) # add current arguments
    eval(FUN, envir = env)
  }
  
  
  
  
  if (parallel == TRUE){
    # parallel setup...
    if(BiocParallel::bpparam()$workers > 10) BiocParallel::register(BiocParallel::MulticoreParam(workers = 2))
    res <- BiocParallel::bplapply(X = lhs, FUN = function(tmp){
      itfun(tmp, dot_args, args, FUN, env)
    }, BPPARAM = BiocParallel::bpparam())
    
  } else if (simplify == TRUE){
    res <- sapply(lhs, function(tmp){
      itfun(tmp, dot_args, args, FUN, env)
    })
    
  } else if (qsub == TRUE){
    res <- qapply(lhs, env = env, function(tmp){
      itfun(tmp, dot_args, args, FUN, env)
    })
    
  } else {
    res <- lapply(lhs, function(tmp){
      itfun(tmp, dot_args, args, FUN, env)
    })
    
  }
  
  
  return(res)
}

