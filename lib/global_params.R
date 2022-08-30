### Annotation colors

require("paletteer")
require("ggraph")

colors <- list()

colors$organoid <- as.character(paletteer::paletteer_d("ggthemes::Tableau_20", n = 16))
colors$treatment <- as.character(paletteer::paletteer_d("RColorBrewer::Dark2", n = 8))
colors$msi_status <- c("#9370DB","#FF7F24")

names(colors$organoid) <- c("CRC01", "CRC02", "CRC03", "CRC04", "CRC11", "CRC12", "CRC13", "CRC15", "CRC17", "CRC18", "CRC19", "CRC20", "CRC21", "CRC22", "CRC26", "CRC26LM")
names(colors$treatment) <- c("control", "BRAFi", "MEKi", "mTORi", "PI3Ki", "TAKi", "TBKi", "TNFalpha")
names(colors$msi_status) <- c("MSI", "MSS")
colors$treatment['TNFa'] <- colors$treatment[['TNFalpha']]


### Network colors

colors$scales <- list("phospho" = c(top = rgb(0.9,0.9,0), zero = rgb(0.93,0.93,0.93), bottom = rgb(0.5,0,0.7), lower = -4, upper = 4, name = "phosphosite LFC"),
                      "kinase" = c(top = rgb(0.8,0.5,0), zero = rgb(0.93,0.93,0.93), bottom = rgb(0,0.55,0.65), lower = -4, upper = 4, name = "kinase activity NES"))

colors$fixed <- list("basenode" = rgb(0.85,0.85,0.85),
                     "baseedge" = rgb(0.7,0.7,0.7),
                     "phospho" = rgb(0.5,0,0.9),
                     "edgepsite" = rgb(0.9,0,0),
                     "both" = rgb(t(matrix(rowMeans(cbind(col2rgb(rgb(0.9,0,0)), col2rgb(rgb(0.5,0,0.9)), ncol = 3)))), maxColorValue = 255),
                     "kinase" = rgb(0.9,0.8,0),
                     "mutation" = "#ff8400")



angle <- c("activation" = 30,
           "inhibition" = 90,
           "unknown" = 15)

ltys <- c("phosphorylation" = "solid",
          "dephosphorylation" = "twodash",
          "transcriptional regulation" = "dashed",
          "binding" = "dotted",
          "other" = "dotdash")

scale_node_fill <- scale_color_gradient2(low = "#008CA6", mid = "#EDEDED", high = "#CC8000", midpoint = 0, na.value = "#adadad",
                                         limits = c(-4,4), oob = scales::squish,
                                         guide = guide_colorbar(ticks.colour = "black"))

scale_node_fill_rna <- scale_color_gradient2(low = "blue", mid = "#EDEDED", high = "red", midpoint = 0, na.value = "#adadad",
                                         limits = c(-3,3),
                                         oob = scales::squish,
                                         guide = guide_colorbar(ticks.colour = "black"))

scale_edge_color <- scale_edge_color_gradient2(low = "#8000B3", mid = "#EDEDED", high = "#E6E600", midpoint = 0, na.value = "#adadad",
                                               limits = c(-4,4), oob = scales::squish,
                                               guide = guide_edge_colorbar(ticks.colour = "black"))

scale_pie_color <- scale_fill_gradient2(low = "#008CA6", mid = "#EDEDED", high = "#CC8000", midpoint = 0, na.value = "#919191", oob = scales::squish,
                                        limits = c(-4,4), guide = guide_colorbar(ticks.colour = "black"))



# treatment targets according to manufacturers
trt_targets <- list(BRAFi = "BRAF",
                    MEKi = c("MAP2K1", "MAP2K2"),
                    mTORi = "MTOR",
                    PI3Ki = "PIK3CB",
                    TAKi = c("MAPK1", "MAP3K7", "MAP2K7", "MAP2K1"),
                    TBKi = c("NUAK1", "TBK1", "MARK4", "AURKB", "IKBKE", "PDK1"),
                    TNFalpha = c("TNF", "TNFRSF1A", "TNFRSF1B"))
