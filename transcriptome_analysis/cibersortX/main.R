# Deconvolution vis
library(pheatmap)
library(svglite)

myColor <- viridis_pal(option = "plasma")(nBins_each*2+1)

Deconv_results <- read.table("./deconvoluted_bulk_KC.txt",sep = "\t",header = T)
Deconv_results$merged <- gsub("_.*","",Deconv_results$Mixture)
Deconv_results <- Deconv_results[order(Deconv_results$merged),]

rownames(Deconv_results) <- Deconv_results$Mixture

colorTheme_wt = c("#c6c6c6","#606060","#c12c38","#e0775f","#f3b694","#fce2d0")
names(colorTheme_wt) <- c("cdcdcd","cdcdhfd","hfdcdcd","hfdcdhfd","hfdhfdcd","hfdhfdhfd")
annoCol <- list(merged = colorTheme_wt)
annotation_row_df <- Deconv_results[,c("merged"),drop=F]


svglite::svgstring("./deconvoluted_bulk_KC_newVis.svg", width = 10, height = 10)
pheatmap::pheatmap(Deconv_results[,c("LSEC","hepatocytes","KC","B.cells","HSC","cholangiocytes","T.cells","DC.monocyte","proliferating.cells")],
                   scale = "none",
                   #color = myColor,
                   annotation_colors = annoCol,
                   annotation_row = annotation_row_df,
                   cellheight = 10,
                   cellwidth = 10,
                   cluster_rows = F,
                   show_rownames=F)

dev.off()
