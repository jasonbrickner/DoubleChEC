GO_plotter <- function(hit_list, outdir, fontsize, termlength){
  if(length(hit_list)>20){
    SGD_YeastGenes <- read.csv("GO SGD_YeastGenesList.csv")
    results.df <- data.frame(Systematic_Name = hit_list)
    #Create a list of gene names for GO-Term Analysis
    #Select genes with significant LFC
    gene_list <- results.df %>% left_join(SGD_YeastGenes, by = "Systematic_Name")
    #Convert annotations to SGD Primary ID for use with package
    gene_list_GO <- gene_list %>%
      pull(Primary_SGID)
    #Individual Gene Enrichment
    ego.BP <- enrichGO(gene          = gene_list_GO,
                       OrgDb         = org.Sc.sgd.db,
                       keyType       = "SGD",
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = FALSE)
    
    ego.BP@result$p.adjust <- signif(ego.BP@result$p.adjust, 0)
    head(ego.BP@result[ego.BP@result$p.adjust < 0.05, c(2,3,4,6)], 20)
    if(nrow(head(ego.BP@result[ego.BP@result$p.adjust < 0.05, c(2,3,4,6)]))>1){
      plot_length <- min(length(ego.BP@result[ego.BP@result$p.adjust < 0.05, ]$p.adjust), termlength)
      GO_bar <- barplot(ego.BP, font.size = fontsize, title=paste0("Top 10 GO terms"), label_format = 21, showCategory=termlength)
      ggsave(paste(outdir, "GO analysis plot.png", sep = "/"), GO_bar, width = 5.5, height = 9)
      message(paste0(plot_length, " GO term enrichments plotted"))
    }else{print(" No GO term enrichments")}
  
  }else{print("No GO term enrichments")}
}
