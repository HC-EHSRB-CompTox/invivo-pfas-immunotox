library('clusterProfiler')
library('enrichplot')
library(org.Mm.eg.db)
library(KEGGREST)

deg_files <- list.files("C:/Users/Admin/Documents/PFAS_invivo_mouse_study/all degs")

sig_deg <- lapply(deg_files, function(file){
  df <- read_tsv(paste0("C:/Users/Admin/Documents/PFAS_invivo_mouse_study/all degs/",file))
  
  df$contrast <- str_replace(df$contrast, "F_", "Female_")
  df$contrast <- str_replace(df$contrast, "M_", "Male_")
  df$contrast <- str_replace(df$contrast, "60day", "56-day")
  df$contrast <- str_replace(df$contrast, "30day", "28-day")
  df$contrast <- str_remove(df$contrast, "_invivo_reseq")
  df$contrast <- str_remove(df$contrast, "_invivo")
  df$contrast <- str_remove(df$contrast, "_one_dose")
  df$contrast <- str_remove(df$contrast, "TempO-Seq_2503_EcclesPFASLiver_")
  df$contrast <- str_remove(df$contrast, "-DESeq_output_ALL.txt")
  df$chemical <- str_extract(df$contrast, "PFOS|PFOA")
  df$dose <- str_extract(df$contrast, "[0-9.]+(?= vs)") 
  df$sample <- str_extract(df$contrast, "^[^ ]+(?= vs )")
  
  df <- df %>%
    separate(contrast, into = c("duration", "sex"), sep = "_") %>%
    filter(duration != "28day") %>%
    mutate(direction = case_when(
      linearFoldChange > 1.5 ~ "up",
      linearFoldChange < -1.5 ~ "down",
      TRUE ~ "no change"
    )) %>%
    mutate(significance = case_when(
      padj <= 0.05 ~ "sig",
      padj > 0.05 ~ "no sig"
    ))
  
  df
}) %>%
  do.call(rbind,.)

sig_deg <- sig_deg %>%
  filter(!is.na(sig_deg$padj) & padj < 0.05 & abs(linearFoldChange) > 1.5)


#df <- read.csv(paste0("C:/Users/Admin/Documents/PFAS_invivo_mouse_study/2025-07-02_all BMDs.csv"))

sig_deg <- df

sex <- c("Male", "Female")
chemical <- c("PFOS", "PFOA")
duration <- c("28-day", "56-day")
dose <- unique(sig_deg$dose)
direction <- c("up", "down")

#KEGG
for(a in chemical) {
  for(b in sex){
 #   for(c in duration){
      #for(d in direction){
      
      gene_symbols <- sig_deg[sig_deg$chemical == a 
                              & sig_deg$sex == b 
                              & sig_deg$duration == c 
                              #& sig_deg$direction == d
                              , ]$Gene_Symbol
      
      gene_list <- bitr(gene_symbols, fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = org.Mm.eg.db)
      
      entrez_ids <- unique(gene_list$ENTREZID)
      
      kegg_enrich <- enrichKEGG(gene         = entrez_ids,
                                organism     = 'mmu',
                                pvalueCutoff = 0.05)
      
      kegg_result <- kegg_enrich[] %>%
        mutate(group = paste0(a,"_",
                        b,"_",
                        c))
      
      write.csv(kegg_result, paste0("KEGG_results_",
                                            a,"_",
                                            b,"_",
                                            c,
                                            ".csv"))
      
      head(kegg_enrich)
      
      # Barplot
      bar_plot <- barplot(kegg_enrich, showCategory=20) +
        scale_fill_gradient(low = "#330597", high = "red") +
        theme_minimal() +
        theme(
          axis.text = element_text(colour ="black",
                                   size = 12))
      
      bar_plot
      
      ggsave(bar_plot, filename = paste0("Bar_plot_",
                                         a,"_",
                                         b,"_",
                                         c,"_",
                                         #d,
                                         ".jpg"), height = 10, width = 8)
      
      # Dotplot
      #dotplot(kegg_enrich, showCategory=10)
      
      # Enrichment map
      emap_plot <- emapplot(pairwise_termsim(kegg_enrich))
      ggsave(emap_plot, filename = paste0("Enrichment map_",
                                          a,"_",
                                          b,"_",
                                          c,"_",
                                          #d,
                                          ".jpg"), height = 8, width = 8)
      
      # cnetplot (gene-concept network)
      #cnetplot(kegg_enrich, categorySize="pvalue", foldChange=NULL)
      
      # List pathways
      #pathways <- keggList("pathway", "mmu")
   
      
      # Get pathway details
      #pathway_detail <- keggGet("mmu04110") 
      
      # Explore genes in a pathway
      #pathway_genes <- pathway_detail[[1]]$GENE
      }
    }
#  }
#}

####Venn Diagram

folder_path <- "C:/Users/Admin/Documents/PFAS_invivo_mouse_study/KEGG_results/"
file_list <- list.files(folder_path)
KEGG_files <- file_list[str_detect(file_list, "KEGG_results")]

KEGG_files_all <- lapply(KEGG_files, function(file){
    dat <- read.csv(paste0(folder_path,file))
    dat
  }) %>%
    do.call(rbind,.)

KEGG_groups <- KEGG_files_all %>%
  group_split(group)

d28_f <- KEGG_groups[[1]]

d28_m <- KEGG_groups[[2]]

d56_f <- KEGG_groups[[3]]

d56_m <- KEGG_groups[[4]]



      