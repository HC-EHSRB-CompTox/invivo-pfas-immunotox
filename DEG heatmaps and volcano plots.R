library(ComplexHeatmap)
library(pheatmap)
library(stringr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(viridis)
library(ggrepel)

#Load DEG files and identify significant degs (linear fold change > 1.5, padj <0.05)
deg_files <- list.files("C:/Users/Admin/Documents/PFAS_invivo_mouse_study/all degs/deg_rev", recursive = TRUE)
output_all_files <- str_subset(deg_files, "output_ALL\\.txt$")

#goi_list <- read_xlsx("C:/Users/Admin/Documents/PFAS_invivo_mouse_study/sex chromosome linked genes.xlsx")

sig_deg <- lapply(output_all_files, function(file){
  df <- read_tsv(paste0("C:/Users/Admin/Documents/PFAS_invivo_mouse_study/all degs/deg_rev/",file))

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
   # filter(duration != "28day") %>%
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
  filter(!is.na(sig_deg$padj))

write.csv(sig_deg, "DEGs_compiled.csv")

##########
sig_deg <- as.data.frame(sig_deg) %>%
  filter(duration== "28-day")

PFOA_f <- sig_deg %>%
  filter(sex == "Female" & chemical == "PFOA" & significance == "sig") %>%
  dplyr::select(Gene_Symbol, dose, direction, linearFoldChange)

PFOS_f <- sig_deg %>%
  filter(sex == "Female" & chemical == "PFOS" & significance == "sig") %>%
  dplyr::select(Gene_Symbol, dose, direction, linearFoldChange)

PFOA_m <- sig_deg %>%
  filter(sex == "Male" & chemical == "PFOA" & significance == "sig") %>%
  dplyr::select(Gene_Symbol, dose, direction, linearFoldChange)

PFOS_m <- sig_deg %>%
  filter(sex == "Male" & chemical == "PFOS" & significance == "sig") %>%
  dplyr::select(Gene_Symbol, dose, direction, linearFoldChange)


##########
#Flag genes with padj = 0  
sig_deg$shape_group <- ifelse(sig_deg$padj == 0, "zero_padj", "other")

#Set padj 0 to E-300
sig_deg$padj_plot <- ifelse(sig_deg$padj == 0, 1e-300, sig_deg$padj)
sig_deg$log10_padj <- -log10(sig_deg$padj_plot)

###### Cut-offs needs to be adjusted
sig_deg$log10_padj_capped <- ifelse(sig_deg$log10_padj > 100, 15, sig_deg$log10_padj)

sig_deg$capped <- ifelse(sig_deg$log10_padj > 100, "padj>0", "padj=0")


fc_thresh <- 1.5
padj_thresh <- 0.05

pval_zero <- sig_deg[sig_deg$padj == 0,]

sex <- c("Male", "Female")
chemical <- c("PFOS", "PFOA")
duration <- c("28-day", "56-day")
dose <- unique(sig_deg$dose)

#################Volcano plot and heatmap for each chemical and sex
for(a in chemical) {
  for(b in sex){
    for(c in duration){
    
      #Volcano plots  
    plot_deg <- sig_deg[sig_deg$sex == b &
                          sig_deg$duration == c & 
                          sig_deg$chemical == a
                         ,]
    
    #label_genes <- plot_deg %>%
    #  filter(abs(linearFoldChange) > fc_thresh, padj < padj_thresh)
    
    top_gene <- plot_deg %>%
      filter(significance == "sig"
             & direction != "no change"
             ) %>%
      group_by(dose, direction) %>%
      arrange(desc(abs(linearFoldChange))) %>%
      slice_head(n = 10) %>%
      ungroup()
    
    volcano_plot <-ggplot(plot_deg, aes(x = log2FoldChange, y = log10_padj_capped, col = direction, shape = capped)) +
      #coord_cartesian(ylim = c(0, 6), xlim = c(-10, 10)) +
      geom_jitter(size = if(length(unique(plot_deg$dose))>2){2}else{4}, stroke = 0, alpha = 0.3) +      
      theme_classic() +
      scale_shape_manual(
        values = c("padj=0" = 16, "padj>0" = 18),
        labels = c("padj > 0", paste("padj \u2248 0"))
        ) +
      scale_color_manual(
        values = c("up" = "red", "down" = "#330597", "no change" = "grey"),
        labels = c("Downregulated", "No Significance","Upregulated")
      ) +
      scale_y_continuous(
        name = expression(-log[10](adjusted~p))
        #####Needs to be adjusted
        breaks = c(0, 2, 4, 6, 8, 10),
        labels = c("0", "2", "4", "6", "8","Inf")
      ) +
      
      guides(shape = guide_legend(title = NULL)) +
      
      labs(title = paste(plot_deg$duration,plot_deg$chemical,"exposure:", plot_deg$sex),
           x = "Log2 Fold Change",
           y = "-Log10(padj)",
           color = "Direction"
      ) +
      
      theme(
        strip.text = element_text(
          face = "bold",     
          size = 7,         
          angle = 0,         
          color = "black"     
        ),
        strip.background = element_rect(
          fill = "white",  
          color = "black",      
          linewidth = 1
        ),
        axis.text = element_text(
          color = "black"
        )
      ) +
      
      geom_hline(yintercept = c(-log10(0.05), 9), linetype = "dashed", color = "azure4") +  # Horizontal line for p-value cutoff
      geom_vline(xintercept = 0, linetype = "dashed", color = "azure4") +        # Vertical lines for log2FC cutoff
      
      geom_text_repel(
        data = top_gene,
        aes(label = Gene_Symbol), 
        size = if(length(unique(plot_deg$dose))>2){2.5}else{4},
        max.overlaps = 30, 
        box.padding = if(length(unique(plot_deg$dose))>2){0.5}else{1.2},
        point.padding = if(length(unique(plot_deg$dose))>2){0.5}else{1.5},
        segment.color = "azure3",
        show.legend = FALSE
      ) +
      facet_wrap(vars(chemical))
  #    if(length(unique(plot_deg$dose))>2){facet_wrap(vars(chemical))}

    volcano_plot
    
    ggsave(volcano_plot, filename = paste0("Volcano_plot","_", unique(plot_deg$duration),"_", unique(plot_deg$chemical),"_", unique(plot_deg$sex),".jpg"), height = 3, width = 5)
      
    }

  ###########################################################
  #DEG heatmaps
 if(FALSE){   
  plot_deg <- sig_deg[!is.na(sig_deg$padj) & sig_deg$sex == a & sig_deg$padj <= 0.05 & sig_deg$chemical == b,]
    
  plot_deg <- plot_deg %>%
    select(Gene_Symbol, log2FoldChange, sample)
  
  dup_genes <- plot_deg %>%
    dplyr::summarise(n = dplyr::n(), .by = c(Gene_Symbol, sample)) %>%
    dplyr::filter(n > 1)
  
  plot_deg <- plot_deg %>%
    group_by(Gene_Symbol, sample) %>%
    summarise(log2FC = mean(log2FoldChange, na.rm = TRUE), .groups = "drop")
  
  wide_df <- plot_deg %>%
    pivot_wider(names_from=sample, values_from=log2FC)
  
  wide_df[is.na(wide_df)] <- 0
  
  row_genes <- wide_df$Gene_Symbol
  
  wide_df <- wide_df %>%
    dplyr::select(-c(Gene_Symbol))
  
  rownames(wide_df) <- row_genes
  
  wide_df <- as.matrix(wide_df)
  
  days <- sub("_.*", "", colnames(wide_df))
  
  doses <- as.numeric(sub(".*_", "", colnames(wide_df)))
  
  ordering_df <- data.frame(
    colname = colnames(wide_df),
    day = factor(days, levels = c("28-day", "56-day")),  # Set order explicitly
    dose = doses,
    stringsAsFactors = FALSE
  )

  ordered_indices <- with(ordering_df, order(day, dose))
  
  df_ordered <- wide_df[, ordered_indices]
  
  ht <- Heatmap(df_ordered,
                col = colorRamp2(c(-2.5,0,2.5),c("#330597", "white", "red")),
                cluster_columns = FALSE,
                #legend
                show_heatmap_legend = TRUE,
                heatmap_legend_param = list(title = "Log2 Fold Change"),
                
                #row names
                show_row_names = FALSE,
                
                #col names
                column_names_rot = 45,
                
                width = unit(ncol(df_ordered) * 1, "cm")
                )
  
  
  jpeg(paste0("Heatmap_", a, "_", b, ".jpeg"), width = 1250, height = 1000, res = 200)
  draw(ht)
  dev.off()     
        
    } }
  }

#####################################################
plot_deg <- sig_deg[!is.na(sig_deg$padj) & sig_deg$sex == "F" & sig_deg$padj <= 0.05 & sig_deg$duration != "28day" & sig_deg$chemical == "PFOA",]

plot_deg <- plot_deg %>%
  select(Gene_Symbol, log2FoldChange, sample)

dup_genes <- plot_deg %>%
    dplyr::summarise(n = dplyr::n(), .by = c(Gene_Symbol, sample)) %>%
    dplyr::filter(n > 1)

plot_deg <- plot_deg %>%
  group_by(Gene_Symbol, sample) %>%
  summarise(log2FC = mean(log2FoldChange, na.rm = TRUE), .groups = "drop")

wide_df <- plot_deg %>%
  pivot_wider(names_from=sample, values_from=log2FC)

wide_df[is.na(wide_df)] <- 0

row_genes <- wide_df$Gene_Symbol

wide_df <- wide_df %>%
  dplyr::select(-c(Gene_Symbol))

rownames(wide_df) <- row_genes

wide_df <- as.matrix(wide_df)

pheatmap(wide_df,
         scale = "row",           # Optional: standardize rows
         cluster_rows = TRUE,     # Cluster genes
         cluster_cols = TRUE,     # Cluster samples
         color = colorRampPalette(c("navy", "white", "maroon"))(100),
         show_rownames = FALSE,
         filename = paste0("Heatmap_F_PFOS.png"),
         width = 10,      # in inches
         height = 10,
         res = 300)# Set TRUE if you want to label genes

