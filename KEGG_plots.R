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
      
      # Enrichment map
      emap_plot <- emapplot(pairwise_termsim(kegg_enrich))
      ggsave(emap_plot, filename = paste0("Enrichment map_",
                                          a,"_",
                                          b,"_",
                                          c,"_",
                                          #d,
                                          ".jpg"), height = 8, width = 8)
  
      }
}


#Trancriptomics
##KEGG Pawthway
files <- list.files(path = "../invivo-pfas-immunotox/Eunnara stuff/KEGG_results/", pattern = "*.csv", full.names = TRUE)

read_with_group <- function(file) {
  df <- read_csv(file)
  df$group <- tools::file_path_sans_ext(basename(file))  
  return(df)
}


all_kegg <- lapply(files, read_with_group) %>%
  bind_rows()

all_kegg <- all_kegg %>%
  mutate(qvalue = as.numeric(qvalue),
         Count = as.numeric(Count),
         Description = str_trunc(Description, 50))
group_order <- c(
  "KEGG_results_PFOA_Female_28-day",
  "KEGG_results_PFOS_Female_28-day",
  "KEGG_results_PFOA_Male_28-day",
  "KEGG_results_PFOS_Male_28-day",
  "KEGG_results_PFOA_Female_56-day",
  "KEGG_results_PFOS_Female_56-day",
  "KEGG_results_PFOA_Male_56-day",
  "KEGG_results_PFOS_Male_56-day"
)
rownames(all_kegg$group)
top_kegg <- all_kegg %>%
  group_by(group) %>%
  arrange(qvalue) %>%
  slice_head(n = 10) %>%
  ungroup()

top_kegg$group <- factor(top_kegg$group, levels = group_order)

top_pathways <- top_kegg %>%
  distinct(Description) %>%
  pull(Description)

plot_kegg <- all_kegg %>%
  filter(Description %in% top_pathways, qvalue < 0.05)

plot_kegg <- plot_kegg %>%
  left_join(
    top_kegg %>% 
      dplyr::select(group, Description) %>%
      mutate(in_top10 = TRUE),
    by = c("group", "Description")
  ) %>%
  mutate(in_top10 = ifelse(is.na(in_top10), FALSE, TRUE))

plot_kegg <- plot_kegg %>%
  mutate(logq = -log10(qvalue + 1e-10))  
plot_kegg$group <- factor(plot_kegg$group, levels = group_order)
kegg_plot <- ggplot(plot_kegg, aes(x = group, y = fct_reorder(Description, Count))) +
  geom_point(aes(size = Count, color = logq, shape = in_top10)) +
  
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c("Other significant", "Top 10 in group"),
    name = "Pathway status"
  ) +
  scale_color_viridis_c(
    option = "plasma",
    name = expression(-log[10](qvalue))
  ) +
  scale_size_continuous(name = "Gene Count") +
  labs(
    x = "Group",
    y = "KEGG Pathway",
    title = "KEGG Pathways: Top 10 and Other Significant Occurrences"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

print(kegg_plot)


ggsave(filename = "../invivo-pfas-immunotox/R input 56d/plots/kegg_top10_vs_all.png",
       plot = kegg_plot, width = 12, height = 6, dpi = 300, bg = "white")



##KEGG similarities
all_kegg <- all_kegg %>%
  extract(group, into = c("prefix", "PFAS", "Sex", "Time"),
          regex = "^(KEGG_results)_(\\w+)_(Female|Male)_(\\d+-day)$",
          remove = FALSE) %>%
  mutate(Time = gsub("-day", "", Time))

kegg_filtered <- all_kegg %>%
  filter(qvalue < 0.05)

pfas_paths <- kegg_filtered %>%
  group_by(PFAS) %>%
  summarise(pathways = list(distinct(cur_data(), Description, qvalue, Count)),
            .groups = "drop")

sex_paths <- kegg_filtered %>%
  group_by(Sex) %>%
  summarise(pathways = list(distinct(cur_data(), Description, qvalue, Count)),
            .groups = "drop")

time_paths <- kegg_filtered %>%
  group_by(Time) %>%
  summarise(pathways = list(distinct(cur_data(), Description, qvalue, Count)),
            .groups = "drop")
plot_euler_sets <- function(sets_df, title, fill_colors, outfile) {
  group_names <- sets_df[[1]]
  
  path_sets <- setNames(
    lapply(sets_df$pathways, function(df) unique(df$Description)),
    group_names
  )
  
  fit <- euler(path_sets)
  
  png(outfile, width = 5, height = 5, units = "in", res = 300)
  print(
    plot(fit,
         fills = list(fill = fill_colors, alpha = 0.5),
         labels = F,
         edges = TRUE,
         quantities = F,
         main = title)
  )
  dev.off()
  
  return(fit)
}

plot_euler_sets(pfas_paths, "DEG Overlap by PFAS Chemical",
                c("#8405a7", "#f68f44"),
                "../invivo-pfas-immunotox/R input 56d/plots/pfas_kegg.png")

plot_euler_sets(sex_paths, "DEG Overlap by Sex",
                c("#8405a7", "#f68f44"),
                "../invivo-pfas-immunotox/R input 56d/plots/sex_kegg.png")

plot_euler_sets(time_paths, "DEG Overlap by Time",
                c("#8405a7", "#f68f44"),
                "../invivo-pfas-immunotox/R input 56d/plots/time_kegg.png")

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

pairwise_jaccard <- function(paths_df, factor_name) {
  sets <- setNames(lapply(paths_df$pathways, function(df) df$Description), paths_df[[factor_name]])
  combos <- combn(names(sets), 2, simplify = FALSE)
  
  map_dfr(combos, function(combo) {
    data.frame(
      Factor = factor_name,
      Group1 = combo[1],
      Group2 = combo[2],
      Jaccard = jaccard(sets[[combo[1]]], sets[[combo[2]]])
    )
  })
}


jaccard_pfas <- pairwise_jaccard(pfas_paths, "PFAS")
jaccard_sex <- pairwise_jaccard(sex_paths, "Sex")
jaccard_time <- pairwise_jaccard(time_paths, "Time")


print(jaccard_pfas)
print(jaccard_sex)
print(jaccard_time)

get_pathway_table <- function(sets_df, factor_name) {
  sets <- setNames(sets_df$pathways, sets_df[[factor_name]])
  if (length(sets) != 2) return(NULL)
  
  group_names <- names(sets)
  set1 <- sets[[1]]
  set2 <- sets[[2]]
  
  shared_paths <- intersect(set1$Description, set2$Description)
  unique1_paths <- setdiff(set1$Description, set2$Description)
  unique2_paths <- setdiff(set2$Description, set1$Description)
  
  get_rows <- function(df, descs) {
    df %>% filter(Description %in% descs)
  }
  
  df_unique1 <- get_rows(set1, unique1_paths) %>%
    mutate(Factor = factor_name, Group = group_names[1], Category = "Unique")
  
  df_unique2 <- get_rows(set2, unique2_paths) %>%
    mutate(Factor = factor_name, Group = group_names[2], Category = "Unique")
  
  df_shared1 <- get_rows(set1, shared_paths) %>%
    mutate(Factor = factor_name, Group = "Both", Category = "Shared")
  
  bind_rows(df_unique1, df_unique2, df_shared1) %>%
    dplyr::rename(Pathway = Description) %>%
    dplyr::select(Factor, Group, Category, Pathway, qvalue, Count)
}

pathway_table <- bind_rows(
  get_pathway_table(pfas_paths, "PFAS"),
  get_pathway_table(sex_paths, "Sex"),
  get_pathway_table(time_paths, "Time")
)

kegg_long <- kegg_filtered %>%
  dplyr::select(Description, qvalue, Count, PFAS, Sex, Time)

meta_kegg <- bind_rows(
  kegg_long %>%
    mutate(Factor = "PFAS", Group = as.character(PFAS)),
  
  kegg_long %>%
    mutate(Factor = "Sex", Group = as.character(Sex)),
  
  kegg_long %>%
    mutate(Factor = "Time", Group = as.character(Time))
)

full_pathway_table <- pathway_table %>%
  mutate(Group = as.character(Group)) %>%
  left_join(meta_kegg, by = c("Pathway" = "Description", "Group", "Factor"))

write_csv(full_pathway_table, "../invivo-pfas-immunotox/Eunnara stuff/KEGG_Pathway_Overlap.csv")

plot_top5_dotplot <- function(df, factor_name) {
  df_subset <- df %>%
    filter(Factor == factor_name, !is.na(qvalue.x), Category == "Unique") %>%
    distinct(Group, Pathway, .keep_all = TRUE) %>%
    group_by(Factor, Group) %>%
    slice_min(order_by = qvalue.x, n = 10, with_ties = FALSE) %>%
    ungroup()
  
  ggplot(df_subset, aes(x = Group, y = fct_reorder(Pathway, qvalue.x))) +
    geom_point(aes(size = Count.x, color = -log10(qvalue.x), shape = Category)) +
    scale_color_viridis_c(option = "plasma", name = expression(-log[10](qvalue))) +
    scale_shape_manual(values = c("Unique" = 16), name = "Pathway Type") +
    scale_size_continuous(name = "Gene Count") +
    labs(
      title = paste("Top 5 KEGG Pathways per", factor_name, "Group (Unique Only)"),
      x = factor_name,
      y = "KEGG Pathway"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
}

plot_top5_dotplot(full_pathway_table, "PFAS")
plot_top5_dotplot(full_pathway_table, "Sex")
plot_top5_dotplot(full_pathway_table, "Time")

ggsave("../invivo-pfas-immunotox/R input 56d/plots/top5_PFAS_pathways.png", plot_top5_dotplot(full_pathway_table, "PFAS"), width = 6, height = 4)
ggsave("../invivo-pfas-immunotox/R input 56d/plots/top5_Sex_pathways.png", plot_top5_dotplot(full_pathway_table, "Sex"), width = 6, height = 4)
ggsave("../invivo-pfas-immunotox/R input 56d/plots/top5_Time_pathways.png", plot_top5_dotplot(full_pathway_table, "Time"), width = 6, height = 4)

##Volcano
deg <- read.csv("../Eunnara stuff/DEGs_compiled.csv")
deg <- deg %>%
  filter(dose == 1.5)
# 1. Find the smallest non-zero padj in the whole dataset
min_nonzero <- min(deg$padj[deg$padj > 0], na.rm = TRUE)

deg2 <- deg %>%
  # Replace padj == 0 with that minimum
  mutate(padj = ifelse(padj == 0, min_nonzero, padj)) %>%
  # Compute –log10(padj) and cap for plotting
  mutate(log10p = -log10(padj)) %>%
  mutate(log10p = pmin(log10p, 50)) %>%
  # Classify regulation
  mutate(regulation = case_when(
    padj < 0.05 & linearFoldChange >  1.5  ~ "Upregulated",
    padj < 0.05 & linearFoldChange < -1.5  ~ "Downregulated",
    TRUE                                ~ "No change"
  )) %>%
  # Ensure factor levels
  mutate(
    regulation = factor(regulation, levels = c("Upregulated","Downregulated","No change")),
    sex        = factor(sex,        levels = c("Male","Female"))
  )

# Colour & shape maps
cols   <- c("Upregulated"   = "red",
            "Downregulated" = "#330597",
            "No change"     = "grey")
shapes <- c("Male" = 17, "Female" = 16)

unique_combos <- unique(deg2[, c("duration","chemical")])

for (i in seq_len(nrow(unique_combos))) {
  dur  <- unique_combos$duration[i]
  chem <- unique_combos$chemical[i]
  
  sub <- filter(deg2, duration == dur, chemical == chem)
  
  p <- ggplot(sub, aes(x = log2FoldChange, y = log10p)) +
    geom_point(aes(color = regulation, shape = sex),
               size = 2, alpha = 0.3,
               position = position_jitter(width = 0, height = 0.1)) +
    scale_color_manual(values = cols) +
    scale_shape_manual(values = shapes) +
    geom_vline(xintercept = c(-0.58,0.58), linetype="dashed") +
    geom_hline(yintercept = -log10(0.05), linetype="dashed") +
    
    # Label only the top 10 by absolute FC among padj<0.05
    geom_label_repel(
      data = sub %>%
        filter(padj < 0.05) %>%
        slice_max(order_by = abs(log2FoldChange), n = 10),
      aes(label = Gene_Symbol, color = regulation),
      fill        = "white",     # white background
      label.size  = 0.5,         # box border thickness
      size        = 4,           # text size
      max.overlaps= 20,
      box.padding = 0.3 ,  
      nudge_y =-1,
      point.padding =0.1 
    ) +
    
    labs(
      title = paste("Volcano —", chem, "@", dur),
      x     = expression(Log[2]~Fold~Change),
      y     = expression(-Log[10]~adjusted~P)
    ) +
    theme_minimal()
  ggsave(
    filename = paste0("../invivo-pfas-immunotox/Eunnara stuff/volcano_", chem, "_", dur, ".png"),
    plot     = p,
    width    = 5, height = 5, dpi = 300,
    bg = "white"
  )
  print(p)
}