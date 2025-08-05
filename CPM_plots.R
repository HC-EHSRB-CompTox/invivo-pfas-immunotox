
cpm_files <- list.files( "C:/Users/Admin/Documents/PFAS_invivo_mouse_study/DEGlists", recursive = TRUE)
cpm_files <- cpm_files[grepl( "_CPM",cpm_files)]

cpm_files_list <- list()

n <- 1

while(n<=length(cpm_files)){
  cpm_files_list[[n]] <- read_tsv(paste0("C:/Users/Admin/Documents/PFAS_invivo_mouse_study/DEGlists/",cpm_files[n]))
  n <- n + 1
  }

d28_f <- as.data.frame(cpm_files_list[1])%>%
mutate(group = "d28_f")
d28_m <- as.data.frame(cpm_files_list[2])%>%
mutate(group = "d28_m")

d28_f <- as.data.frame(cpm_files_list[3]) %>%
  mutate(group = "d28_f")
d28_m <- as.data.frame(cpm_files_list[4])%>%
  mutate(group = "d28_m")
d56_f <- as.data.frame(cpm_files_list[5])%>%
  mutate(group = "d56_f")
d56_m <- as.data.frame(cpm_files_list[6])%>%
  mutate(group = "d56_m")

fem <- full_join(d28_f, d56_f, by = "genes", keep = FALSE)
male <- full_join(d28_m, d56_m, by = "genes")

pattern <- paste0(goi_list$Gene_Symbol, collapse = "|")

d28_f <- d28_f %>%
#  filter(grepl(pattern, Gene_Symbol)) %>%
  select(Gene_Symbol, c(6:ncol(d28_f)))

d28_m <- d28_m %>%
#  filter(grepl(pattern, Gene_Symbol)) %>%
  select(Gene_Symbol, c(6:ncol(d28_m)))

d30_f <- d30_f %>%
#  filter(grepl(pattern, Gene_Symbol)) %>%
  select(Gene_Symbol, c(6:ncol(d30_f)))

d30_m <- d30_m %>%
#  filter(grepl(pattern, Gene_Symbol)) %>%
  select(Gene_Symbol, c(6:ncol(d30_m)))

d60_f <- d60_f %>%
#  filter(grepl(pattern, Gene_Symbol)) %>%
  select(Gene_Symbol, c(6:ncol(d60_f)))

d60_m <- d60_m %>%
#  filter(grepl(pattern, Gene_Symbol)) %>%
  select(Gene_Symbol, c(6:ncol(d60_m)))

#########################
#transpose

dat_list <- list(d28_f, d28_m) #, d56_f, d56_m)
names(dat_list) <- c("d28_f", "d28_m") #, "d56_f", "d56_m")
list_names <- names(dat_list)

for(i in list_names){
  
  genes <- dat_list[[i]]$Gene_Symbol
  cpm_dat <- dat_list[[i]]
  cpm_dat <- cpm_dat[,-c(1)]
  cpm_dat[is.na(cpm_dat)] <- 0
  
  cpm_dat_t <- t(cpm_dat)
  colnames(cpm_dat_t) <- genes
  
  nzv <- nearZeroVar(cpm_dat_t)
  if(is.integer(nzv)){
   dat_filtered <- cpm_dat_t
   }else{
      dat_filtered <- cpm_dat_t[, -nzv]
      }
  
  #PCA of samples
  scaled_data <- scale(dat_filtered)
  pca_result <- prcomp(scaled_data, center = TRUE, scale. = TRUE)
  
  pca_df <- as.data.frame(pca_result$x) %>%
    mutate(Sample = rownames(as.data.frame(pca_result$x))) %>%
    separate(Sample, into = c("ID", "Sex", "Treatment", "Dose", "Date"), sep = "_", remove = TRUE) %>%
    mutate(Dose = ifelse(Dose == "DSMO", "VC", Dose)) %>%
    mutate(Dose = ifelse(Dose == "DMSO", "VC", Dose)) %>%
    mutate(Dose = ifelse(Dose == "Oct2024", "VC", Dose)) %>%
    mutate(Dose = ifelse(Treatment == "VC", "VC", Dose)) %>%
    mutate(Treatment = ifelse(Treatment == "Vehicle", "VC", Treatment)) %>%
    mutate(Dose = ifelse(Dose == "D1", "0.166", Dose)) %>%
    mutate(Dose = ifelse(Dose == "D2", "0.5", Dose)) %>%
    mutate(Dose = ifelse((i == "d30_f"|i == "d30_m") & Dose == "D3","1.5", Dose)) %>%
    mutate(Dose = ifelse((i == "d60_f"|i == "d60_m") & Dose == "D3","1", Dose)) %>%
    mutate(Dose = ifelse(Dose == "D4", "1.5", Dose))
  
  # Basic PCA scatter plot (PC1 vs PC2)
  CPM_PCA <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Treatment, fill = Treatment, label = Dose)) +
    geom_point(size = 3, alpha = 0.5) +
    scale_color_manual(values = plasma(4)) +
    scale_fill_manual(values = viridis(4)) +
    xlab(paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "% variance)")) +
    ylab(paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "% variance)")) +
    geom_text_repel(size = 4,
                    max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                    show.legend = FALSE) +
    theme_minimal()
  
  CPM_PCA
  
  ggsave(filename = paste0(i, "_PCA.jpeg"), CPM_PCA, height = 5, width = 7)

}
#########################
clean_join_result <- function(df) {
  x_cols <- grep("\\.x$", names(df), value = TRUE)
  
  for (x_col in x_cols) {
    base_name <- sub("\\.x$", "", x_col)
    y_col <- paste0(base_name, ".y")
    
    if (y_col %in% names(df)) {
      df[[base_name]] <- coalesce(df[[x_col]], df[[y_col]])
      df[[x_col]] <- NULL
      df[[y_col]] <- NULL
    }
  }
  return(df)
}


fem <- full_join(d28_f, d56_f, by = "genes")
  
fem <- clean_join_result(fem) %>%
  select(-c(genes, description, Ensembl_Gene_ID, "...1"))

fem[is.na(fem)] <- 0

male <- full_join(d28_m, d56_m, by = "genes")

male <- clean_join_result(male) %>%
  select(-c(genes, description, Ensembl_Gene_ID, "...1"))

male[is.na(male)] <- 0


all_dat <- full_join(fem, male, by = "Gene_Symbol")
all_dat <- clean_join_result(all_dat)
all_dat[is.na(all_dat)] <- 0

PFOS_dat <- all_dat %>%
  select(Gene_Symbol, contains("PFOS"), contains("VC"), contains("Vehicle"))

PFOA_dat <- all_dat %>%
  select(Gene_Symbol, contains("PFOA"), contains("VC"), contains("Vehicle"))

genes <- PFOA_dat$Gene_Symbol

cpm_dat <- PFOA_dat %>%
  select(-c(Gene_Symbol))

cpm_dat_t <- t(cpm_dat)
colnames(cpm_dat_t) <- genes

nzv <- nearZeroVar(cpm_dat_t)
if(!is.integer(nzv)){
  dat_filtered <- cpm_dat_t
}else{
  dat_filtered <- cpm_dat_t[, -nzv]
}

#PCA of samples
scaled_data <- scale(dat_filtered)
pca_result <- prcomp(scaled_data, center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(pca_result$x) %>%
  mutate(Sample = rownames(as.data.frame(pca_result$x))) %>%
  separate(Sample, into = c("ID", "Sex", "Treatment", "Dose", "Date"), sep = "_", remove = TRUE) %>%
  mutate(Date = ifelse(Dose == "Oct2024", "Oct2024", Date)) %>%
  mutate(Dose = ifelse(Dose == "Oct2024", "VC", Dose)) %>%
  mutate(Dose = ifelse(Dose == "DSMO", "VC", Dose)) %>%
  mutate(Dose = ifelse(Dose == "DMSO", "VC", Dose)) %>%
  mutate(Dose = ifelse(Treatment == "VC", "VC", Dose)) %>%
  mutate(Treatment = ifelse(Treatment == "Vehicle", "VC", Treatment)) %>%
  mutate(Dose = ifelse(Dose == "D1", "0.166", Dose)) %>%
  mutate(Dose = ifelse(Dose == "D2", "0.5", Dose)) %>%
  mutate(Dose = ifelse(Date == "Oct2024" & Dose == "D3","1.5", Dose)) %>%
  mutate(Dose = ifelse(Date == "Sept2024" & Dose == "D3","1", Dose)) %>%
  mutate(Dose = ifelse(Dose == "D4", "1.5", Dose)) %>%
  mutate(Duration = ifelse(Date == "Oct2024","28-day", "56-day"))


# Basic PCA scatter plot (PC1 vs PC2)
Dose_colors_56 <- c("VC" = "#fbd524", "0.166" = "#d5546e", "0.5" =   "#8606A6", "1" = "#6300A7", "1.5" = "#0D0887")
Dose_colors_28 <- c("VC" = "#fbd524", "1.5" = "#0D0887")
pca_df$fill_color <- ifelse(pca_df$Duration == "28-day",Dose_colors_28[as.character(pca_df$Dose)] , Dose_colors_56[as.character(pca_df$Dose)])

CPM_PCA <- ggplot(pca_df, aes(x = PC1, y = PC2,
                              colour = Duration,
                              fill = Dose,
                              #label = Dose,
                              shape = Sex)) +
  geom_point(size = 3, alpha = 0.7, stroke = 1.5) +
  scale_shape_manual(name = "Sex", 
                     values = c("Female" = 21, "Male" = 24)
                     #labels = c("28-day (hollow)", "56-day (filled)")
  ) +
  
  
  scale_fill_manual(name = "Dose",
                    values = c("VC" = "#fbd524",
                               "0.166" = "#d5546e",
                               "0.5" =   "#8606A6",
                               "1" = "#6300A7",
                               "1.5" = "#0D0887"
                    )) +  
  
  scale_color_manual(name = "Duration",
                     values = c("28-day" = "#f68d45", "56-day" = "#7e03a8")) +
  
  xlab(paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "% variance)")) +
  ylab(paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "% variance)")) +
  
  geom_text_repel(aes(label = Dose),
                  colour = pca_df$fill_color,
                  size = 4,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  show.legend = FALSE) +
  theme_minimal() +
  theme(legend.position = "none")

CPM_PCA

ggsave(filename = paste0("PFOA_PCA.jpeg"), CPM_PCA, height = 5, width = 7)



#########################
dat_list <- list(d30_f, d30_m, d60_f, d60_m)
names(dat_list) <- c("d30_f", "d30_m", "d60_f", "d60_m")
list_names <- names(dat_list)

dat_list <- lapply(dat_list, unique)

for(i in list_names){
  
  cpm_dat <- dat_list[[i]] %>%
    filter(grepl("ppar", Gene_Symbol, ignore.case = TRUE))
  
  genes <- cpm_dat$Gene_Symbol
  
  cpm_dat_long <- pivot_longer(
    cpm_dat,
    cols = -Gene_Symbol,       # All columns except Gene_Symbol
    names_to = "Sample",       # New column for sample names
    values_to = "CPM"          # New column for CPM values
  )
  
  cpm_dat_long <- cpm_dat_long %>%
    separate(Sample, into = c("ID", "Sex", "Treatment", "Dose", "Date"), sep = "_", remove = TRUE) %>%
    mutate(Dose = ifelse(Dose == "DSMO", "VC", Dose)) %>%
    mutate(Dose = ifelse(Dose == "DMSO", "VC", Dose)) %>%
    mutate(Dose = ifelse(Dose == "Oct2024", "VC", Dose)) %>%
    mutate(Dose = ifelse(Treatment == "VC", "VC", Dose)) %>%
    mutate(Treatment = ifelse(Treatment == "Vehicle", "VC", Treatment)) %>%
    mutate(Dose = ifelse(Dose == "D1", "0.166", Dose)) %>%
    mutate(Dose = ifelse(Dose == "D2", "0.5", Dose)) %>%
    mutate(Dose = ifelse((i == "d30_f"|i == "d30_m") & Dose == "D3","1.5", Dose)) %>%
    mutate(Dose = ifelse((i == "d60_f"|i == "d60_m") & Dose == "D3","1", Dose)) %>%
    mutate(Dose = ifelse(Dose == "D4", "1.5", Dose))

CPM_dot <-  ggplot(cpm_dat_long, aes(x = Gene_Symbol, y = Dose)) +
    geom_point(aes(size = CPM, colour = Dose)) +
    scale_size(range = c(1, 10)) +
    scale_colour_manual(values = plasma(5)) +
    theme(axis.title = element_text(size = 20, vjust = 1.3),
          axis.text.y = element_text(size = 13, colour = "black"),
          axis.text.x = element_text(size = 13, angle = 45, hjust = 1, colour = "black"),
          panel.background = element_blank(),     # Remove plot background
          plot.background = element_blank(),      # Remove outer background
          panel.grid = element_blank()
          ) +
    labs(
      x = "Gene", y = "Treatment",
      color = "Dose (mg/kg/day)", size = "CPM"
    ) +
    if(length(unique(cpm_dat_long$Dose))>2){facet_wrap(vars(Treatment))}
  
  ggsave(filename = paste0(i, "_Dotplot_CPM.jpeg"), CPM_dot, height = 5, width = 7)
  
}



df_long <- pivot_longer(
  d28_m,
  cols = -Gene_Symbol,       # All columns except Gene_Symbol
  names_to = "Sample",       # New column for sample names
  values_to = "CPM"          # New column for CPM values
)

df_long <- df_long %>% 
  separate(Sample, into = c("ID", "Sex", "Treatment", "Dose", "Date"), sep = "_", remove = TRUE) %>%
  select(Gene_Symbol, ID, Sex, Treatment, Dose, CPM)

df_long <- df_long %>%
  mutate(Dose = ifelse(Dose == "DSMO", "DMSO", Dose)) %>%
  mutate(Treatment = ifelse(Treatment == "Vehicle", "VC", Treatment)) %>%
  mutate(Dose = ifelse(Dose == "Oct2024", "DMSO", Dose)) %>%
  mutate(Duration = "28")

d28_m <- df_long


d28_f_long <- d28_f
d28_m_long <- d28_m
d30_f_long <- d30_f
d30_m_long <- d30_m
d60_f_long <- d60_f
d60_m_long <- d60_m

all_dat <- rbind(d28_f, d28_m, d30_f, d30_m, d60_f, d60_m) %>%
  filter(CPM > 5)

all_dat <- rbind(d30_f, d28_f)

ggplot(all_dat, aes(x = Gene_Symbol, y = CPM, colour = Duration, label = Dose)) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  facet_grid(cols=vars(Treatment)) +  # Rows = Sex, Columns = Treatment
  theme_bw() +
  xlab("Gene") +
  ylab("CPM") + 
  geom_text_repel(size = 3,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 10),
                  show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
