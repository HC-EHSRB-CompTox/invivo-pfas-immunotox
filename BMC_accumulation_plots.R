library(readxl)
library(dplyr)
library(stats)
library(tidyverse)
library(stringr)
library(purrr)
require(caret)
library(data.table)
library(viridis)
library(ggrepel)

bmd_56d <- read.csv("C:/Users/Admin/Documents/PFAS_invivo_mouse_study/BMDExpress output/Combined Analyses_filtered_56day.csv", skip = 8)

df <- bmd_56d

analysis <- data.frame(analysis = df$Analysis)

df$Analysis <- analysis %>%
        separate(analysis,
               into = c("chemical", "tool", "input", "duration", "model", "sex", "method", "pvalue", "filtering", "foldfilter", "analysis"),
               sep = "_") %>%
        tidyr::unite("unique_id", c(chemical, duration, sex), sep = "_") %>%
        pull(unique_id)

group <- df$Analysis

df <- df %>%
  separate(Analysis,
           into = c("chemical", "duration", "sex"),
           sep = "_")

df$sex <- ifelse(df$sex == "F", "Female", "Male")

df$group <- group

df <- df %>%
  dplyr::select(group, c(1:6), contains("Best")) %>%
  group_by(chemical, duration, sex) %>%
  arrange(Best.BMD, .by_group = TRUE) %>%
  mutate(bmd_no = row_number())

df <- df %>%
  dplyr::select(c(1:4), Genes.Symbols, Best.BMD, Best.BMDL, Best.BMDU, bmd_no, Best.adverseDirection)

colnames(df)[c(5:8,10)] <- c("Gene", "BMD", "BMDL", "BMDU", "Direction")

df$Direction <- str_replace(df$Direction, "-1", "down")
df$Direction <- str_replace(df$Direction, "1", "up")

#df used in KEGG 
####################

#Top 400 BMD
df_400_BMD <- df %>%
  filter(bmd_no<=400)

write.csv(df_400_BMD,  paste0(Sys.Date(),"_400 BMDs.csv"))

write.csv(df,  paste0(Sys.Date(),"_all BMDs.csv"))


df_25 <- df %>%
  filter(bmd_no == 25)

df_25 <- df_25 %>%
  select(c(1:4), Gene, BMD, BMDL, BMDU)

#All BMD
plot <- ggplot(df, aes(x = BMD, y = bmd_no, colour = chemical)) +
  geom_point(size = 3, alpha= 0.7) +
  scale_x_log10() +
  scale_colour_manual(values = c("PFOA" = "#330597", "PFOS" = "#b12a90" )) +
  ylab("Gene accumulation") +
  xlab("BMD (mg/kg/day)") +
  theme_minimal() +
  theme(strip.text = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        axis.text = element_text(size= 15, colour = "black"),
        axis.title.y = element_text(size= 20, vjust = 3, colour = "black"),
        axis.title.x = element_text(size= 20, vjust = -0.7, colour = "black")
        ) +
  facet_wrap(vars(sex))

plot

#Save all BMD plot
ggsave(paste0(Sys.Date(), "_Gene accumulation_all.jpeg"), plot, height = 10, width = 10)

#Top 25 BMDs
plot_25 <- ggplot(df_25_BMD, aes(x = BMD, y = bmd_no, colour = chemical, fill = chemical, shape = Direction, label = Gene)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_shape_manual(values = c("up" = 24, "down" = 25)) +
  #scale_fill_manual(values = c("Cytokines"="#D55E00", "Chemokine Receptors"="#0072B2", "Chemokines"="#CC79A7", "Cytokine Receptors"="#009E73")) +
  #scale_color_manual(values = c("Cytokines"="#D55E00", "Chemokine Receptors"="#0072B2", "Chemokines"="#CC79A7", "Cytokine Receptors"="#009E73")) +  
  #geom_errorbar(aes(xmin = BMDL, xmax = BMDU), width = 5) +
  geom_text_repel(size = 6,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 10),
                  show.legend = FALSE) +
  scale_x_log10() +
  scale_color_manual(values = c("PFOA" = "#330597", "PFOS" = "#b12a90" )) +
  scale_fill_manual(values = c("PFOA" = "#330597", "PFOS" = "#b12a90" )) +
  ylab("BMD ranking") +
  xlab("BMD (mg/kg/day)") +
  theme_minimal() +
  theme(strip.text = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size= 20, vjust = 3),
        axis.title.x = element_text(size= 20, vjust = -0.7),
        axis.text = element_text(size = 20, colour = "black")
        ) +
  facet_wrap(vars(sex))

plot_25

ggsave(paste0(Sys.Date(), "_Gene accumulation_lowest25.jpeg"), plot_25, height = 12, width = 10)

#25th BMD
plot_25 <- ggplot(df_25, aes(x = BMD, y = chemical, colour = chemical, label = `Gene/Pathway`)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = BMDL, xmax = BMDU), width = 0.05) +
  scale_x_log10() +
  ylab("Treatment group") +
  xlab("BMD (mg/kg/day)") +
  geom_text_repel(size = 6, show.legend = FALSE) +
  theme_bw() +
  theme(strip.text = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.title = element_text(size= 20),
        axis.text = element_text(size = 15)) +
  facet_wrap(vars(sex))

plot_25

ggsave("25th gene BMDs.jpeg", plot_25, height = 10, width = 15)

############################################

df_25$type <- "25th Gene BMD"
lowest_BMD$type <- "Lowest pathway BMD"

BMD_combined <- rbind(df_25, lowest_BMD)

All_BMD <- ggplot(BMD_combined, aes(x = BMD, y = type, colour = chemical, label = `Gene/Pathway`)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = BMDL, xmax = BMDU), width = 0.05) +
  scale_x_log10() +
  ylab("BMD type") +
  xlab("BMD (mg/kg/day)") +
  geom_text_repel(size = 4, show.legend = FALSE) +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.title = element_text(size= 15),
        axis.text = element_text(size = 10)) +
  facet_grid(cols = vars(duration), rows = vars(sex))

All_BMD

ggsave("Lowest pathway and 25th gene BMDs.jpeg", All_BMD, height = 10, width = 15)

###########################################
#Top 50 and 25 BMDs by Sex
ggplot(df_50[df_50$sex == "F",], aes(x = Best.BMD, y = bmd_no, colour = chemical)) +
  geom_point() +
  scale_x_log10() +
  ylab("Accumulation") +
  xlab("BMD (mg/kg/day)") +
  geom_text(aes(label = Genes.Symbols), size = 3, vjust = -1, hjust = 0.5) +
  facet_wrap(vars(duration))

ggplot(df_50[df_50$sex == "M",], aes(x = Best.BMD, y = bmd_no, colour = chemical)) +
  geom_point() +
  scale_x_log10() +
  ylab("Accumulation") +
  xlab("BMD (mg/kg/day)") +
  geom_text(aes(label = Genes.Symbols), size = 3, vjust = -1, hjust = 0.5) +
  facet_wrap(vars(duration))

ggplot(df_25[df_25$sex == "F",], aes(x = Best.BMD, y = bmd_no, colour = chemical)) +
  geom_point() +
  scale_x_log10() +
  ylab("Accumulation") +
  xlab("BMD (mg/kg/day)") +
  geom_text(aes(label = Genes.Symbols), size = 3, vjust = -1, hjust = 0.5) +
  facet_wrap(vars(duration))

ggplot(df_25[df_25$sex == "M",], aes(x = Best.BMD, y = bmd_no, colour = chemical)) +
  geom_point() +
  scale_x_log10() +
  ylab("Accumulation") +
  xlab("BMD (mg/kg/day)") +
  geom_text(aes(label = Genes.Symbols), size = 3, vjust = -1, hjust = 0.5) +
  facet_wrap(vars(duration))

plot <- ggplot(df_25_BMD, aes(x = Best.BMD, y = bmd_no, colour = chemical, label = Genes.Symbols)) +
  geom_point() +
  scale_x_log10() +
  ylab("Accumulation") +
  xlab("BMD (mg/kg/day)") +
  theme_bw() +
  geom_text_repel(size = 5,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 15),
                  show.legend = FALSE) +
  theme(strip.text = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.title = element_text(size= 20),
        axis.text = element_text(size = 15)) +
  facet_grid(rows = vars(sex), cols = vars(duration))

plot

ggsave("Lowest 25 BMDs.jpeg", plot, height = 15, width = 20)

##################################################
#BMDs of specified genes of interest
df <- df %>% ungroup()

df$duration <- as.numeric(gsub("[^0-9]", "", df$duration))

df$Direction <- as.character(df$Direction)

df$Direction <- str_replace(df$Direction, "-1", "down")
df$Direction <- str_replace(df$Direction, "1", "up")

#colnames(df) <- gsub("Best.adverseDirection", "Direction", colnames(df) )

###grepl input
pattern <- paste0(path_genes, collapse = "|")

#pattern <- "cyp"

GOI_bmd <- df %>%
  filter(chemical == "PFOS" & sex == "F") %>%
  filter(grepl(pattern, Gene))
  

#"cdk|p53|gadd|53b|rad|parp|xrcc|chtf|atad"

plot_GOI <- ggplot(GOI_bmd, aes(x = BMD, y = bmd_no, shape = chemical, label = Gene, colour = Direction)) +
  geom_point() +
  scale_x_log10() +
  scale_colour_manual(values=plasma(3)) +
  ylab("BMD ranking") +
  xlab("BMD (mg/kg/day)") +
  theme_bw() +
  geom_text_repel(size = 3,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                  show.legend = FALSE) +
  theme(strip.text = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.title = element_text(size= 20),
        axis.text = element_text(size = 15)) +
  facet_grid(rows = vars(sex), cols = vars(duration))

plot_GOI

ggsave("GOI_BMD_PFOS_F_immune.jpeg", plot_GOI, height = 10, width = 10)


plot_GOI <- ggplot(GOI_bmd, aes(x = Best.BMD, y = Genes.Symbols, shape = chemical, label = Genes.Symbols, colour = Direction)) +
  geom_point() +
  scale_x_log10() +
  ylab("Genes") +
  xlab("BMD (mg/kg/day)") +
  theme_bw() +
  geom_text_repel(size = 3,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                  show.legend = FALSE) +
  theme(strip.text = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.title = element_text(size= 20),
        axis.text = element_text(size = 9)) +
  facet_wrap(vars(duration))

plot_GOI

