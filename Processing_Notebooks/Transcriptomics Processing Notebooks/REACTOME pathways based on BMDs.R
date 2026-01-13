library(readxl)
library(dplyr)
library(stats)
library(tidyverse)
library(stringr)
library(purrr)
require(caret)
library(data.table)
library(ggrepel)

#Load REACTOME output from BMDExpress
sig_path <- read_tsv("C:/Users/Admin/Documents/PFAS_invivo_mouse_study/BMDExpress output/Combined Analyses_filtered_REACTOME_July2025", skip = 4)
sig_path <- read.csv("C:/Users/Admin/Documents/PFAS_invivo_mouse_study/BMDExpress output/GO_biological_processes_all.csv", row.names = NULL)

dat <- sig_path

dat$Analysis<- str_remove(dat$Analysis, "_reseq")

analysis <- data.frame(analysis = dat$Analysis)

dat$Analysis <- analysis %>%
  separate(analysis,
           into = c("chemical", "tool", "input", "duration", "model", "sex", "method", "pvalue", "filtering", "foldfilter", "analysis"),
           sep = "_") %>%
  unite("unique_id", c(chemical, duration, sex), sep = "_") %>%
  pull(unique_id)

group <- dat$Analysis

dat <- dat %>%
  separate(Analysis,
           into = c("chemical", "duration", "sex"),
           sep = "_")

dat$group <- group

dat <- dat %>%
  filter(duration == "60day")

dat <- dat %>%
  mutate(duration = ifelse(dat$duration == "60day", "56-day", NA))

write.csv(dat, "GO_data_56-day.csv")

#Apply filters to REACTOME hits
dat_f <- dat %>%
  filter(dat$`Fisher's Exact Two Tail` <= 0.1 & dat$`Genes That Passed All Filters` > 5 & dat$`Input Genes` >= 3 & dat$Percentage >= 5)

dat_f <- dat_f %>%
  group_by(chemical, duration, sex) %>%
  arrange(dat_f$'BMD Median', .by_group = TRUE) %>%
  mutate(bmd_no = row_number()) %>%
  ungroup()

dat_top <- dat_f %>%
  filter(bmd_no <= 20)

lowest_BMD <- dat_f %>%
  filter(bmd_no < 15) %>%
  select(group, c(1:3), `GO/Pathway/Gene Set/Gene Name`, `BMD Median`, `BMDL Median`, `BMDU Median`)
  
  colnames(lowest_BMD)[5:8] <- c("Gene/Pathway", "BMD", "BMDL", "BMDU")

#Pathways with the lowest BMDs
lowest_bmd_plot <- ggplot(lowest_BMD, aes(x = BMD, y = duration, colour = chemical, label = `Gene/Pathway`)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = BMDL, xmax = BMDU), width = 0.05) +
  scale_x_log10() +
  ylab("Treatment group") +
  xlab("BMD (mg/kg/day)") +
  geom_text_repel(size = 6, show.legend = FALSE) +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.title = element_text(size= 15),
        axis.text = element_text(size = 10)) +
  facet_wrap(vars(sex))

lowest_bmd_plot

lowest_BMD <- lowest_BMD %>%
  arrange(BMD, .by_group = FALSE)

lowest_bmd_plot_all <- ggplot(lowest_BMD, aes(x = BMD, y = reorder(group, BMD), colour = chemical, label = `Gene/Pathway`)) +
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
        axis.text = element_text(size = 15))
# facet_wrap(vars(sex))

lowest_bmd_plot_all

ggsave("All lowest pathway BMDs.jpeg", lowest_bmd_plot_all, height = 10, width = 10)

#########################

dat_f$duration <- gsub("[^0-9]", "", dat_f$duration)

dat_top <- dat_f %>%
  filter(bmd_no < 10, duration != "28day")

write.csv(dat_f, "Reactome results_58day.csv")
#Genes in a select pathway

#dat_top <- dat_f %>%
#  filter(grepl("phago|immun|antig|complem|neutro", dat_f$'GO/Pathway/Gene Set/Gene Name', ignore.case = TRUE))

path_genes <- dat_top %>%
  select('Gene Symbols') %>%
  tidyr::separate_rows('Gene Symbols', sep = ";") %>%
  distinct() %>%
  filter('Gene Symbols' != "") %>%
  pull('Gene Symbols')

plot <- ggplot(dat_top, aes(x = dat_top$'BMD Median', y = bmd_no, colour = chemical,label = dat_top$'GO/Pathway/Gene Set/Gene Name' )) +
  geom_point(size = 2.5, alpha = 0.6) +
  scale_x_log10() +
  scale_colour_manual(values = c("PFOA" = "#330597", "PFOS" = "#b12a90" )) +
  ylab("Pathway BMD Ranking") +
  xlab("BMD (mg/kg/day)") +
  theme_minimal() +
  geom_text_repel(size = 4.5,
                  max.overlaps = 3,
                  show.legend = FALSE,
                  box.padding = 1,
                  point.padding = 2) +
  theme(strip.text = element_text(size = 20),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        axis.title = element_text(size= 20),
        axis.text = element_text(size = 15, colour = "black"),
        strip.text.y = element_text(colour = "black")
        ) +
  #facet_wrap(vars(duration))
  facet_wrap(~sex, labeller = as_labeller(c("F" = "Female", "M" = "Male")))

plot

ggsave("Top10 pathway_BMD.jpeg", plot, height = 10, width = 13)


