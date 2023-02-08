# Extremophilic Fungi Metagenome Analysis
# For review paper with Quandt Mycology Lab, University of Colorado Boulder
# by Cliff Bueno de Mesquita, JGI, Summer 2022, Updated Spring 2023
# This script is based on Extreme_Fungi.R, but updated for the new set of samples
# Click the "Show document outline" button in the top right corner to view document outline
# Sections are:
# 1. Retrieving Data
# 2. Setup
# 3. Postprocessing
# 4. Taxonomic analyses
# 5. Functional analyses
# 6. Other



#### Retrieving Data ####
# Info on assembling the IMG dataset is in the Data Acquisition.text file
# Ecosystems of interest were searched by ecosystem or name
# Domain = *Microbiome
# GOLD Analysis Project Type = Metagenome Analysis
# Preprocessing yielded 1560 samples ("Extreme_Combined_Prefilt.txt")
# Postprocessing:
# -Keep Sequence center != DOE Joint Genome Institute (JGI) and Is Public = Yes
# -Keep Sequence center = DOE Joint Genome Institute (JGI) and JGI Data Utilization Status = Unrestricted
# -Filter out duplicates that were reassembled with SPAdes and contain (SPAdes) at the end of the Genome Name/ Sample Name
# -Filter out engineered ecosystems (e.g. bioreactors)
# -Check that there are no duplicate taxonoids
# -Make comma separated list of taxonoids to search for in Genome Search
# Updates:
# Quandt Lab reviewed samples - deleted some samples that were not relevant (e.g. too far from mine or glacier)
# Also found some more samples - 140 desert samples are added in this analysis



#### Setup ####
# Libraries
suppressWarnings(suppressMessages(library(readxl))) # For read_xlsx
suppressWarnings(suppressMessages(library(janitor))) # For cleaning
suppressWarnings(suppressMessages(library(cowplot))) # For multipanel
suppressWarnings(suppressMessages(library(plyr))) # For data manipulation
suppressWarnings(suppressMessages(library(tidyverse))) # For data manipulation
suppressWarnings(suppressMessages(library(reshape2))) # For melting
suppressWarnings(suppressMessages(library(vegan))) # For analysis
suppressWarnings(suppressMessages(library(car))) # For leveneTest
suppressWarnings(suppressMessages(library(PMCMRplus))) # For Nemenyi posthoc test
suppressWarnings(suppressMessages(library(indicspecies))) # For multipatt
suppressWarnings(suppressMessages(library(scales))) # For muted
suppressWarnings(suppressMessages(library(DESeq2))) # For normalization
suppressWarnings(suppressMessages(library(FSA))) # For standard error
suppressWarnings(suppressMessages(library(mctoolsr))) # For taxonomic analysis
suppressWarnings(suppressMessages(library(cowplot))) # For multipanel graphs
suppressWarnings(suppressMessages(library(plotly))) # For interactive graphs
suppressWarnings(suppressMessages(library(RColorBrewer))) # For color palettes
suppressWarnings(suppressMessages(library(dendextend))) # For dendrogram plots
suppressWarnings(suppressMessages(library(viridis))) # For viridis palette
suppressWarnings(suppressMessages(library(gplots))) # For heatmaps
suppressWarnings(suppressMessages(library(maps))) # For geographic maps
suppressWarnings(suppressMessages(library(mapproj))) # For geographic maps
suppressWarnings(suppressMessages(library(magrittr))) # For setting column names
suppressWarnings(suppressMessages(library(writexl))) # For writing Excel file
suppressWarnings(suppressMessages(library(plotrix))) # For standard error
suppressWarnings(suppressMessages(library(emmeans))) # For Tukey
suppressWarnings(suppressMessages(library(multcomp))) # For Tukey
suppressWarnings(suppressMessages(library(RCurl))) # For KEGG
suppressWarnings(suppressMessages(library(KEGGREST))) # For KEGG
suppressWarnings(suppressMessages(library(multcompView))) # For significance letters
suppressWarnings(suppressMessages(library(rcompanion))) # For significance letters
suppressWarnings(suppressMessages(library(pheatmap))) # For heatmaps
suppressWarnings(suppressMessages(library(qvalue))) # For qvalue

# Working directory (GitHub repository)
setwd("~/Documents/GitHub/Extremophilic_Fungi/")

# Functions
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
find_hullj <- function(df) df[chull(df$Axis01j, df$Axis02j),]
`%notin%` <- Negate(`%in%`)
save_pheatmap_pdf <- function(x, filename, width = 7, height = 5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

source("plot_multipatt.R")

# So far used Krusal-Wallis + Nemenyi posthoc
# May want to update some stats to zero-inflated negative binomial regression instead
# If using ANOVA and Tukey, can use this code to get sig. letters
# Example from genus richness
#tuk <- emmeans(object = m, specs = "Environment") %>%
#  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
#  mutate(name = "rich",
#         y = max(input_fungi$map_loaded$rich)+(max(input_fungi$map_loaded$rich)-min(input_fungi$map_loaded$rich))/20)


#### ..........................####
#### Taxonomic ####
# Metadata ("mapping file") downloaded from IMG

# Tax table for mctoolsr (only need to do once, then skip to Import)
#t <- read.delim("Extreme_Updated_1149/UI_data_output.txt") %>%
#  dplyr::select(-c(3:9)) %>%
#  dplyr::rename(taxonomy = FeatureName) %>%
#  dplyr::rename(ASV_ID = Feature) %>%
#  dplyr::select(ASV_ID, 3:ncol(.), taxonomy)
#names(t) <- abbreviate(names(t), minlength = 11)
#table.fp <- "~/Documents/GitHub/Extremophilic_Fungi"
#out_fp <- paste0(table.fp, "/genus_table_mctoolsr_updated.txt")
#names(t)[1] = "#ASV_ID"
#write("#Exported for mctoolsr", out_fp)
#suppressWarnings(write.table(t, out_fp, sep = "\t", row.names = FALSE, append = TRUE))



#### _Setup ####
# Import with mctoolsr (matches sampleIDs, 1141 samples)
tax_table_fp <- file.path("genus_table_mctoolsr_updated.txt")
map_fp <- file.path("metadata_updated.txt")
input = load_taxa_table(tax_table_fp, map_fp)
new_tab <- read_excel("Extremophilic_fungi_dataset_final.xlsx", sheet = 1) %>%
  dplyr::select(taxon_oid, Environment)

# Update map_loaded - sampleID, GenomeSize, Environment, Assembly Method, Year
# Filter out samples from before 2012. 1116 remaining
dim(input$map_loaded)
input$map_loaded <- input$map_loaded %>%
  mutate(sampleID = paste("X", taxon_oid, sep = ""),
         GenomeSize = `Genome Size   * assembled`) %>%
  left_join(., new_tab, by = "taxon_oid") %>%
  separate(`Add Date`, into = c("Day", "Month", "Year"), sep = "/", remove = F) %>%
  mutate(Year = as.integer(Year)) %>%
  filter(Year > 11) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(sampleID2 = sampleID) %>%
  column_to_rownames(var = "sampleID2")
dim(input$map_loaded)
input$map_loaded$`Assembly Method`[input$map_loaded$`Assembly Method` == ""] <- "Unknown"

# Filter out samples with no genus level reads (removes 7 samples, 1109 remaining)
# Note: This filters 5 unassembled samples, taxonoids 3300002080-84
count <- as.data.frame(sort(colSums(input$data_loaded))) %>%
  filter(`sort(colSums(input$data_loaded))` == 0)
input <- filter_data(input,
                     filter_cat = "sampleID",
                     filter_vals = rownames(count))

# Check sequencing depth 
sort(colSums(input$data_loaded))
mean(colSums(input$data_loaded)) # 333051.3
se(colSums(input$data_loaded)) # 16129.89
input$map_loaded$count <- colSums(input$data_loaded)
ggplot(input$map_loaded, aes(reorder(`Environment`, count, mean), count)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.25, width = 0.25) +
  labs(x = "Environment", 
       y = "# Reads") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

pdf("FigsUpdated/AssembledMetagenomeSizes_AllGenusReads.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(GenomeSize, count)) +
  geom_point(size = 1.5, alpha = 0.25) +
  geom_smooth(method = "lm") +
  labs(x = "Assembled genome size", 
       y = "Assigned genus reads") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10))
dev.off()





#### __Check Euks ####
# View Domain level taxa
# This shows that euks are not very relatively abundant in these environments
tax_sum_Domain <- summarize_taxonomy(input, level = 1, report_higher_tax = T, relative = TRUE)
plot_taxa_bars(tax_sum_Domain,
               input$map_loaded,
               type_header = "Environment",
               num_taxa = 100,
               data_only = F) +
  theme_classic() +
  labs(x = "Environment",
       y = "Relative abundance",
       fill = "Domain") +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

# Range of % abundance in individual samples
range(tax_sum_Domain[3,]) # 0 to 0.636, so can be high in some samples

# Euks
euks <- tax_sum_Domain[3,]
plot_taxa_bars(euks,
               input$map_loaded,
               type_header = "Environment",
               num_taxa = 100,
               data_only = F) +
  scale_fill_manual(values = c("grey")) +
  theme_classic() +
  labs(x = "Environment",
       y = "Relative abundance") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

euks_t <- as.data.frame(t(euks))
input$map_loaded$Euks <- euks_t$Eukaryota
ggplot(input$map_loaded, aes(reorder(sampleID, Euks, mean), Euks, fill = Environment)) +
  geom_bar(stat = "identity", color = NA) +
  labs(x = NULL, 
       y = "Relative abundance") +
  theme_classic() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Too many, filter top ones (> 5%), want to see which environments the top samples are from
eukabund <- data.frame("Euks" = input$map_loaded$Euks,
                       "sampleID" = input$map_loaded$sampleID) %>%
  filter(Euks > 0.05) %>%
  column_to_rownames(var = "sampleID")
topeuk <- filter_data(input,
                      filter_cat = "sampleID",
                      keep_vals = rownames(eukabund))
pdf("FigsUpdated/EukTopSamples.pdf", width = 7, height = 5)
ggplot(topeuk$map_loaded, aes(reorder(sampleID, Euks, mean), Euks, fill = Environment)) +
  geom_bar(stat = "identity", color = NA) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x = NULL, 
       y = "Relative abundance") +
  ggtitle("Samples with Eukaryota > 5% (n = 30)") +
  theme_classic() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
dev.off()



#### __Check Fungi ####
# First do same as done above for euk but extract fungal phyla
tax_sum_phyla <- summarize_taxonomy(input, level = 2, report_higher_tax = T, relative = TRUE)
fungal_phyla <- tax_sum_phyla[grep("Ascomycota|Basidiomycota|Blastocladiomycota|Chytridiomycota|
                  Cryptomycota|Microsporidia|Mucoromycota|Nephridiophagidae|Olpidiomycota|
                  Sanchytriomycota|Zoopagomycota", rownames(tax_sum_phyla)),]
fungal_phyla <- fungal_phyla[!grepl("Plasmid", rownames(fungal_phyla)),]
rownames(fungal_phyla) <- substring(rownames(fungal_phyla), 12)
plot_taxa_bars(fungal_phyla,
               input$map_loaded,
               type_header = "Environment",
               num_taxa = 100,
               data_only = F) +
  geom_hline(yintercept = 0.01, linetype = "dashed", color = "grey") +
  theme_classic() +
  labs(x = "Environment",
       y = "Relative abundance",
       fill = "Phylum") +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

# Get summed relative abundance of fungi and plot by sample and environment
input$map_loaded$Fungi <- colSums(fungal_phyla)
ggplot(input$map_loaded, aes(reorder(sampleID, Fungi, mean), Fungi, fill = Environment)) +
  geom_bar(stat = "identity", color = NA) +
  labs(x = NULL, 
       y = "Relative abundance") +
  theme_classic() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Too many, filter top ones (> 1%), want to see which environments the top samples are from
funabund <- data.frame("Fungi" = input$map_loaded$Fungi,
                       "sampleID" = input$map_loaded$sampleID) %>%
  filter(Fungi > 0.01) %>%
  column_to_rownames(var = "sampleID")
topfun <- filter_data(input,
                      filter_cat = "sampleID",
                      keep_vals = rownames(funabund))
pdf("FigsUpdated/FungalTopSamples.pdf", width = 7, height = 5)
ggplot(topfun$map_loaded, aes(reorder(sampleID, Fungi, mean), Fungi, fill = Environment)) +
  geom_bar(stat = "identity", color = NA) +
  geom_hline(yintercept = 0.01, linetype = "dashed") +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x = NULL, 
       y = "Relative abundance") +
  ggtitle("Samples with Fungi > 1% (n = 38)") +
  theme_classic() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
dev.off()

# By environment
leveneTest(input$map_loaded$Fungi ~ input$map_loaded$Environment) # Bad
m1 <- aov(Fungi ~ Environment, data = input$map_loaded)
shapiro.test(m1$residuals) # Bad
summary(m1)
TukeyHSD(m1)
kruskal.test(Fungi ~ Environment, data = input$map_loaded)
nyi1 <- kwAllPairsNemenyiTest(Fungi ~ Environment, data = input$map_loaded)
nyi_table1 <- fullPTable(nyi1$p.value)
nyi_list1 <- multcompLetters(nyi_table1)
nyi_let1 <- as.data.frame(nyi_list1$Letters) %>%
  mutate(label = `nyi_list1$Letters`,
         y = rep(0.12, nrow(.))) %>%
  rownames_to_column(var = "Environment")
pdf("FigsUpdated/FungalRel.pdf", width = 5, height = 4)
ggplot(input$map_loaded, aes(reorder(Environment, Fungi, mean), Fungi)) +
  #geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1, alpha = 0.2, width = 0.4) +
  geom_text(data = nyi_let1, aes(Environment, y, label = label), 
            size = 4, color = "black") +
  labs(x = NULL, y = "Fungal relative abundance") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        strip.text = element_text(size = 10))
dev.off()



#### __Filter ####
# Now filter to euk and then fungi (extract fungal phyla)
input_euk <- filter_taxa_from_input(input,
                                    taxa_to_keep = "Eukaryota",
                                    at_spec_level = 1)

# Check phyla
table(input_euk$taxonomy_loaded$taxonomy2)

# List of RefSeq fungal phyla
fungal_phyla_names <- c("Ascomycota", "Basidiomycota", "Blastocladiomycota", 
                        "Chytridiomycota","Cryptomycota", "Microsporidia", 
                        "Mucoromycota", "Nephridiophagidae", "Olpidiomycota", 
                        "Sanchytriomycota", "Zoopagomycota")
input_fungi <- filter_taxa_from_input(input_euk,
                                      taxa_to_keep = fungal_phyla_names,
                                      at_spec_level = 2)
# Filter out plasmids (2 removed)
input_fungi <- filter_taxa_from_input(input_fungi,
                                      taxa_to_remove = "Plasmid:Eukaryota",
                                      at_spec_level = 1)

nrow(input$taxonomy_loaded) # 4038 total
nrow(input_euk$taxonomy_loaded) # 400 euks
nrow(input_fungi$taxonomy_loaded) # 304 fungi

# Now check reads again
sort(colSums(input_fungi$data_loaded))
# Note lots of samples with 0 or very few fungi
# Purposefully not filtering those out those as 0's are interesting in this analysis
# These are extreme environments, some may have few to no fungi
# Further below, however, some analyses will be done with zeroes removed
mean(colSums(input_fungi$data_loaded)) # 505
se(colSums(input_fungi$data_loaded)) # 50
input_fungi$map_loaded$fung_count <- colSums(input_fungi$data_loaded)
input_fungi$map_loaded$present <- ifelse(input_fungi$map_loaded$fung_count > 0,
                                         1,
                                         0)

# Note, will do this again but with CPM normalization
ggplot(input_fungi$map_loaded, aes(reorder(`Environment`, fung_count, mean), fung_count)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.25, width = 0.25) +
  labs(x = "Environment", 
       y = "# Fungal Reads") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

# Check genome size vs fungal genus reads
pdf("FigsUpdated/AssembledMetagenomeSizes_FungalGenusReads.pdf", width = 7, height = 5)
ggplot(input_fungi$map_loaded, aes(GenomeSize, fung_count)) +
  geom_point(size = 1.5, alpha = 0.25) +
  geom_smooth(method = "lm") +
  labs(x = "Assembled genome size", 
       y = "Assigned fungal genus reads") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10))
dev.off()

# Get and plot fraction of samples with 0 versus > 0 fungal reads, and also n
env_prev <- input_fungi$map_loaded %>%
  group_by(Environment) %>%
  summarise(num_present = sum(present),
            num_samples = n(),
            prevalence = round(num_present/num_samples * 100, digits = 2)) %>%
  mutate(num_absent = num_samples - num_present)
sum(env_prev$num_present)

# Melt for stacked bar
env_prev_long <- melt(env_prev,
                      id.vars = "Environment",
                      measure.vars = c("num_present", "num_absent"))

pdf("FigsUpdated/FungalPrevalence.pdf", width = 5, height = 4)
ggplot(env_prev, aes(reorder(Environment, prevalence, mean), prevalence)) +
  geom_bar(stat = "identity") +
  labs(x = NULL, 
       y = "% prevalence of fungi") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

pdf("FigsUpdated/SampleSize.pdf", width = 5, height = 4.5)
ggplot(env_prev_long, aes(reorder(Environment, value, mean), value, 
                          group = Environment, fill = variable)) +
  geom_bar(stat = "identity") +
  geom_text(data = env_prev, 
            aes(reorder(Environment, num_samples, mean), num_samples+10, 
                label = num_samples), inherit.aes = F) +
  scale_fill_manual(values = c("#F8766D", "#619CFF"),
                    breaks = c("num_absent", "num_present"),
                    labels = c("Absent", "Present")) +
  labs(x = NULL, 
       y = "Sample size",
       fill = "Fungi") +
  ggtitle("Total sample size = 1109\nSamples with Fungi = 921") +
  theme_bw() +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_blank(),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



#### _CPM ####
input_fungi_CPM <- input_fungi
# Replace counts in "data_loaded" with CPM transformed counts
# This is CPM assembled metagenomic reads
# Do (count*1000000)/GenomeSize
# i is samples 1 to 1142
# j is taxa 1 to 304
for (i in 1:ncol(input_fungi$data_loaded)) {
  for (j in 1:nrow(input_fungi$data_loaded)) {
    input_fungi_CPM$data_loaded[j, i] <- (input_fungi$data_loaded[j, i]*1000000)/input_fungi$map_loaded$GenomeSize[i]
  }
}

# Make stacked bar plots by taxonomic level
# Resort so unassigned and other are on the top
# Use "Paired" palette from RColorBrewer
# Note: "Set2" palette could be another option
# Unassigned gets grey75, Other gets grey90
# Show all phyla; for others show top 15 - use colorRampPalette to expand colors
# Can easily update this later if people want a different number of taxa shown
ntax <- 15
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(ntax)

# Phyla
tax_sum_Phyla <- summarize_taxonomy(input_fungi_CPM, level = 2, report_higher_tax = F, relative = F)
barsPhyla <- plot_taxa_bars(tax_sum_Phyla,
                            input_fungi_CPM$map_loaded,
                            type_header = "Environment",
                            num_taxa = ntax,
                            data_only = T) %>%
  mutate(taxon = fct_rev(taxon))

pdf("FigsUpdated/CPM_Phyla.pdf", width = 7, height = 5)
ggplot(barsPhyla, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Abundance (CPM)", fill = "Phylum") +
  scale_fill_manual(values = brewer.pal(12, "Paired")[7:1]) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.5, unit = "cm"))
dev.off()
taxa_summary_by_sample_type(tax_sum_Phyla, input_fungi_CPM$map_loaded, 'Environment', 0.0001, 'KW')

# Class
tax_sum_Class <- summarize_taxonomy(input_fungi_CPM, level = 3, report_higher_tax = T, relative = F)
rownames(tax_sum_Class) <- substring(rownames(tax_sum_Class), 12)
plot_taxa_bars(tax_sum_Class,
               input_fungi_CPM$map_loaded,
               type_header = "Environment",
               num_taxa = nrow(tax_sum_Class),
               data_only = F) +
  theme_classic() +
  labs(x = "Environment", y = "Abundance (CPM)", fill = "Class") +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

tax_sum_Class <- summarize_taxonomy(input_fungi_CPM, level = 3, report_higher_tax = F, relative = F)
barsClass <- plot_taxa_bars(tax_sum_Class,
               input_fungi_CPM$map_loaded,
               type_header = "Environment",
               num_taxa = ntax,
               data_only = T) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))

pdf("FigsUpdated/CPM_Class.pdf", width = 7, height = 5)
ggplot(barsClass, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Abundance (CPM)", fill = "Class") +
  scale_fill_manual(values = c("grey90", mycolors[15:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.5, unit = "cm"))
dev.off()
taxa_summary_by_sample_type(tax_sum_Class, input_fungi_CPM$map_loaded, 'Environment', 0.0001, 'KW')

# Order
tax_sum_Order <- summarize_taxonomy(input_fungi_CPM, level = 4, report_higher_tax = F, relative = F)
barsOrder <- plot_taxa_bars(tax_sum_Order,
                            input_fungi_CPM$map_loaded,
                            type_header = "Environment",
                            num_taxa = ntax,
                            data_only = T) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))

pdf("FigsUpdated/CPM_Order.pdf", width = 7, height = 5)
ggplot(barsOrder, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Abundance (CPM)", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", mycolors[14:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.5, unit = "cm"))
dev.off()
taxa_summary_by_sample_type(tax_sum_Order, input_fungi_CPM$map_loaded, 'Environment', 0.0001, 'KW')

# Family
tax_sum_Family <- summarize_taxonomy(input_fungi_CPM, level = 5, report_higher_tax = F, relative = F)
barsFamily <- plot_taxa_bars(tax_sum_Family,
                            input_fungi_CPM$map_loaded,
                            type_header = "Environment",
                            num_taxa = ntax,
                            data_only = T) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))

pdf("FigsUpdated/CPM_Family.pdf", width = 7, height = 5)
ggplot(barsFamily, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Abundance (CPM)", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", mycolors[14:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.5, unit = "cm"))
dev.off()
taxa_summary_by_sample_type(tax_sum_Family, input_fungi_CPM$map_loaded, 'Environment', 0.0001, 'KW')

# Genus
tax_sum_Genus <- summarize_taxonomy(input_fungi_CPM, level = 6, report_higher_tax = F, relative = F)
barsGenus <- plot_taxa_bars(tax_sum_Genus,
                            input_fungi_CPM$map_loaded,
                            type_header = "Environment",
                            num_taxa = ntax,
                            data_only = T) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  # mutate(taxon = fct_relevel(taxon, "unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))

pdf("FigsUpdated/CPM_Genus.pdf", width = 7, height = 5)
ggplot(barsGenus, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Abundance (CPM)", fill = "Genus") +
  scale_fill_manual(values = c("grey90", mycolors[15:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.5, unit = "cm"))
dev.off()
taxa_summary_by_sample_type(tax_sum_Genus, input_fungi_CPM$map_loaded, 'Environment', 0.0001, 'KW')

# Also plot total fungal CPM by environment
input_fungi_CPM$map_loaded$totalFun <- colSums(input_fungi_CPM$data_loaded)
leveneTest(input_fungi_CPM$map_loaded$totalFun ~ input_fungi_CPM$map_loaded$Environment) # Bad
m2 <- aov(totalFun ~ Environment, data = input_fungi_CPM$map_loaded)
shapiro.test(m2$residuals) # Bad
summary(m2)
TukeyHSD(m2)
kruskal.test(totalFun ~ Environment, data = input_fungi_CPM$map_loaded)
nyi2 <- kwAllPairsNemenyiTest(totalFun ~ Environment, data = input_fungi_CPM$map_loaded)
nyi_table2 <- fullPTable(nyi2$p.value)
nyi_list2 <- multcompLetters(nyi_table2)
nyi_let2 <- as.data.frame(nyi_list2$Letters) %>%
  mutate(label = `nyi_list2$Letters`,
         y = rep(175, nrow(.))) %>%
  rownames_to_column(var = "Environment")
pdf("FigsUpdated/FungalCPM.pdf", width = 5, height = 4)
ggplot(input_fungi_CPM$map_loaded, aes(reorder(Environment, totalFun, mean), totalFun)) +
  # geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1, alpha = 0.2, width = 0.4) +
  geom_text(data = nyi_let2, aes(Environment, y, label = label), 
            size = 4, color = "black") +
  labs(x = NULL, y = "Fungal abundance (CPM)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        strip.text = element_text(size = 10))
dev.off()

# Could also try plot by location with environment facets
# mctoolsR only accepts one factor, so need to combine then split
input_fungi_CPM$map_loaded$EnvGeo <- paste(input_fungi_CPM$map_loaded$Environment,
                                           input_fungi_CPM$map_loaded$`Geographic Location`,
                                           sep = "_")
tax_sum_Phyla <- summarize_taxonomy(input_fungi_CPM, level = 2, report_higher_tax = F, relative = F)
barsPhyla <- plot_taxa_bars(tax_sum_Phyla,
                            input_fungi_CPM$map_loaded,
                            type_header = c("EnvGeo"),
                            num_taxa = ntax,
                            data_only = T) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Environment", "Location"), sep = "_")

pdf("FigsUpdated/CPM_Phyla_EnvGeo.pdf", width = 12, height = 8)
ggplot(barsPhyla, aes(reorder(Location, -mean_value, mean), mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Abundance (CPM)", fill = "Phylum") +
  scale_fill_manual(values = brewer.pal(12, "Paired")[7:1]) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(~ Environment, scales = "free_x", space = "free") +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 2, angle = 45, hjust = 1, vjust = 1,
                                   margin = unit(c(-1,0,0,0), "pt")),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"),
        legend.position = c(0.9, 0.5),
        strip.text = element_text(size = 8, angle = 90, hjust = 0),
        strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.5), "cm"))
dev.off()



#### _Relative ####
# Remove zeroes (188 removed, 921 remaining)
countFun <- as.data.frame(sort(colSums(input_fungi_CPM$data_loaded))) %>%
  filter(`sort(colSums(input_fungi_CPM$data_loaded))` == 0)
input_fungi_nz <- filter_data(input_fungi,
                              filter_cat = "sampleID",
                              filter_vals = rownames(countFun))

# Calculate relative abundance of fungi, only for samples with fungi (n = 921)
# Make relative abundance stacked bar plots by taxonomic level
# Re-sort so unassigned and other are on the top
# Use "Paired" palette from RColorBrewer
# Unassigned gets grey75, Other gets grey90
# Show all phyla; for others show top 15 - use colorRampPalette to expand colors
# Can easily update this later if people want a different number of taxa shown
ntax <- 15
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(ntax)

# Phyla
tax_sum_Phyla <- summarize_taxonomy(input_fungi_nz, level = 2, report_higher_tax = F, relative = T)
barsPhyla <- plot_taxa_bars(tax_sum_Phyla,
                            input_fungi_nz$map_loaded,
                            type_header = "Environment",
                            num_taxa = ntax,
                            data_only = T) %>%
  mutate(taxon = fct_rev(taxon))

pdf("FigsUpdated/Rel_Phyla.pdf", width = 7, height = 5)
ggplot(barsPhyla, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = brewer.pal(12, "Paired")[7:1]) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))
dev.off()
taxa_summary_by_sample_type(tax_sum_Phyla, input_fungi_nz$map_loaded, 'Environment', 0.0001, 'KW')

# Class
tax_sum_Class <- summarize_taxonomy(input_fungi_nz, level = 3, report_higher_tax = T, relative = T)
rownames(tax_sum_Class) <- substring(rownames(tax_sum_Class), 12)
plot_taxa_bars(tax_sum_Class,
               input_fungi_nz$map_loaded,
               type_header = "Environment",
               num_taxa = nrow(tax_sum_Class),
               data_only = F) +
  theme_classic() +
  labs(x = "Environment", y = "Relative abundance", fill = "Class") +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

tax_sum_Class <- summarize_taxonomy(input_fungi_nz, level = 3, report_higher_tax = F, relative = T)
barsClass <- plot_taxa_bars(tax_sum_Class,
                            input_fungi_nz$map_loaded,
                            type_header = "Environment",
                            num_taxa = ntax,
                            data_only = T) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%  
  mutate(taxon = fct_relevel(taxon, "unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))

pdf("FigsUpdated/Rel_Class.pdf", width = 7, height = 5)
ggplot(barsClass, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Relative abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", mycolors[14:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))
dev.off()
taxa_summary_by_sample_type(tax_sum_Class, input_fungi_nz$map_loaded, 'Environment', 0.0001, 'KW')

# Order
tax_sum_Order <- summarize_taxonomy(input_fungi_nz, level = 4, report_higher_tax = F, relative = T)
barsOrder <- plot_taxa_bars(tax_sum_Order,
                            input_fungi_nz$map_loaded,
                            type_header = "Environment",
                            num_taxa = ntax,
                            data_only = T) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  #mutate(taxon = fct_relevel(taxon, "unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))

pdf("FigsUpdated/Rel_Order.pdf", width = 7, height = 5)
ggplot(barsOrder, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Relative abundance", fill = "Order") +
  scale_fill_manual(values = c("grey90", mycolors[15:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))
dev.off()
taxa_summary_by_sample_type(tax_sum_Order, input_fungi_nz$map_loaded, 'Environment', 0.0001, 'KW')

# Family
tax_sum_Family <- summarize_taxonomy(input_fungi_nz, level = 5, report_higher_tax = F, relative = T)
barsFamily <- plot_taxa_bars(tax_sum_Family,
                             input_fungi_nz$map_loaded,
                             type_header = "Environment",
                             num_taxa = ntax,
                             data_only = T) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))

pdf("FigsUpdated/Rel_Family.pdf", width = 7, height = 5)
ggplot(barsFamily, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Relative abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", mycolors[14:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))
dev.off()
taxa_summary_by_sample_type(tax_sum_Family, input_fungi_nz$map_loaded, 'Environment', 0.0001, 'KW')

# Genus
tax_sum_Genus <- summarize_taxonomy(input_fungi_nz, level = 6, report_higher_tax = F, relative = T)
barsGenus <- plot_taxa_bars(tax_sum_Genus,
                            input_fungi_nz$map_loaded,
                            type_header = "Environment",
                            num_taxa = ntax,
                            data_only = T) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  # mutate(taxon = fct_relevel(taxon, "unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))

pdf("FigsUpdated/Rel_Genus.pdf", width = 7, height = 5)
ggplot(barsGenus, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey90", mycolors[15:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))
dev.off()
taxa_summary_by_sample_type(tax_sum_Genus, input_fungi_nz$map_loaded, 'Environment', 0.0001, 'KW')



#### _PCoAs ####
# PCoA and PERMANOVA by Environment
# Also check Habitat, Ecosystem.Category, Ecosystem.Subtype, Ecosystem.Type, Specific.Ecosystem

# Filter out fungal zeroes from CPM data (filters 188 samples, 921 remaining)
input_fungi_CPM_nz <- filter_data(input_fungi_CPM,
                                  filter_cat = "sampleID",
                                  filter_vals = rownames(countFun))

Phylum_nz <- summarize_taxonomy(input_fungi_CPM_nz, level = 2, report_higher_tax = F, relative = F)
Class_nz <- summarize_taxonomy(input_fungi_CPM_nz, level = 3, report_higher_tax = F, relative = F)
Order_nz <- summarize_taxonomy(input_fungi_CPM_nz, level = 4, report_higher_tax = F, relative = F)
Family_nz <- summarize_taxonomy(input_fungi_CPM_nz, level = 5, report_higher_tax = F, relative = F)

# Do at each taxonomic level
# Genus (lowest level)
bc <- calc_dm(input_fungi_CPM_nz$data_loaded)

# Check some PERMANOVA models to get a sense of R2 values
set.seed(1150)
adonis2(bc ~ Environment, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.088, p = 0.001
adonis2(bc ~ Habitat, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.38, p = 0.001
adonis2(bc ~ `Ecosystem Category`, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.016, p = 0.001
adonis2(bc ~ `Ecosystem Subtype`, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.138, p = 0.001
adonis2(bc ~ `Ecosystem Type`, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.064, p = 0.001
adonis2(bc ~ `Specific Ecosystem`, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.097, p = 0.001
anova(betadisper(bc, input_fungi_CPM_nz$map_loaded$Environment)) # Dispersion not homogeneous
pcoa <- cmdscale(bc, k = nrow(input_fungi_CPM_nz$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
input_fungi_CPM_nz$map_loaded$Axis01 <- scores(pcoa)[,1]
input_fungi_CPM_nz$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(input_fungi_CPM_nz$map_loaded, c("Environment"), find_hull)
pdf("FigsUpdated/PCoA_Genus.pdf", width = 7, height = 5)
g <- ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 2, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
g
dev.off()

# Interactive plot (can hover mouse over points)
ggplotly(g)

# There's basically no clustering of fungal composition by environment at genus level, lots of overlap!
# Extract legend to plot separately as its own panel later
g_leg <- get_legend(g)

# Family
bc_Family <- calc_dm(Family_nz)
set.seed(1150)
adonis2(bc_Family ~ Environment, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.096 , p = 0.001
anova(betadisper(bc_Family, input_fungi_CPM_nz$map_loaded$Environment)) # Dispersion not homogeneous
pcoa_Family <- cmdscale(bc_Family, k = nrow(input_fungi_CPM_nz$map_loaded) - 1, eig = T)
pcoaA1F <- round((eigenvals(pcoa_Family)/sum(eigenvals(pcoa_Family)))[1]*100, digits = 1)
pcoaA2F <- round((eigenvals(pcoa_Family)/sum(eigenvals(pcoa_Family)))[2]*100, digits = 1)
input_fungi_CPM_nz$map_loaded$Axis01 <- scores(pcoa_Family)[,1]
input_fungi_CPM_nz$map_loaded$Axis02 <- scores(pcoa_Family)[,2]
micro.hulls <- ddply(input_fungi_CPM_nz$map_loaded, c("Environment"), find_hull)
g1 <- ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F, size = 0.25) +
  geom_point(size = 1, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1F, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2F, "%", sep = "")) +
  ggtitle("Family") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, vjust = -1))
g1
# There's basically no clustering of fungal composition by environment at family level, lots of overlap!

# Order
bc_Order <- calc_dm(Order_nz)
set.seed(1150)
adonis2(bc_Order ~ Environment, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.105, p = 0.001 
anova(betadisper(bc_Order, input_fungi_CPM_nz$map_loaded$Environment)) # Dispersion not homogeneous
pcoa_Order <- cmdscale(bc_Order, k = nrow(input_fungi_CPM_nz$map_loaded) - 1, eig = T)
pcoaA1O <- round((eigenvals(pcoa_Order)/sum(eigenvals(pcoa_Order)))[1]*100, digits = 1)
pcoaA2O <- round((eigenvals(pcoa_Order)/sum(eigenvals(pcoa_Order)))[2]*100, digits = 1)
input_fungi_CPM_nz$map_loaded$Axis01 <- scores(pcoa_Order)[,1]
input_fungi_CPM_nz$map_loaded$Axis02 <- scores(pcoa_Order)[,2]
micro.hulls <- ddply(input_fungi_CPM_nz$map_loaded, c("Environment"), find_hull)
g2 <- ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F, size = 0.25) +
  geom_point(size = 1, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1O, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2O, "%", sep = "")) +
  ggtitle("Order") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, vjust = -1))
g2
# There's basically no clustering of fungal composition by environment at Order level, lots of overlap!

# Class
bc_Class <- calc_dm(Class_nz)
set.seed(1150)
adonis2(bc_Class ~ Environment, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.108, p = 0.001 
anova(betadisper(bc_Class, input_fungi_CPM_nz$map_loaded$Environment)) # Dispersion not homogeneous
pcoa_Class <- cmdscale(bc_Class, k = nrow(input_fungi_CPM_nz$map_loaded) - 1, eig = T)
pcoaA1C <- round((eigenvals(pcoa_Class)/sum(eigenvals(pcoa_Class)))[1]*100, digits = 1)
pcoaA2C <- round((eigenvals(pcoa_Class)/sum(eigenvals(pcoa_Class)))[2]*100, digits = 1)
input_fungi_CPM_nz$map_loaded$Axis01 <- scores(pcoa_Class)[,1]
input_fungi_CPM_nz$map_loaded$Axis02 <- scores(pcoa_Class)[,2]
micro.hulls <- ddply(input_fungi_CPM_nz$map_loaded, c("Environment"), find_hull)
g3 <- ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F, size = 0.25) +
  geom_point(size = 1, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1C, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2C, "%", sep = "")) +
  ggtitle("Class") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, vjust = -1))
g3
# There's basically no clustering of fungal composition by environment at Class level, lots of overlap!

# Phylum
bc_Phylum <- calc_dm(Phylum_nz)
set.seed(1150)
adonis2(bc_Phylum ~ Environment, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.128, p = 0.001 
anova(betadisper(bc_Phylum, input_fungi_CPM_nz$map_loaded$Environment)) # Dispersion not homogeneous
pcoa_Phylum <- cmdscale(bc_Phylum, k = nrow(input_fungi_CPM_nz$map_loaded) - 1, eig = T)
pcoaA1P <- round((eigenvals(pcoa_Phylum)/sum(eigenvals(pcoa_Phylum)))[1]*100, digits = 1)
pcoaA2P <- round((eigenvals(pcoa_Phylum)/sum(eigenvals(pcoa_Phylum)))[2]*100, digits = 1)
input_fungi_CPM_nz$map_loaded$Axis01 <- scores(pcoa_Phylum)[,1]
input_fungi_CPM_nz$map_loaded$Axis02 <- scores(pcoa_Phylum)[,2]
micro.hulls <- ddply(input_fungi_CPM_nz$map_loaded, c("Environment"), find_hull)
g4 <- ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F, size = 0.25) +
  geom_point(size = 1, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1P, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2P, "%", sep = "")) +
  ggtitle("Phylum") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, vjust = -1))
g4
# There's basically no clustering of fungal composition by environment at Phylum level, lots of overlap!
# More % variation explained at phylum level though

# Combine plot
# Remake first one with title and no legend
input_fungi_CPM_nz$map_loaded$Axis01 <- scores(pcoa)[,1]
input_fungi_CPM_nz$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(input_fungi_CPM_nz$map_loaded, c("Environment"), find_hull)
g <- ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F, size = 0.25) +
  geom_point(size = 1, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  ggtitle("Genus") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, vjust = -1))
g

# Multipanel
pdf("FigsUpdated/PCoA_AllLevels.pdf", width = 8, height = 5)
plot_grid(g,g1,g2,g3,g4,g_leg, ncol = 3, hjust = "hv")
dev.off()



#### __Prok ####
# Now that we've seen not much clustering of the fungal communities, I'm curious to see how archaeal and bacterial communities look, just for the sake of comparison.

# Archaea
input_arc <- filter_taxa_from_input(input,
                                    taxa_to_keep = "Archaea",
                                    at_spec_level = 1)
nrow(input_arc$data_loaded) # 150
input_arc_CPM <- input_arc
for (i in 1:ncol(input_arc$data_loaded)) {
  for (j in 1:nrow(input_arc$data_loaded)) {
    input_arc_CPM$data_loaded[j, i] <- (input_arc$data_loaded[j, i]*1000000)/input_arc$map_loaded$GenomeSize[i]
  }
}

# Remove zeroes (55 removed, 1054 remaining)
countArc <- as.data.frame(sort(colSums(input_arc_CPM$data_loaded))) %>%
  filter(`sort(colSums(input_arc_CPM$data_loaded))` == 0)
input_arc_CPM_nz <- filter_data(input_arc_CPM,
                                filter_cat = "sampleID",
                                filter_vals = rownames(countArc))

# BC, PERMANOVA, PERMDISP, PCoA
bc_arc <- calc_dm(input_arc_CPM_nz$data_loaded)
set.seed(1150)
adonis2(bc_arc ~ Environment, data = input_arc_CPM_nz$map_loaded) # R2 = 0.19, p = 0.001
anova(betadisper(bc_arc, input_arc_CPM_nz$map_loaded$Environment)) # Dispersion not homogeneous
pcoa_arc <- cmdscale(bc_arc, k = nrow(input_arc_CPM_nz$map_loaded) - 1, eig = T)
pcoaA1arc <- round((eigenvals(pcoa_arc)/sum(eigenvals(pcoa_arc)))[1]*100, digits = 1)
pcoaA2arc <- round((eigenvals(pcoa_arc)/sum(eigenvals(pcoa_arc)))[2]*100, digits = 1)
input_arc_CPM_nz$map_loaded$Axis01 <- scores(pcoa_arc)[,1]
input_arc_CPM_nz$map_loaded$Axis02 <- scores(pcoa_arc)[,2]
micro.hulls <- ddply(input_arc_CPM_nz$map_loaded, c("Environment"), find_hull)
g_arc <- ggplot(input_arc_CPM_nz$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F, size = 0.25) +
  geom_point(size = 1, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1arc, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2arc, "%", sep = "")) +
  ggtitle("Archaea") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, vjust = -1))
g_arc

# Bacteria
input_bac <- filter_taxa_from_input(input,
                                    taxa_to_keep = "Bacteria",
                                    at_spec_level = 1)
nrow(input_bac$data_loaded) # 2682
input_bac_CPM <- input_bac
for (i in 1:ncol(input_bac$data_loaded)) {
  for (j in 1:nrow(input_bac$data_loaded)) {
    input_bac_CPM$data_loaded[j, i] <- (input_bac$data_loaded[j, i]*1000000)/input_bac$map_loaded$GenomeSize[i]
  }
}

# Check zeroes - none, all samples had bacteria
countbac <- as.data.frame(sort(colSums(input_bac_CPM$data_loaded))) %>%
  filter(`sort(colSums(input_bac_CPM$data_loaded))` == 0)

# BC, PERMANOVA, PERMDISP, PCoA
bc_bac <- calc_dm(input_bac_CPM$data_loaded)
set.seed(1150)
adonis2(bc_bac ~ Environment, data = input_bac_CPM$map_loaded) # R2 = 0.17, p = 0.001
anova(betadisper(bc_bac, input_bac_CPM$map_loaded$Environment)) # Dispersion not homogeneous
pcoa_bac <- cmdscale(bc_bac, k = nrow(input_bac_CPM$map_loaded) - 1, eig = T)
pcoaA1bac <- round((eigenvals(pcoa_bac)/sum(eigenvals(pcoa_bac)))[1]*100, digits = 1)
pcoaA2bac <- round((eigenvals(pcoa_bac)/sum(eigenvals(pcoa_bac)))[2]*100, digits = 1)
input_bac_CPM$map_loaded$Axis01 <- scores(pcoa_bac)[,1]
input_bac_CPM$map_loaded$Axis02 <- scores(pcoa_bac)[,2]
micro.hulls <- ddply(input_bac_CPM$map_loaded, c("Environment"), find_hull)
g_bac <- ggplot(input_bac_CPM$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F, size = 0.25) +
  geom_point(size = 1, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1bac, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2bac, "%", sep = "")) +
  ggtitle("Bacteria") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, vjust = -1))
g_bac

# Remake first one with fungi title
input_fungi_CPM_nz$map_loaded$Axis01 <- scores(pcoa)[,1]
input_fungi_CPM_nz$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(input_fungi_CPM_nz$map_loaded, c("Environment"), find_hull)
g <- ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F, size = 0.25) +
  geom_point(size = 1, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  ggtitle("Fungi") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, vjust = -1))
g

# Multipanel
pdf("FigsUpdated/PCoA_ArcBacFun.pdf", width = 8, height = 5)
plot_grid(g_arc,g_bac,g,g_leg, ncol = 2, hjust = "hv")
dev.off()



#### _Subset ####
# Lopsided sample sizes may be throwing things off
# Randomly sample 22 samples from each environment and redo
table(input_fungi_CPM_nz$map_loaded$Environment)
table(input_bac_CPM$map_loaded$Environment)
table(input_arc_CPM$map_loaded$Environment)

set.seed(1210)
subset22 <- input_fungi_CPM_nz$map_loaded %>%
  group_by(Environment) %>%
  slice_sample(n = 22)
table(subset22$Environment)
sum(table(subset22$Environment)) # 176

sub_arc <- filter_data(input_arc_CPM_nz,
                       filter_cat = "sampleID",
                       keep_vals = subset22$sampleID) # 173 (3 had zero Archaea)
sub_bac <- filter_data(input_bac_CPM,
                       filter_cat = "sampleID",
                       keep_vals = subset22$sampleID) # 176, good
sub_fun <- filter_data(input_fungi_CPM_nz,
                       filter_cat = "sampleID",
                       keep_vals = subset22$sampleID) # 176, good

bc_arc2 <- calc_dm(sub_arc$data_loaded)
set.seed(1150)
adonis2(bc_arc2 ~ Environment, data = sub_arc$map_loaded) # R2 = 0.28, p = 0.001
anova(betadisper(bc_arc2, sub_arc$map_loaded$Environment)) # Dispersion not homogeneous
pcoa_arc2 <- cmdscale(bc_arc2, k = nrow(sub_arc$map_loaded) - 1, eig = T)
pcoaA1arc2 <- round((eigenvals(pcoa_arc2)/sum(eigenvals(pcoa_arc2)))[1]*100, digits = 1)
pcoaA2arc2 <- round((eigenvals(pcoa_arc2)/sum(eigenvals(pcoa_arc2)))[2]*100, digits = 1)
sub_arc$map_loaded$Axis01 <- scores(pcoa_arc2)[,1]
sub_arc$map_loaded$Axis02 <- scores(pcoa_arc2)[,2]
micro.hulls <- ddply(sub_arc$map_loaded, c("Environment"), find_hull)
g_arc2 <- ggplot(sub_arc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F, size = 0.25) +
  geom_point(size = 1, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1arc2, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2arc2, "%", sep = "")) +
  ggtitle("Archaea") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, vjust = -1))
g_arc2

bc_bac2 <- calc_dm(sub_bac$data_loaded)
set.seed(1150)
adonis2(bc_bac2 ~ Environment, data = sub_bac$map_loaded) # R2 = 0.28, p = 0.001
anova(betadisper(bc_bac2, sub_bac$map_loaded$Environment)) # Dispersion not homogeneous
pcoa_bac2 <- cmdscale(bc_bac2, k = nrow(sub_bac$map_loaded) - 1, eig = T)
pcoaA1bac2 <- round((eigenvals(pcoa_bac2)/sum(eigenvals(pcoa_bac2)))[1]*100, digits = 1)
pcoaA2bac2 <- round((eigenvals(pcoa_bac2)/sum(eigenvals(pcoa_bac2)))[2]*100, digits = 1)
sub_bac$map_loaded$Axis01 <- scores(pcoa_bac2)[,1]
sub_bac$map_loaded$Axis02 <- scores(pcoa_bac2)[,2]
micro.hulls <- ddply(sub_bac$map_loaded, c("Environment"), find_hull)
g_bac2 <- ggplot(sub_bac$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F, size = 0.25) +
  geom_point(size = 1, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1bac2, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2bac2, "%", sep = "")) +
  ggtitle("Bacteria") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, vjust = -1))
g_bac2

bc_fun2 <- calc_dm(sub_fun$data_loaded)
set.seed(1150)
adonis2(bc_fun2 ~ Environment, data = sub_fun$map_loaded) # R2 = 0.13, p = 0.001
anova(betadisper(bc_fun2, sub_fun$map_loaded$Environment)) # Dispersion not homogeneous
pcoa_fun2 <- cmdscale(bc_fun2, k = nrow(sub_fun$map_loaded) - 1, eig = T)
pcoaA1fun2 <- round((eigenvals(pcoa_fun2)/sum(eigenvals(pcoa_fun2)))[1]*100, digits = 1)
pcoaA2fun2 <- round((eigenvals(pcoa_fun2)/sum(eigenvals(pcoa_fun2)))[2]*100, digits = 1)
sub_fun$map_loaded$Axis01 <- scores(pcoa_fun2)[,1]
sub_fun$map_loaded$Axis02 <- scores(pcoa_fun2)[,2]
micro.hulls <- ddply(sub_fun$map_loaded, c("Environment"), find_hull)
g_fun2 <- ggplot(sub_fun$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F, size = 0.25) +
  geom_point(size = 1, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1fun2, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2fun2, "%", sep = "")) +
  ggtitle("Fungi") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, vjust = -1))
g_fun2

pdf("FigsUpdated/PCoA_ArcBacFun_n22.pdf", width = 8, height = 5)
plot_grid(g_arc2,g_bac2,g_fun2,g_leg, ncol = 2, hjust = "hv")
dev.off()



#### _Methods ####
# Check to see if methods differences are obscuring beta-diversity patterns
# Key metadata columns are Add.Date, Assembly.Method
bc <- calc_dm(input_fungi_CPM_nz$data_loaded)
pcoa <- cmdscale(bc, k = nrow(input_fungi_CPM_nz$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
input_fungi_CPM_nz$map_loaded$Axis01 <- scores(pcoa)[,1]
input_fungi_CPM_nz$map_loaded$Axis02 <- scores(pcoa)[,2]

# Year
input_fungi_CPM_nz$map_loaded$Year <- as.factor(as.character(input_fungi_CPM_nz$map_loaded$Year))
table(input_fungi_CPM_nz$map_loaded$Year)
set.seed(1150)
adonis2(bc ~ Year, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.15, p = 0.001
anova(betadisper(bc, input_fungi_CPM_nz$map_loaded$Year)) # Dispersion not homogeneous
micro.hulls <- ddply(input_fungi_CPM_nz$map_loaded, c("Year"), find_hull)
g_year <- ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Year, fill = Year),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 1, alpha = 0.5, aes(colour = Year),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  ggtitle("Year: R2 = 0.15") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "right",
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))
g_year

# Assembly method
# Need to make new column for method type or multiple or unknown
table(input_fungi_CPM_nz$map_loaded$`Assembly Method`)
input_fungi_CPM_nz$map_loaded$Assembler <- dplyr::recode_factor(input_fungi_CPM_nz$map_loaded$`Assembly Method`,
                                   "AbySS v1.5.0" = "AbySS",
                                   "canu v. 1.7" = "Canu",
                                   "canu v. 1.9" = "Canu",
                                   "Celera WGS Assembler v5.3" = "Celera",
                                   "Celera, Phrap" = "Multiple",
                                   "CLC Bio package" = "CLC",
                                   "CLC Genomics" = "CLC",
                                   "CLC genomics wb7" = "CLC",
                                   "CLC genomics workbench" = "CLC",
                                   "CLC genomics workbench, v. 7.0" = "CLC",
                                   "Custom JGI assembly, Nielsen et. al." = "Custom JGI",
                                   "Custom JGI assembly." = "Custom JGI",
                                   "IDBA" = "IDBA",
                                   "IDBA 1.1.1 PRE correction" = "IDBA",
                                   "idba v. 1.1.1" = "IDBA",
                                   "IDBA v. 1.1.1" = "IDBA",
                                   "IDBA_1.1.1 PRE_correction" = "IDBA",
                                   "IDBA_UD" = "IDBA_UD",
                                   "idba_ud 1.1.1" = "IDBA_UD",
                                   "IDBA_UD k52-k92 step10" = "IDBA_UD",
                                   "IDBA_UD k52-k92 step10 - contigs greater than or equal to 1kb and less than 4kb" = "IDBA_UD",
                                   "IDBA_UD k52-k92 step10 - contigs greater than or equal to 4kb" = "IDBA_UD",
                                   "IDBA_UD k52-k92 step10 - contigs less than 1kb" = "IDBA_UD",
                                   "idba_ud v 1.1.1" = "IDBA_UD",
                                   "idba_ud v. 1.1.1" = "IDBA_UD",
                                   "IDBA_UD v. 1.1.1" = "IDBA_UD",
                                   "IDBA-UD" = "IDBA_UD",
                                   "IDBA-UD 1.1.3" = "IDBA_UD",
                                   "IDBA-ud v. 1.1.1" = "IDBA_UD",
                                   "IDBA-UD v. 1.1.1" = "IDBA_UD",
                                   "idba.1.1.1" = "IDBA",
                                   "lucy / pga" = "Multiple",
                                   "Megahit" = "MEGAHIT",
                                   "MEGAHit" = "MEGAHIT",
                                   "MegaHit v. 1.02" = "MEGAHIT",
                                   "MEGAHIT v. 1.1.1" = "MEGAHIT",
                                   "megahit v. 1.1.3" = "MEGAHIT",
                                   "Megahit v. 1.1.3" = "MEGAHIT",
                                   "MegaHIT v. 1.2.9" = "MEGAHIT",
                                   "MEGAHIT v. MEGAHIT v0.2.0" = "MEGAHIT",
                                   "MEGAHIT v. MEGAHIT v1.0.3" = "MEGAHIT",
                                   "MEGAHIT v. MEGAHIT v1.0.6" = "MEGAHIT",
                                   "MEGAHIT v.1.0.3" = "MEGAHIT",
                                   "MEGAN" = "MEGAN",
                                   "metaSPAdes" = "metaSPAdes",
                                   "metaSPAdes v. 3.10.0" = "metaSPAdes",
                                   "metaspades v. 3.13.0" = "metaSPAdes",
                                   "metaspades v. 3.14.1" = "metaSPAdes",
                                   "metaSPAdes v. 3.7.1" = "metaSPAdes",
                                   "metaSPAdes v3.10, CLC genomic workbench v7.5.1" = "metaSPAdes",
                                   "metaSPAdes v3.10.1" = "metaSPAdes",
                                   "Metavelvet" = "MetaVelvet",
                                   "MetaVelvet - v1.2.01" = "MetaVelvet",
                                   "MetaVelvet 1.2.02" = "MetaVelvet",
                                   "metavelvet v. 1.2.02" = "MetaVelvet",
                                   "Metavelvet v. 1.2.02" = "MetaVelvet",
                                   "mira 3.0.4" = "MIRA",
                                   "MIRA 4.9.5_2" = "MIRA",
                                   "Newbler" = "Newbler",
                                   "Newbler and/or Velvet" = "Multiple",
                                   "Newbler v. 2.5" = "Newbler",
                                   "Newbler v. 2.5.3" = "Newbler",
                                   "Newbler v2.7" = "Newbler",
                                   "pga" = "PGA",
                                   "Ray 2.3.1" = "Ray",
                                   "Ray 2.3.1 (no min length)" = "Ray",
                                   "reassembled with IDBA_UD" = "IDBA_UD",
                                   "reassembly by IDBA-UD" = "IDBA_UD",
                                   "reassembly with IDBA_UD" = "IDBA_UD",
                                   "reassembly with IDBA-UD" = "IDBA_UD",
                                   "SAPDES" = "SPAdes",
                                   "Soap denovo and minimus" = "Multiple",
                                   "SOAPdenovo v2.04" = "SOAPdenovo",
                                   "SOAPdenovo,newbler,minimus2 v. Version 1.05: testing... 2010,(v2.8 (20120726_1306)),AMOS/3.1.0" = "SOAPdenovo",
                                   "Spades" = "SPAdes",
                                   "SPADES" = "SPAdes",
                                   "spades 3.0" = "SPAdes",
                                   "Spades 3.6.1" = "SPAdes",
                                   "SPAdes 3.7.1" = "SPAdes",
                                   "SPAdes 3.8.0" = "SPAdes",
                                   "Spades v. 3.10" = "SPAdes",
                                   "SPAdes v. 3.11.0" = "SPAdes",
                                   "spades v. 3.11.1" = "SPAdes",
                                   "Spades v. 3.11.1" = "SPAdes",
                                   "SPAdes v. 3.11.1" = "SPAdes",
                                   "spades v. 3.12.0" = "SPAdes",
                                   "Spades v. 3.12.0" = "SPAdes",
                                   "spades v. 3.13.0" = "SPAdes",
                                   "SPades v. 3.13.0" = "SPAdes",
                                   "SPAdes v. 3.6.0" = "SPAdes",
                                   "SPAdes v. 3.9.0" = "SPAdes",
                                   "spades v. SPAdes version: 3.10.1" = "SPAdes",
                                   "spades v. SPAdes version: 3.11.1-check" = "SPAdes",
                                   "SPADES v3.6.1" = "SPAdes",
                                   "Spades_3.6" = "SPAdes",
                                   "SPAdes3.1.0" = "SPAdes",
                                   "Unknown" = "Unknown",
                                   "Unkown" = "Unknown",
                                   "Velvet" = "Velvet",
                                   "Velvet + MetaVelvet at multiple Kmers" = "Multiple",
                                   "Velvet + MetaVelvet at multiple Kmers followed by Minimus2 on all assemblies" = "Multiple",
                                   "Velvet, MetaVelvet, Minimus" = "Multiple")
levels(input_fungi_CPM_nz$map_loaded$Assembler)
table(input_fungi_CPM_nz$map_loaded$Assembler)

# Get these columns and sample ID to add to any metadata
methods <- input_fungi_CPM_nz$map_loaded %>%
  dplyr::select(sampleID, Year, Assembler)

set.seed(1150)
adonis2(bc ~ Assembler, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.16, p = 0.001
anova(betadisper(bc, input_fungi_CPM_nz$map_loaded$Assembler)) # Dispersion not homogeneous
micro.hulls <- ddply(input_fungi_CPM_nz$map_loaded, c("Assembler"), find_hull)
g_assem <- ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Assembler, fill = Assembler),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 1, alpha = 0.5, aes(colour = Assembler),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  ggtitle("Assembler: R2 = 0.16") +
  guides(colour = guide_legend(override.aes = list(alpha = 1),
                               ncol = 1)) +
  theme_bw() +  
  theme(legend.position = "right",
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))
g_assem

pdf("FigsUpdated/PCoA_Year_Assembler.pdf", width = 9, height = 4)
plot_grid(g_year, g_assem, ncol = 2, rel_widths = c(0.4, 0.6))
dev.off()

# Plot year and assembler numbers
year_df <- as.data.frame(table(input_fungi_CPM_nz$map_loaded$Year))
g_numyear <- ggplot(year_df, aes(Var1, Freq)) +
  geom_bar(stat = "identity") +
  labs(x = "Year", 
       y = "Number of samples") +
  ggtitle("Sample size by year") +
  scale_y_continuous(limits = c(0, 350),
                     expand = c(0.01, 0.01)) +
  theme_bw() +  
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        plot.title = element_text(vjust = -1))
g_numyear

assem_df <- as.data.frame(table(input_fungi_CPM_nz$map_loaded$Assembler))
g_numassem <- ggplot(assem_df, aes(reorder(Var1, Freq, mean), Freq)) +
  geom_bar(stat = "identity") +
  labs(x = "Assembler", 
       y = "Number of samples") +
  ggtitle("Sample size by assembler") +
  scale_y_continuous(limits = c(0, 350),
                                expand = c(0.01, 0.01)) +
  theme_bw() +  
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        plot.title = element_text(vjust = -1))
g_numassem

pdf("FigsUpdated/SampleSize_Year_Assembler.pdf", width = 9, height = 4)
plot_grid(g_numyear, g_numassem, align = "hv", ncol = 2)
dev.off()



#### _by Ecosystem ####
# Subset the data into each environment
# Make pie charts of taxa
# Use fungi only input data - get relative abundances of fungal phyla out of just fungi
# Use samples with at least 1 fungal count

# Need to do several for loops
# For subseting, summarizing, plotting, store dfs in a list to enable for loop/indexing
env <- list()
studyN <- list()
df <- list()
phy <- list()
p_colors <- list()
p <- list()

# Subset by Environment (8 data frames)
for (i in 1:length(levels(input_fungi_nz$map_loaded$Environment))) {
  env[[i]] <- filter_data(input_fungi_nz,
                          filter_cat = "Environment",
                          keep_vals = levels(input_fungi_nz$map_loaded$Environment)[i])
}

# Do one for each, or by individual studies? How many studies in each env?
length(levels(env[[1]]$map_loaded$`Study Name`)) # 5
length(levels(env[[2]]$map_loaded$`Study Name`)) # 13
length(levels(env[[3]]$map_loaded$`Study Name`)) # 7
length(levels(env[[4]]$map_loaded$`Study Name`)) # 1
length(levels(env[[5]]$map_loaded$`Study Name`)) # 26
length(levels(env[[6]]$map_loaded$`Study Name`)) # 44
length(levels(env[[7]]$map_loaded$`Study Name`)) # 31
length(levels(env[[8]]$map_loaded$`Study Name`)) # 5

# Probably best to show some of the variability within environment type
# For example, Shu and Huang 2022 have 4-5 sites for each environment
# Here let's do 1-6, for env. with more than 6 studies, take 6 with greatest sample size
# Could also decide to show sites with most fungi instead of most samples
for (i in 1:length(env)) {
  studyN[[i]] <- as.data.frame(table(env[[i]]$map_loaded$`Study Name`)) %>%
    arrange(desc(Freq)) %>%
    slice_head(n = 6)
}

studies <- rbind(studyN[[1]], studyN[[2]], studyN[[3]], studyN[[4]],
                 studyN[[5]], studyN[[6]], studyN[[7]], studyN[[8]]) %>%
  mutate(Var1 = as.character(Var1))

# Subset Environments by Studies
counter <- 1
for (i in 1:length(env)) {
  k <- nrow(studyN[[i]])
  for (l in 1:k) {
  df[[counter]] <- filter_data(env[[i]],
                               filter_cat = "Study Name",
                               keep_vals = studyN[[i]]$Var1[l])
  counter <- counter + 1
  }
}

# Location and n (for plot titles), do manually (hard to automate this part)
loc_n <- as.data.frame(matrix(NA, nrow = length(df), ncol = 3)) %>%
  set_names(c("Location", "n", "Environment"))
for (i in 1:length(df)) {
  loc_n$Location[i] <- levels(df[[i]]$map_loaded$`Geographic Location`)[1]
  loc_n$n[i] <- nrow(df[[i]]$map_loaded)
  loc_n$Environment[i] <- levels(df[[i]]$map_loaded$Environment)[1]
}

df[[1]]$map_loaded$Location <- "Richmond Mine\nUSA\n(n = 17)"
df[[2]]$map_loaded$Location <- "Malanjkhand copper mine\nIndia\n(n = 10)"
df[[3]]$map_loaded$Location <- "Los Rueldos mercury mine\nSpain\n(n = 3)"
df[[4]]$map_loaded$Location <- "Fankou Mines\nChina\n(n = 1)"
df[[5]]$map_loaded$Location <- "Akron\nUSA\n(n = 1)"
df[[6]]$map_loaded$Location <- "Glacial meltwater/mats\nAntarctica\n(n = 17)"
df[[7]]$map_loaded$Location <- "Cryoconites\nGreenland\n(n = 12)"
df[[8]]$map_loaded$Location <- "Glacier sediment\nAntarctica\n(n = 6)"
df[[9]]$map_loaded$Location <- "Glacial ice\nCanada\n(n = 4)"
df[[10]]$map_loaded$Location <- "Glacial sediment\nCanada/Iceland\n(n = 3)"
df[[11]]$map_loaded$Location <- "Cryoconite\nItaly\n(n = 2)"
df[[12]]$map_loaded$Location <- "Temperate Desert soil\nCalifornia\n(n = 56)"
df[[13]]$map_loaded$Location <- "Polar Desert soil\nAntarctica\n(n = 32)"
df[[14]]$map_loaded$Location <- "Temperate Desert soil\nUtah\n(n = 25)"
df[[15]]$map_loaded$Location <- "Temperate Desert soil\nNew Mexico\n(n = 18)"
df[[16]]$map_loaded$Location <- "Temperate Desert soil\nUtah\n(n = 8)"
df[[17]]$map_loaded$Location <- "Temperate Desert soil\nUtah\n(n = 6)"
df[[18]]$map_loaded$Location <- "Glacial forefield soil\nSweden/Norway/Greenland\n(n = 49)"
df[[19]]$map_loaded$Location <- "Yellowstone hot springs\nUSA\n(n = 34)"
df[[20]]$map_loaded$Location <- "Various hot springs\nUSA/Canada/China/South Africa\n(n = 33)"
df[[21]]$map_loaded$Location <- "Waikite Valley hot spring mats\nNew Zealand\n(n = 18)"
df[[22]]$map_loaded$Location <- "Yellowstone hot spring mats\nUSA\n(n = 33)"
df[[23]]$map_loaded$Location <- "Great Boiling Spring sediment\nUSA\n(n = 7)"
df[[24]]$map_loaded$Location <- "Yellowstone hot spring sediment\nUSA\n(n = 7)"
df[[25]]$map_loaded$Location <- "Guaymas Basin sediment/mats\nMexico\n(n = 47)"
df[[26]]$map_loaded$Location <- "Mid Cayman Rise\n(n = 27)"
df[[27]]$map_loaded$Location <- "Various vents\nPacific/Atlantic\n(n = 21)"
df[[28]]$map_loaded$Location <- "Mid Cayman Rise plume\n(n = 14)"
df[[29]]$map_loaded$Location <- "Various vents\nPacific\n(n = 13)"
df[[30]]$map_loaded$Location <- "Axial seamount\n(n = 12)"
df[[31]]$map_loaded$Location <- "Various lakes\nAustralia\n(n = 117)"
df[[32]]$map_loaded$Location <- "Organic Lake\nAntarctica\n(n = 34)"
df[[33]]$map_loaded$Location <- "Salton Sea\nUSA\n(n = 10)"
df[[34]]$map_loaded$Location <- "Salterns\nNamibia\n(n = 7)"
df[[35]]$map_loaded$Location <- "Bras del Port saltern\nSpain\n(n = 6)"
df[[36]]$map_loaded$Location <- "Bras del Port saltern\nSpain\n(n = 6)"
df[[37]]$map_loaded$Location <- "Alkaline sediment\nRussia/Germany\n(n = 11)"
df[[38]]$map_loaded$Location <- "Alkaline water\nItaly/Philippines/Costa Rica\n(n = 8)"
df[[39]]$map_loaded$Location <- "Bras del Port saltern\nSpain\n(n = 1)"
df[[40]]$map_loaded$Location <- "Tuz Lake\nTurkey\n(n = 1)"
df[[41]]$map_loaded$Location <- "Diamante Lake biofilm\nAgrgentina\n(n = 1)"

# Phyla
for (i in 1:nrow(studies)) {
  phy[[i]] <- summarize_taxonomy(df[[i]], level = 2, report_higher_tax = F, relative = T)
}

# Colors
phyla_present <- as.data.frame(table(input_fungi_nz$taxonomy_loaded$taxonomy2)) %>%
  mutate(Var1 = as.character(Var1))
for (i in 1:length(df)) {
  p_colors[[i]] <- brewer.pal(7, "Paired")
  p_colors[[i]] <- subset(p_colors[[i]],
                         phyla_present$Var1 %in% df[[i]]$taxonomy_loaded$taxonomy2)
}

# Pies
for (i in 1:nrow(studies)) {
  t <- df[[i]]$map_loaded$Location
  p[[i]] <- plot_taxa_bars(phy[[i]], df[[i]]$map_loaded, "Study Name", 20) +
    # geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = p_colors[[i]]) +
    theme_void() +
    ggtitle(t) +
    theme(legend.position = "none",
          plot.title = element_text(size = 5, hjust = 0.5),
          plot.margin = margin(0, 0, 0, 0, "pt"))
}

# Get legend - use one with all 7 phyla
t <- df[[7]]$map_loaded$Location
p_forleg <- plot_taxa_bars(phy[[7]], df[[7]]$map_loaded, "Study Name", 20) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start=0) +
  scale_fill_manual(values = p_colors[[7]]) +
  labs(fill = "Phylum") +
  theme_void() +
  ggtitle(t) +
  theme(legend.position = "right",
        # legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(size = 10, hjust = 0.5, vjust = -15),
        plot.margin = margin(0,0,0,0, "pt"))
p_forleg
pie_leg <- get_legend(p_forleg)

# Make huge multipanel
# 5, 6, 6, 1, 6, 6, 6, 5
pies <- plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], NULL,
                  p[[6]], p[[7]], p[[8]], p[[9]], p[[10]], p[[11]],
                  p[[12]], p[[13]], p[[14]], p[[15]], p[[16]], p[[17]],
                  p[[18]], NULL, NULL, NULL, NULL, NULL,
                  p[[19]], p[[20]], p[[21]], p[[22]], p[[23]], p[[24]],
                  p[[25]], p[[26]], p[[27]], p[[28]], p[[29]], p[[30]],
                  p[[31]], p[[32]], p[[33]], p[[34]], p[[35]], p[[36]],
                  p[[37]], p[[38]], p[[39]], p[[40]], p[[41]], NULL,
                  ncol = 6)
pies

# Add legend
pdf("FigsUpdated/Pies_Phyla.pdf", width = 8.5, height = 6.5)
plot_grid(pies, pie_leg, rel_widths = c(4, 1))
dev.off()

# Now do at class level
cla <- list()
c_colors <- list()
c <- list()

# Classes
for (i in 1:nrow(studies)) {
  cla[[i]] <- summarize_taxonomy(df[[i]], level = 3, report_higher_tax = F, relative = T)
}

# Colors
classes_present <- as.data.frame(table(input_fungi_nz$taxonomy_loaded$taxonomy3)) %>%
  mutate(Var1 = as.character(Var1))
for (i in 1:length(df)) {
  c_colors[[i]] <- colorRampPalette(brewer.pal(12, "Paired"))(37)
  c_colors[[i]] <- subset(c_colors[[i]],
                          classes_present$Var1 %in% df[[i]]$taxonomy_loaded$taxonomy3)
}

# Pies
for (i in 1:nrow(studies)) {
  t <- df[[i]]$map_loaded$Location
  c[[i]] <- plot_taxa_bars(cla[[i]], df[[i]]$map_loaded, "Study Name", 20) +
    # geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = c_colors[[i]]) +
    theme_void() +
    ggtitle(t) +
    theme(legend.position = "none",
          plot.title = element_text(size = 5, hjust = 0.5),
          plot.margin = margin(0, 0, 0, 0, "pt"))
}

# Get legend - use one with all 37 classes
t <- df[[7]]$map_loaded$Location
c_forleg <- plot_taxa_bars(cla[[7]], df[[7]]$map_loaded, "Study Name", 37) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start=0) +
  scale_fill_manual(values = c_colors[[7]]) +
  labs(fill = "Classes") +
  theme_void() +
  ggtitle(t) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(size = 10, hjust = 0.5, vjust = -15),
        plot.margin = margin(0,0,0,0, "pt"))
c_forleg
pie_leg_c <- get_legend(c_forleg)

# Make huge multipanel
pies_c <- plot_grid(c[[1]], c[[2]], c[[3]], c[[4]], c[[5]], NULL,
                    c[[6]], c[[7]], c[[8]], c[[9]], c[[10]], c[[11]],
                    c[[12]], c[[13]], c[[14]], c[[15]], c[[16]], c[[17]],
                    c[[18]], NULL, NULL, NULL, NULL, NULL,
                    c[[19]], c[[20]], c[[21]], c[[22]], c[[23]], c[[24]],
                    c[[25]], c[[26]], c[[27]], c[[28]], c[[29]], c[[30]],
                    c[[31]], c[[32]], c[[33]], c[[34]], c[[35]], c[[36]],
                    c[[37]], c[[38]], c[[39]], c[[40]], c[[41]], NULL,
                    ncol = 6)
pies_c

# Add legend
pdf("FigsUpdated/Pies_Classes.pdf", width = 8.5, height = 6.5)
plot_grid(pies_c, pie_leg_c, rel_widths = c(4, 1))
dev.off()



#### _Alpha ####
# Look at number of different taxa levels by environments
# OTU Richness
input_fungi$map_loaded$rich <- specnumber(input_fungi$data_loaded, 
                                          MARGIN = 2)

# Shannon diversity
input_fungi$map_loaded$shannon <- diversity(input_fungi$data_loaded, 
                                            index = "shannon", 
                                            MARGIN = 2)

# Stats and graphs
leveneTest(input_fungi$map_loaded$rich ~ input_fungi$map_loaded$Environment) # Bad
m3 <- aov(rich ~ Environment, data = input_fungi$map_loaded)
summary(m3)
TukeyHSD(m3)
shapiro.test(m3$residuals) # Bad
kruskal.test(rich ~ Environment, data = input_fungi$map_loaded)
nyi3 <- kwAllPairsNemenyiTest(rich ~ Environment, data = input_fungi$map_loaded)
nyi_table3 <- fullPTable(nyi3$p.value)
nyi_list3 <- multcompLetters(nyi_table3)
nyi_let3 <- as.data.frame(nyi_list3$Letters) %>%
  mutate(label = `nyi_list3$Letters`,
         y = rep(300, nrow(.)),
         name = "rich") %>%
  dplyr::select(-`nyi_list3$Letters`) %>%
  rownames_to_column(var = "Environment")

leveneTest(input_fungi$map_loaded$shannon ~ input_fungi$map_loaded$Environment) # Bad
m4 <- aov(shannon ~ Environment, data = input_fungi$map_loaded)
shapiro.test(m4$residuals) # Bad
summary(m4)
TukeyHSD(m4)
nyi4 <- kwAllPairsNemenyiTest(shannon ~ Environment, data = input_fungi$map_loaded)
nyi_table4 <- fullPTable(nyi4$p.value)
nyi_list4 <- multcompLetters(nyi_table4)
nyi_let4 <- as.data.frame(nyi_list4$Letters) %>%
  mutate(label = `nyi_list4$Letters`,
         y = rep(5.5, nrow(.)),
         name = "shannon") %>%
  dplyr::select(-`nyi_list4$Letters`) %>%
  rownames_to_column(var = "Environment")

label_df <- rbind(nyi_let3, nyi_let4)
facet_df <- c("rich" = "(a) Genus richness",
              "shannon" = "(b) Genus Shannon")
alpha_long <- input_fungi$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
pdf("FigsUpdated/AlphaDiversity.pdf", width = 6, height = 3)
ggplot(alpha_long, aes(reorder(Environment, value, mean), value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.5, alpha = 0.2, width = 0.3) +
  geom_text(data = label_df, aes(Environment, y, label = label), 
            size = 4, color = "black") +
  labs(x = NULL, y = NULL) +
  facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        strip.text = element_text(size = 10))
dev.off()



#### _Drivers/Indicators ####
# MULTIPATT (list taxa associated with each group)
# Phyla
tax_sum_phyla <- summarize_taxonomy(input_fungi_CPM, level = 2, report_higher_tax = F, relative = F)
set.seed(425)
mp_phyla <- multipatt(t(tax_sum_phyla), 
                      input_fungi_CPM$map_loaded$Environment, 
                      func = "r.g", 
                      control = how(nperm=999))
summary(mp_phyla) # 5 soda lake

pdf("FigsUpdated/mp_Phyla.pdf", width = 5, height = 5)
plot_multipatt(mp_obj = mp_phyla, 
               input = input_fungi_CPM,
               tax_sum = tax_sum_phyla,
               group = "Environment")
dev.off()

# Class
tax_sum_class <- summarize_taxonomy(input_fungi_CPM, level = 3, report_higher_tax = F, relative = F)
set.seed(425)
mp_class <- multipatt(t(tax_sum_class), 
                      input_fungi_CPM$map_loaded$Environment, 
                      func = "r.g", 
                      control = how(nperm=999))
summary(mp_class) # 23 soda lake
pdf("FigsUpdated/mp_Class.pdf", width = 5, height = 5)
plot_multipatt(mp_obj = mp_class, 
               input = input_fungi_CPM,
               tax_sum = tax_sum_class,
               group = "Environment")
dev.off()

# Order
tax_sum_order <- summarize_taxonomy(input_fungi_CPM, level = 4, report_higher_tax = F, relative = F)
set.seed(425)
mp_order <- multipatt(t(tax_sum_order), 
                      input_fungi_CPM$map_loaded$Environment, 
                      func = "r.g", 
                      control = how(nperm=999))
summary(mp_order) # Acid mine 4, Cyrosphere 12, Soda lake 41, 
pdf("FigsUpdated/mp_Order.pdf", width = 5, height = 7)
plot_multipatt(mp_obj = mp_order, 
               input = input_fungi_CPM,
               tax_sum = tax_sum_order,
               group = "Environment")
dev.off()

# Family
tax_sum_family <- summarize_taxonomy(input_fungi_CPM, level = 5, report_higher_tax = F, relative = F)
set.seed(425)
mp_family <- multipatt(t(tax_sum_family), 
                      input_fungi_CPM$map_loaded$Environment, 
                      func = "r.g", 
                      control = how(nperm=999))
summary(mp_family) # Acid mine 4, crysosphere 31, desert 2, glacial forefield 1, soda lake 69
pdf("FigsUpdated/mp_Family.pdf", width = 5, height = 10)
plot_multipatt(mp_obj = mp_family, 
               input = input_fungi_CPM,
               tax_sum = tax_sum_family,
               group = "Environment")
dev.off()

# Genus
tax_sum_genus <- summarize_taxonomy(input_fungi_CPM, level = 6, report_higher_tax = F, relative = F)
set.seed(425)
mp_genus <- multipatt(t(input_fungi_CPM$data_loaded), 
                      input_fungi_CPM$map_loaded$Environment, 
                      func = "r.g", 
                      control = how(nperm=999))
summary(mp_genus) # Acid mine 8, cyrosphere 44, desert 7, glacial forefield 1, hot spring 2, soda lake 112
pdf("FigsUpdated/mp_Genus.pdf", width = 5, height = 12)
plot_multipatt(mp_obj = mp_genus, 
               input = input_fungi_CPM,
               tax_sum = tax_sum_genus,
               group = "Environment")
dev.off()


#### ..........................####
#### Functional ####
# Sent Dongying Wu of IMG staff list of taxonoids and list of fungal phyla
# Dongying ran custom python scripts on JGI super computer to pull out KOs of only fungal phyla
# Folder FungalKOs has a file for each metagenome with the KO hits of the fungal phyla scaffolds
# Already deleted 318 blank files (no fungal KOs); 837 had at least 1 KO
# Note that there is bias in eukaryote gene calling/KO assignment
# We could also get COG or Pfam profiles if we want

# Run a for loop to read in the file for each metagenome and combine into 1
setwd("FungalKOs/")
list.files()
ko <- list()
ko_input <- data.frame(V1 = "NA",
                       V2 = "NA",
                       V3 = "NA")
ko_table <- ko_input
for (i in 1:length(list.files())) {
  ko[[i]] <- read.delim(list.files()[i], header = F)
  ko_table <- ko_table %>%
    rbind(ko[[i]])
}

setwd("~/Documents/GitHub/Extremophilic_Fungi/")

# Add new desert samples
new_kos <- read.table("fungi.ko.txt")
ko_table <- rbind(ko_table, new_kos)

# Clean up table
ko_table_wTax <- ko_table %>%
  filter(V1 != "NA") %>%
  separate(V1, into = c("Junk", "KO"), sep = ":") %>%
  dplyr::select(-Junk, -V2) %>%
  mutate(taxon_oid = substring(V3, first = 1, last = 10)) %>%
  separate(V3, into = c("Junk", "taxonomy"), sep = "Eukaryota;") %>%
  dplyr::select(-Junk) %>%
  separate(taxonomy, 
           into = c("Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = ";") %>%
  dplyr::select(taxon_oid, KO, everything())

# KO count (abundance) by metagenome
ko_table_MGcount <- ko_table_wTax %>%
  dplyr::select(taxon_oid, KO) %>%
  group_by(taxon_oid, KO) %>%
  summarise(KO_abund = n()) %>%
  pivot_wider(id_cols = KO, names_from = taxon_oid, values_from = KO_abund) %>%
  column_to_rownames(var = "KO")
ko_table_MGcount[is.na(ko_table_MGcount) == TRUE] <- 0

# Export KO list to upload to KEGG and get classifications
# For Definitions, can also use the keggFind function below
# write.csv(ko_table$KO, file = "KOlist.csv", row.names = F)

# List of KOs sorted by overall abundance
ko_list <- ko_table %>%
  filter(V1 != "NA") %>%
  dplyr::select(V1) %>%
  group_by(V1) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(Definition = "NA") %>%
  separate(V1, into = c("Junk", "KO"), sep = ":", remove = F) %>%
  dplyr::select(-Junk)

# Add definitions to list for only new KOs and append to old list
old_ko_list <- read.csv("KOlist_wDefinitions.csv")
new_ko_list <- ko_list %>%
  filter(KO %notin% old_ko_list$KO)
for (i in 1:nrow(new_ko_list)) {
  def <- keggFind(database = "ko", query = new_ko_list$KO[i])
  if (length(def) != 0) {
    new_ko_list$Definition[i] <- def
  }
}
#write.csv(new_ko_list, file = "KOlist_wDefinitions_new.csv", row.names = F)
new_ko_list <- read.csv("KOlist_wDefinitions_new.csv")
ko_list <- rbind(old_ko_list, new_ko_list)
ko_list$KO_def <- paste(ko_list$KO, ko_list$Definition, sep = " ")

# Make community style table and metadata, match IDs
ko_comm <- ko_table_MGcount %>%
  t() %>%
  as.data.frame() %>%
  filter(rownames(.) %in% input_fungi$map_loaded$taxon_oid) %>%
  arrange(rownames(.))

ko_meta <- input_fungi$map_loaded %>%
  filter(taxon_oid %in% rownames(ko_comm)) %>%
  arrange(taxon_oid) %>%
  left_join(., methods, by = "sampleID") %>%
  mutate(rn = sampleID) %>%
  column_to_rownames(var = "rn")

# Check match (should be zero)
sum(rownames(ko_comm) != ko_meta$taxon_oid)

# Check environment sample size
table(ko_meta$Environment)

# Get richness and Shannon
ko_meta$richness_KO = specnumber(ko_comm)
ko_meta$shannon_KO = diversity(ko_comm, index = "shannon")
range(ko_meta$richness_KO)
range(ko_meta$shannon_KO)

pdf("FigsUpdated/KO_Genus_richness.pdf", width = 6, height = 4)
ggplot(ko_meta, aes(rich, richness_KO)) +
  geom_point(size = 1.5, alpha = 0.25, aes(colour = Environment)) +
  geom_smooth() +
  labs(x = "Fungal genus richness", 
       y = "Fungal KO richness") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
dev.off()

# DESeq Normalization (run once, then reload from saved)
#dds <- DESeqDataSetFromMatrix(countData = t(ko_comm) + 1,
#                              colData = ko_meta,
#                              design = ~ 1)
#dds <- estimateSizeFactors(dds)
#dds <- estimateDispersions(dds)
#ko_comm_DESeq <- as.data.frame(t(counts(dds, normalized = T)))
# Save so you don't have to redo the DESeq (takes a while)
# saveRDS(ko_comm_DESeq, "ko_comm_DESeq_updated.rds")
ko_comm_DESeq <- readRDS("ko_comm_DESeq_updated.rds")



#### _KO Alpha ####
leveneTest(richness_KO ~ Environment, data = ko_meta) # Bad
m5 <- aov(richness_KO ~ Environment, data = ko_meta)
summary(m5)
shapiro.test(m5$residuals) # Bad
TukeyHSD(m5)
kruskal.test(richness_KO ~ Environment, data = ko_meta)
nyi5 <- kwAllPairsNemenyiTest(richness_KO ~ Environment, data = ko_meta)
nyi_table5 <- fullPTable(nyi5$p.value)
nyi_list5 <- multcompLetters(nyi_table5)
nyi_let5 <- as.data.frame(nyi_list5$Letters) %>%
  mutate(label = `nyi_list5$Letters`,
         y = rep(2510, nrow(.)),
         name = "richness_KO") %>%
  dplyr::select(-`nyi_list5$Letters`) %>%
  rownames_to_column(var = "Environment")

leveneTest(shannon_KO ~ Environment, data = ko_meta) # Bad
m6 <- aov(shannon_KO ~ Environment, data = ko_meta)
summary(m6)
shapiro.test(m6$residuals) # Bad
TukeyHSD(m6)
kruskal.test(shannon_KO ~ Environment, data = ko_meta)
nyi6 <- kwAllPairsNemenyiTest(shannon_KO ~ Environment, data = ko_meta)
nyi_table6 <- fullPTable(nyi6$p.value)
nyi_list6 <- multcompLetters(nyi_table6)
nyi_let6 <- as.data.frame(nyi_list6$Letters) %>%
  mutate(label = `nyi_list6$Letters`,
         y = rep(8, nrow(.)),
         name = "shannon_KO") %>%
  dplyr::select(-`nyi_list6$Letters`) %>%
  rownames_to_column(var = "Environment")

label_df_ko <- rbind(nyi_let5, nyi_let6)
facet_df_ko <- c("richness_KO" = "(a) KO richness",
                 "shannon_KO" = "(b) KO Shannon")
alpha_long_ko <- ko_meta %>%
  pivot_longer(cols = c("richness_KO", "shannon_KO"))
pdf("FigsUpdated/KO_AlphaDiversity.pdf", width = 6, height = 3)
ggplot(alpha_long_ko, aes(reorder(Environment, value, mean), value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.5, alpha = 0.2, width = 0.3) +
  geom_text(data = label_df_ko, aes(Environment, y, label = label), 
            size = 4, color = "black") +
  labs(x = NULL, y = NULL) +
  facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df_ko)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        strip.text = element_text(size = 10))
dev.off()



#### _KO Beta ####
#### __KO PCoAs ####
bc_ko <- vegdist(ko_comm_DESeq, method = "bray")
pcoa_ko <- cmdscale(bc_ko, k = nrow(ko_meta) - 1, eig = T)
pcoaA1 <- round((eigenvals(pcoa_ko)/sum(eigenvals(pcoa_ko)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa_ko)/sum(eigenvals(pcoa_ko)))[2]*100, digits = 1)
ko_meta$Axis01 <- scores(pcoa_ko)[,1]
ko_meta$Axis02 <- scores(pcoa_ko)[,2]
micro.hulls <- ddply(ko_meta, c("Environment"), find_hull)
g1_ko <- ggplot(ko_meta, aes(Axis01, Axis02, colour = Environment)) +
  geom_polygon(data = micro.hulls, aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Site",
       title = "KO Bray-Curtis") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        plot.title = element_text(vjust = 0))
g1_ko

jac_ko <- vegdist(ko_comm_DESeq, method = "jaccard")
pcoa1_ko <- cmdscale(jac_ko, k = nrow(ko_meta) - 1, eig = T)
pcoa1A1 <- round((eigenvals(pcoa1_ko)/sum(eigenvals(pcoa1_ko)))[1]*100, digits = 1)
pcoa1A2 <- round((eigenvals(pcoa1_ko)/sum(eigenvals(pcoa1_ko)))[2]*100, digits = 1)
ko_meta$Axis01j <- scores(pcoa1_ko)[,1]
ko_meta$Axis02j <- scores(pcoa1_ko)[,2]
micro.hullsj <- ddply(ko_meta, c("Environment"), find_hullj)
g2_ko <- ggplot(ko_meta, aes(Axis01j, Axis02j, colour = Environment)) +
  geom_polygon(data = micro.hullsj, aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = paste("PC1: ", pcoa1A1, "%", sep = ""), 
       y = paste("PC2: ", pcoa1A2, "%", sep = ""),
       colour = "Environment",
       title = "KO Jaccard") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        plot.title = element_text(vjust = 0))
g2_ko

plot_grid(g1_ko, g2_ko, ncol = 2, rel_widths = c(1,1.515))

# These are incredibly similar and also show a dense cluster with barely any dissimilarity
# I suspect this is driven by KO count; metagenomes with only 1 KO can't be too dissimilar
# Check number of ones
ko_richness <- ko_meta %>%
  dplyr::select(richness_KO) %>%
  arrange(desc(richness_KO)) %>%
  mutate(index = seq(1:nrow(.)))
plot(rownames(ko_list), ko_list$n)
plot(ko_richness$index, ko_richness$richness_KO)
sum(ko_richness$richness_KO == 1) # 66 with just 1 KO
sum(ko_richness$richness_KO == 2) # 47 with just 2 KOs
sum(ko_richness$richness_KO == 3) # 50 with just 3 KOs

# Need to find good cutoff with some KOs and still high sample size
sum(ko_richness$richness_KO > 10) # 543 with > 10

# Subset the metadata and community datasets
ko_meta_filt <- ko_meta %>%
  filter(richness_KO > 10)
ko_comm_DESeq_filt <- ko_comm_DESeq %>%
  filter(rownames(.) %in% rownames(ko_meta_filt))
sum(rownames(ko_comm_DESeq_filt) != rownames(ko_meta_filt))
table(ko_meta_filt$Environment)

# Redo the above with the filtered dataset
bc_ko <- vegdist(ko_comm_DESeq_filt, method = "bray")
pcoa_ko <- cmdscale(bc_ko, k = nrow(ko_meta_filt) - 1, eig = T)
pcoaA1 <- round((eigenvals(pcoa_ko)/sum(eigenvals(pcoa_ko)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa_ko)/sum(eigenvals(pcoa_ko)))[2]*100, digits = 1)
ko_meta_filt$Axis01 <- scores(pcoa_ko)[,1]
ko_meta_filt$Axis02 <- scores(pcoa_ko)[,2]
micro.hulls <- ddply(ko_meta_filt, c("Environment"), find_hull)
g3_ko <- ggplot(ko_meta_filt, aes(Axis01, Axis02, colour = Environment)) +
  geom_polygon(data = micro.hulls, aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Site",
       title = "KO Bray-Curtis") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        plot.title = element_text(vjust = 0))
g3_ko

jac_ko <- vegdist(ko_comm_DESeq_filt, method = "jaccard")
pcoa1_ko <- cmdscale(jac_ko, k = nrow(ko_meta_filt) - 1, eig = T)
pcoa1A1 <- round((eigenvals(pcoa1_ko)/sum(eigenvals(pcoa1_ko)))[1]*100, digits = 1)
pcoa1A2 <- round((eigenvals(pcoa1_ko)/sum(eigenvals(pcoa1_ko)))[2]*100, digits = 1)
ko_meta_filt$Axis01j <- scores(pcoa1_ko)[,1]
ko_meta_filt$Axis02j <- scores(pcoa1_ko)[,2]
micro.hullsj <- ddply(ko_meta_filt, c("Environment"), find_hullj)
g4_ko <- ggplot(ko_meta_filt, aes(Axis01j, Axis02j, colour = Environment)) +
  geom_polygon(data = micro.hullsj, aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = paste("PC1: ", pcoa1A1, "%", sep = ""), 
       y = paste("PC2: ", pcoa1A2, "%", sep = ""),
       colour = "Environment",
       title = "KO Jaccard") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        plot.title = element_text(vjust = 0))
g4_ko

plot_grid(g3_ko, g4_ko, align = "hv", ncol = 2, rel_widths = c(1,1.515))

# Still looks weird, increase the cutoff and rerun
sum(ko_richness$richness_KO > 100) # 209 with > 100
ko_meta_filt <- ko_meta %>%
  filter(richness_KO > 100)
ko_comm_DESeq_filt <- ko_comm_DESeq %>%
  filter(rownames(.) %in% rownames(ko_meta_filt))
sum(rownames(ko_comm_DESeq_filt) != rownames(ko_meta_filt))
table(ko_meta_filt$Environment)

bc_ko <- vegdist(ko_comm_DESeq_filt, method = "bray")
pcoa_ko <- cmdscale(bc_ko, k = nrow(ko_meta_filt) - 1, eig = T)
pcoaA1 <- round((eigenvals(pcoa_ko)/sum(eigenvals(pcoa_ko)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa_ko)/sum(eigenvals(pcoa_ko)))[2]*100, digits = 1)
ko_meta_filt$Axis01 <- scores(pcoa_ko)[,1]
ko_meta_filt$Axis02 <- scores(pcoa_ko)[,2]
micro.hulls <- ddply(ko_meta_filt, c("Environment"), find_hull)
g5_ko <- ggplot(ko_meta_filt, aes(Axis01, Axis02, colour = Environment)) +
  geom_polygon(data = micro.hulls, aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Site",
       title = "KO Bray-Curtis") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        plot.title = element_text(vjust = 0))
g5_ko

jac_ko <- vegdist(ko_comm_DESeq_filt, method = "jaccard")
pcoa1_ko <- cmdscale(jac_ko, k = nrow(ko_meta_filt) - 1, eig = T)
pcoa1A1 <- round((eigenvals(pcoa1_ko)/sum(eigenvals(pcoa1_ko)))[1]*100, digits = 1)
pcoa1A2 <- round((eigenvals(pcoa1_ko)/sum(eigenvals(pcoa1_ko)))[2]*100, digits = 1)
ko_meta_filt$Axis01j <- scores(pcoa1_ko)[,1]
ko_meta_filt$Axis02j <- scores(pcoa1_ko)[,2]
micro.hullsj <- ddply(ko_meta_filt, c("Environment"), find_hullj)
g6_ko <- ggplot(ko_meta_filt, aes(Axis01j, Axis02j, colour = Environment)) +
  geom_polygon(data = micro.hullsj, aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = paste("PC1: ", pcoa1A1, "%", sep = ""), 
       y = paste("PC2: ", pcoa1A2, "%", sep = ""),
       colour = "Environment",
       title = "KO Jaccard") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        plot.title = element_text(vjust = 0))
g6_ko
ko_pcoa_leg <- get_legend(ggplot(ko_meta_filt, aes(Axis01j, Axis02j, colour = Environment)) +
  geom_polygon(data = micro.hullsj, aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = paste("PC1: ", pcoa1A1, "%", sep = ""), 
       y = paste("PC2: ", pcoa1A2, "%", sep = ""),
       colour = "Environment",
       title = "KO Jaccard") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        plot.title = element_text(vjust = 0)))

pdf("FigsUpdated/KO_PCoA_min100KOs.pdf", width = 8.5, height = 3.5)
plot_grid(g5_ko, g6_ko, ko_pcoa_leg, ncol = 3, rel_widths = c(2.5,2.5,1))
dev.off()

# Also do with the same subset of the data as taxonomy (22 samples from each env)
# Note numbers are slightly different here because of 0 fungal KOs in some samples
ko_meta_subset22 <- ko_meta %>%
  filter(sampleID %in% subset22$sampleID)
ko_comm_DESeq_subset22 <- ko_comm_DESeq %>%
  filter(rownames(.) %in% rownames(ko_meta_subset22))
sum(rownames(ko_comm_DESeq_subset22) != rownames(ko_meta_subset22))
table(ko_meta_subset22$Environment)

bc_ko <- vegdist(ko_comm_DESeq_subset22, method = "bray")
pcoa_ko <- cmdscale(bc_ko, k = nrow(ko_meta_subset22) - 1, eig = T)
pcoaA1 <- round((eigenvals(pcoa_ko)/sum(eigenvals(pcoa_ko)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa_ko)/sum(eigenvals(pcoa_ko)))[2]*100, digits = 1)
ko_meta_subset22$Axis01 <- scores(pcoa_ko)[,1]
ko_meta_subset22$Axis02 <- scores(pcoa_ko)[,2]
micro.hulls <- ddply(ko_meta_subset22, c("Environment"), find_hull)
g7_ko <- ggplot(ko_meta_subset22, aes(Axis01, Axis02, colour = Environment)) +
  geom_polygon(data = micro.hulls, aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Site",
       title = "KO Bray-Curtis") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        plot.title = element_text(vjust = 0))
g7_ko

jac_ko <- vegdist(ko_comm_DESeq_subset22, method = "jaccard")
pcoa1_ko <- cmdscale(jac_ko, k = nrow(ko_meta_subset22) - 1, eig = T)
pcoa1A1 <- round((eigenvals(pcoa1_ko)/sum(eigenvals(pcoa1_ko)))[1]*100, digits = 1)
pcoa1A2 <- round((eigenvals(pcoa1_ko)/sum(eigenvals(pcoa1_ko)))[2]*100, digits = 1)
ko_meta_subset22$Axis01j <- scores(pcoa1_ko)[,1]
ko_meta_subset22$Axis02j <- scores(pcoa1_ko)[,2]
micro.hullsj <- ddply(ko_meta_subset22, c("Environment"), find_hullj)
g8_ko <- ggplot(ko_meta_subset22, aes(Axis01j, Axis02j, colour = Environment)) +
  geom_polygon(data = micro.hullsj, aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = paste("PC1: ", pcoa1A1, "%", sep = ""), 
       y = paste("PC2: ", pcoa1A2, "%", sep = ""),
       colour = "Environment",
       title = "KO Jaccard") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        plot.title = element_text(vjust = 0))
g8_ko

pdf("FigsUpdated/KO_PCoA_n22.pdf", width = 8.5, height = 3.5)
plot_grid(g7_ko, g8_ko, ko_pcoa_leg, ncol = 3, rel_widths = c(2.5,2.5,1))
dev.off()



#### __Stats ####
# Rerun whichever richness filter you want, then run this
bc_ko <- vegdist(ko_comm_DESeq_filt, method = "bray")
jac_ko <- vegdist(ko_comm_DESeq_filt, method = "jaccard")
set.seed(308)
adonis2(bc_ko ~ Environment, data = ko_meta_filt) # R2 = 0.17, p = 0.001
anova(betadisper(bc_ko, ko_meta_filt$Environment)) # Dispersion not homogeneous

set.seed(308)
adonis2(jac_ko ~ ko_meta_filt$Environment) # R2 = 0.16, p = 0.001
anova(betadisper(jac_ko, ko_meta_filt$Environment)) # Dispersion not homogeneous

# Subset 22 (rerun that section first)
bc_ko <- vegdist(ko_comm_DESeq_subset22, method = "bray")
jac_ko <- vegdist(ko_comm_DESeq_subset22, method = "jaccard")
set.seed(308)
adonis2(bc_ko ~ Environment, data = ko_meta_subset22) # R2 = 0.10, p = 0.002
anova(betadisper(bc_ko, ko_meta_subset22$Environment)) # Dispersion not homogeneous

set.seed(308)
adonis2(jac_ko ~ ko_meta_subset22$Environment) # R2 = 0.10, p = 0.002
anova(betadisper(jac_ko, ko_meta_subset22$Environment)) # Dispersion not homogeneous

set.seed(308)
adonis2(bc_ko ~ ko_meta_subset22$Assembler) # R2 = 0.16, p = 0.084
anova(betadisper(bc_ko, ko_meta_subset22$Assembler)) # Dispersion not homogeneous

ko_meta_subset22 <- ko_meta_subset22 %>%
  separate(`Add Date`, into = c("Day", "Month", "Year"), sep = "/", remove = F)
set.seed(308)
adonis2(bc_ko ~ as.factor(ko_meta_subset22$Year)) # R2 = 0.05, p = 0.06
anova(betadisper(bc_ko, ko_meta_subset22$Year)) # Dispersion homogeneous



#### __Drivers/Indicators ####
# KO MULTIPATT (list KOs associated with each group)
set.seed(425)
mp <- multipatt(ko_comm_DESeq, 
                ko_meta$Environment, 
                func = "IndVal.g", 
                control = how(nperm=999))
summary(mp) # None!!



#### _KO Individual ####
# Extract and analyze list of KOs of interest
# Need to get list of KOs from Lara/Quandt Lab
# For now make barplot and heatmap of top 10 KOs to develop script
head(ko_list$KO, n = 10)
# To rerun with another list of KOs, just update the lines below
# Data frame
gene_plot <- data.frame("Environment" = as.factor(ko_meta$Environment),
                        K03164 = ko_comm_DESeq$K03164,
                        K10955 = ko_comm_DESeq$K10955,
                        K00777 = ko_comm_DESeq$K00777,
                        K06867 = ko_comm_DESeq$K06867,
                        K03006 = ko_comm_DESeq$K03006,
                        K03010 = ko_comm_DESeq$K03010,
                        K05658 = ko_comm_DESeq$K05658,
                        K00698 = ko_comm_DESeq$K00698,
                        K15503 = ko_comm_DESeq$K15503,
                        K12811 = ko_comm_DESeq$K12811)



#### __ (I) Stats ####
# Run a loop 
kruskal_results_genes <- as.data.frame(matrix(data = NA, 11, 3)) %>%
  set_names(c("Gene", "X2", "P"))
for (i in 2:11) {
  k <- kruskal.test(gene_plot[[i]] ~ gene_plot$Environment)
  kruskal_results_genes$Gene[i] <- names(gene_plot)[i]
  kruskal_results_genes$X2[i] <- round(k$statistic, digits = 2)
  kruskal_results_genes$P[i] <- k$p.value
}
kruskal_results_genes <- kruskal_results_genes %>%
  filter(is.na(Gene) == F) %>%
  mutate(Pfdr = p.adjust(P, method = "fdr"))
# All significant



#### __ (II) Graphs ####
# Barplot
table(ko_meta$Environment)
gene_plot_long <- gene_plot %>%
  pivot_longer(c("K03164", "K10955", "K00777", "K06867", "K03006",
                 "K03010", "K05658", "K00698", "K15503", "K12811"), 
               names_to = "Gene", values_to = "Abundance") %>%
  mutate(Gene = factor(Gene,
                       levels = c("K03164", "K10955", "K00777", "K06867", "K03006",
                                  "K03010", "K05658", "K00698", "K15503", "K12811"))) %>%
  droplevels() %>%
  mutate(Environment = dplyr::recode_factor(Environment,
                                     "Acid mine drainage" = "Acid mine drainage (n = 33)",
                                     "Cryosphere" = "Cryosphere (n = 47)",
                                     "Desert" = "Desert (n = 33)",
                                     "Glacial forefield" = "Glacial forefield (n = 60)",
                                     "Hot spring" = "Hot spring (n = 167)",
                                     "Hydrothermal vent" = "Hydrothermal vent (n = 237)",
                                     "Hypersaline" = "Hypersaline (n = 225)",
                                     "Soda lake" = "Soda lake (n = 25)"))

gene_plot_summary <- gene_plot_long %>%
  group_by(Environment, Gene) %>%
  summarise(mean = mean(Abundance),
            se = std.error(Abundance))

pdf("FigsUpdated/KO_Barplot_Top.pdf", width = 8, height = 5)
ggplot(gene_plot_summary, aes(Environment, mean, fill = Gene, group = Gene)) +
  geom_bar(stat = "identity", position = position_dodge(0.75)) +
  geom_linerange(aes(x = Environment, ymin = mean - se, ymax = mean + se, 
                     group = Gene),
                 position = position_dodge(0.75)) +
  labs(x = NULL, 
       y = "Abundance (DESeq2 normalized)",
       fill = "KO") +
  scale_fill_manual(values = brewer.pal(10, "Paired"),
                    labels = ko_list$KO_def[1:10]) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.7), "cm"))
dev.off()

# Heatmap
# Use pheatmap (pretty heatmap) package
gene_hm <- data.frame("sampleID" = ko_meta$sampleID,
                       K03164 = ko_comm_DESeq$K03164,
                       K10955 = ko_comm_DESeq$K10955,
                       K00777 = ko_comm_DESeq$K00777,
                       K06867 = ko_comm_DESeq$K06867,
                       K03006 = ko_comm_DESeq$K03006,
                       K03010 = ko_comm_DESeq$K03010,
                       K05658 = ko_comm_DESeq$K05658,
                       K00698 = ko_comm_DESeq$K00698,
                       K15503 = ko_comm_DESeq$K15503,
                       K12811 = ko_comm_DESeq$K12811) %>%
  column_to_rownames(var = "sampleID") %>%
  t() %>%
  as.matrix()

ann_cols <- data.frame(row.names = colnames(gene_hm), 
                       Environment = ko_meta$Environment)
ann_colors <- list(Environment = c("Acid mine drainage" = hue_pal()(8)[1],
                                   "Cryosphere" = hue_pal()(8)[2],
                                   "Desert" = hue_pal()(8)[3],
                                   "Glacial forefield" = hue_pal()(8)[4],
                                   "Hot spring" = hue_pal()(8)[5],
                                   "Hydrothermal vent" = hue_pal()(8)[6],
                                   "Hypersaline" = hue_pal()(8)[7],
                                   "Soda lake" = hue_pal()(8)[8]))
phm1 <- pheatmap(gene_hm,
         scale = "row",
         show_colnames = F,
         cluster_rows = F,
         annotation_col = ann_cols,
         annotation_colors = ann_colors)
save_pheatmap_pdf(phm1, "FigsUpdated/KO_heatmap_Top.pdf")



#### ..........................####
#### Other ####
#### _Map ####
# Sample map with ggplot, color by environment
# Need to adjust color scheme here and throughout
world <- map_data("world")
nrow(input$map_loaded)
input$map_loaded$Latitude <- as.numeric(input$map_loaded$Latitude)
input$map_loaded$Longitude <- as.numeric(input$map_loaded$Longitude)

# 21 missing
pdf("FigsUpdated/SampleMap.pdf", width = 8, height = 5)
ggplot() +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "white", fill = "lightgray", size = 0.1) +
  geom_point(data = input$map_loaded, aes(x = Longitude, y = Latitude,
                                          colour = Environment),
             size = 1, alpha = 0.5) +
  theme_void() +
  labs(x = NULL,
       y = NULL) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()



#### End Script ####