# Extremophilic Fungi Metagenome Analysis
# For review paper with Quandt Mycology Lab, University of Colorado Boulder
# by Cliff Bueno de Mesquita, JGI, Summer 2022, Updated Spring 2023, Finalized Summer 2023
# 
# This script is based on Extreme_Fungi.R, but updated for the new set of samples
# Then updated again with 9 additional cryosphere samples and updated environment type
# Click the "Show document outline" button in the top right corner to view document outline
# Sections are:
# 1. Retrieving Data
## - Notes on processing
# 2. Setup
# 3. Taxonomic analyses
# 4. Functional analyses
# 5. Other



#### 1. Retrieving Data ####
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
# Here also adding 9 more cryosphere samples by Lara



#### 2. Setup ####
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
suppressWarnings(suppressMessages(library(writexl))) # For writing to Excel
suppressWarnings(suppressMessages(library(ggpubr))) # For density plots
suppressWarnings(suppressMessages(library(PMCMRplus))) # For Nemenyi posthoc
suppressWarnings(suppressMessages(library(DirichletReg))) # For analyzing proportions
suppressWarnings(suppressMessages(library(MASS))) # For zinf reg
suppressWarnings(suppressMessages(library(pscl))) # For zinf reg
suppressWarnings(suppressMessages(library(boot))) # For zinf reg
suppressWarnings(suppressMessages(library(gamlss))) # For zinf reg
suppressWarnings(suppressMessages(library(writexl))) # Write Excel spreadsheets

# Working directory (GitHub repository)
setwd("~/Documents/GitHub/Extremophilic_Fungi/")

# Note directory is organized into folders for code and data

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

source("code/plot_multipatt.R")
source("code/compareBC.R")

# Custom color palette from Benjamin Young
color_mapping <- c(
  "Acid mine drainage" = "darkorange",
  "Cryosphere - soil" = "lightskyblue3",
  "Cryosphere - water" = "royalblue",
  "Desert" = "burlywood2",
  "Glacial forefield" = "grey70",
  "Hot spring" = "red",
  "Hydrothermal vent" = "firebrick4",
  "Hypersaline" = "plum2",
  "Soda lake" = "magenta3"
)

# So far used Kruskal-Wallis + Nemenyi posthoc
# May want to update some stats to zero-inflated negative binomial regression instead
# If using ANOVA and Tukey, can use this code to get sig. letters
# Example from genus richness
#tuk <- emmeans(object = m, specs = "Environment") %>%
#  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
#  mutate(name = "rich",
#     y = max(input_fungi$map_loaded$rich)+(max(input_fungi$map_loaded$rich)-min(input_fungi$map_loaded$rich))/20)


#### ..........................####
#### 3. Taxonomic ####

# Metadata ("mapping file") downloaded from IMG

# Tax table for mctoolsr (only need to do once, then skip to Import)
# This time a bit more complicated because need to full join the 1149 and the cryo9
# Should full join on the full taxonomy string and then make up an ASV variable
# Or reextract the genus level - but can do that later
c9 <- read.delim("data/Extreme_Cryo9/UI_data_output.txt") %>%
  dplyr::select(-c(3:8)) %>% # N.B.! Check columns! Don't delete column 9!
  dplyr::rename(taxonomy = FeatureName) %>%
  dplyr::rename(ASV_ID = Feature) %>%
  dplyr::select(ASV_ID, 3:ncol(.), taxonomy) %>%
  dplyr::select(-ASV_ID)
names(c9) <- abbreviate(names(c9), minlength = 11)

t <- read.delim("data/Extreme_Updated_1149/UI_data_output.txt") %>%
  dplyr::select(-c(3:9)) %>%
  dplyr::rename(taxonomy = FeatureName) %>%
  dplyr::rename(ASV_ID = Feature) %>%
  dplyr::select(ASV_ID, 3:ncol(.), taxonomy) %>%
  dplyr::select(-ASV_ID) %>%
  dplyr::full_join(., c9, by = "taxonomy") %>%
  dplyr::mutate(ASV_ID = paste("ASV_", rownames(.), sep = "")) %>%
  dplyr::select(ASV_ID, everything()) %>%
  dplyr::select(-taxonomy, taxonomy) %>%
  replace(is.na(.), 0)
names(t) <- abbreviate(names(t), minlength = 11)
table.fp <- "~/Documents/GitHub/Extremophilic_Fungi/data"
out_fp <- paste0(table.fp, "/genus_table_mctoolsr_final.txt")
names(t)[1] = "#ASV_ID"
write("#Exported for mctoolsr", out_fp)
suppressWarnings(write.table(t, out_fp, sep = "\t", row.names = FALSE, append = TRUE))



#### _Setup ####
# Import with mctoolsr (n = 932)

# metadata_final.txt = metadata_updated.txt + metadata_cryo9.txt - samplesToRemove (LV)
metadata_updated <- read.delim("data/metadata_updated.txt")

metadata_cryo9 <- read.delim("data/metadata_cryo9.txt") %>%
  mutate("IMG.Release.Pipeline.Version" = NA,
         "IMG.Submission.ID" = NA,
         "GOLD.Analysis.Project.ID" = NA,
         "GOLD.Sequencing.Project.ID" = NA,
         "GOLD.Study.ID" = NA,
         "NCBI.Assembly.Accession" = NA,
         "NCBI.Bioproject.Accession" = NA,
         "NCBI.Biosample.Accession" = NA,
         "NCBI.GenBank.ID" = NA,
         "pH" = NA,
         "Salinity" = NA) %>%
  dplyr::select(names(metadata_updated))

names(metadata_cryo9)
names(metadata_updated)
sum(names(metadata_updated) %notin% names(metadata_cryo9))

info <- read_excel("data/Extremophilic_fungi_dataset_final.xlsx") %>%
  filter(Use == "Use")

info_LV <- read_excel("data/Extremophilic_fungi_dataset_final_LV.xlsx", sheet = 3) %>%
  filter(Lara_remove == "remove")

# metadata_final_use.txt = delete samples according to Lara's review
metadata_final_use <- rbind(metadata_updated, metadata_cryo9) %>%
  filter(sampleID %in% info$sampleID) %>%
  filter(taxon_oid %notin% info_LV$taxon_oid)

nrow(metadata_final_use) # N.B.! n = 932 after Lara's review

# write
# write.table(metadata_final_use, 
#             file = "data/metadata_final_use.txt", 
#             sep = "\t",
#             row.names = F)



#### _Start Here ####
tax_table_fp <- file.path("data/genus_table_mctoolsr_final.txt")
map_fp <- file.path("data/metadata_final_use.txt")
input = load_taxa_table(tax_table_fp, map_fp) # 932
new_tab <- read_excel("data/Extremophilic_fungi_dataset_final.xlsx", sheet = 1) %>%
  mutate(Location2 = Geographic.Location) %>%
  dplyr::select(taxon_oid, Study.Name2, Location2, Environment)
table(new_tab$Environment)

# Update map_loaded - sampleID, GenomeSize, Environment, Assembly Method, Year
# Filter out samples from before 2012. 904 remaining
dim(input$map_loaded)
input$map_loaded <- input$map_loaded %>%
  mutate(sampleID = paste("X", taxon_oid, sep = ""),
         GenomeSize = `Genome.Size.....assembled`) %>%
  left_join(., new_tab, by = "taxon_oid") %>%
  separate(`Add.Date`, into = c("Day", "Month", "Year"), sep = "/", remove = F) %>%
  mutate(Year = as.integer(Year)) %>%
  filter(Year > 11) %>%
  mutate(sampleID2 = sampleID) %>%
  column_to_rownames(var = "sampleID2")
dim(input$map_loaded) # Now down to 904 samples
for (i in 1:nrow(input$map_loaded)) {
  if (input$map_loaded$`Assembly.Method`[1] == "") {
    input$map_loaded$`Assembly.Method`[i] <- "Unknown"
  }
  if (input$map_loaded$`Sequencing.Method`[i] == "") {
    input$map_loaded$`Sequencing.Method`[i] <- "Unknown"
  }
}
input$map_loaded <- input$map_loaded %>%
  mutate_if(is.character, as.factor)

# Filter out samples with no genus level reads (7 removed, 897 remaining)
count <- as.data.frame(sort(colSums(input$data_loaded))) %>%
  filter(`sort(colSums(input$data_loaded))` == 0)
input <- filter_data(input,
                     filter_cat = "sampleID",
                     filter_vals = rownames(count))

# Filter out 454 datasets (removes 42 samples. 855 remaining)
table(input$map_loaded$`Sequencing.Method`)
input <- filter_data(input,
                     filter_cat = "Sequencing.Method",
                     filter_vals = c("454", "454 GS FLX", "454 GS FLX Titanium",
                                     "454 GS FLX Titanium, Illumina GAIIx",
                                     "454 GS FLX Titanium, Illumina HiSeq 2000",
                                     "Illumina HiSeq 2000, 454 GS FLX Titanium"))

# Filter by Assembly method? For now no. Also many unknowns...
table(input$map_loaded$`Assembly.Method`)

# Filter out Plasmid:Bacteria, Plasmid:Eukaryota, Viruses, unknown at Domain level
input <- filter_taxa_from_input(input,
                                taxa_to_remove = c("Plasmid:Bacteria",
                                                   "Plasmid:Eukaryota",
                                                   "Viruses",
                                                   "unknown"),
                                at_spec_level = 1) # 1243 taxa removed

#write_xlsx(input$map_loaded, "data/TableS1.xlsx", format_headers = F)

# Check sequencing depth (# fungal assigned genes)
sort(colSums(input$data_loaded))
mean(colSums(input$data_loaded)) # 361621.3
se(colSums(input$data_loaded)) # 19512.22
input$map_loaded$count <- colSums(input$data_loaded)

# Genus reads
ggplot(input$map_loaded, aes(reorder(`Environment`, count, mean), count)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.25, width = 0.25) +
  labs(x = "Environment", 
       y = "# Reads") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

# Genome size
ggplot(input$map_loaded, aes(reorder(`Environment`, GenomeSize, median), log10(GenomeSize))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.25, width = 0.25) +
  labs(x = "Environment", 
       y = "Assembled metagenome size") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

# Log 10
ggplot(input$map_loaded, aes(reorder(`Environment`, GenomeSize, median), log10(GenomeSize))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.25, width = 0.25) +
  labs(x = "Environment", 
       y = "Assembled metagenome size") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

# Log axis - use this, old Figure 1d, now Figure S1
# To sort by "Fungi" need to run some code below first (line 426)
#### ___Figure S1 ####
gs <- ggplot(input$map_loaded, aes(reorder(`Environment`, Fungi, median), GenomeSize,
                                   colour = Environment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.5, width = 0.25, shape = 16) +
  scale_colour_manual(values = color_mapping) +
  labs(x = NULL, 
       y = "Assembly size (bp)") +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides = "l", outside = T, short = unit(1,"mm"), mid = unit(1,"mm"), 
                      long = unit(2,"mm")) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
gs
png("FinalFigs/FigureS1.png", width = 7, height = 5, units = "in", res = 300)
gs
dev.off()

max(input$map_loaded$GenomeSize)
min(input$map_loaded$GenomeSize)

# Size vs genus reads (actually genus genes)
pdf("FigsUpdated2/AssembledMetagenomeSizes_AllGenusReads.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(GenomeSize, count)) +
  geom_point(size = 1.5, alpha = 0.25) +
  geom_smooth(method = "lm") +
  labs(x = "Assembled genome size", 
       y = "Assigned genus genes") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10))
dev.off()
summary(lm(input$map_loaded$count ~ input$map_loaded$GenomeSize)) # R2 = 0.86



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
pdf("FigsUpdated2/EukTopSamples.pdf", width = 7, height = 5)
ggplot(topeuk$map_loaded, aes(reorder(sampleID, Euks, mean), Euks, fill = Environment)) +
  geom_bar(stat = "identity", color = NA) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x = NULL, 
       y = "Relative abundance") +
  ggtitle("Samples with Eukaryota > 5% (n = 33)") +
  theme_classic() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
dev.off()



#### __Check Fungi ####
# First do same as done above for euk but extract fungal phyla
tax_sum_phyla <- summarize_taxonomy(input, level = 2, report_higher_tax = T, relative = T)
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
pdf("FigsUpdated2/FungalTopSamples.pdf", width = 7, height = 5)
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
hist(input$map_loaded$Fungi)
hist(log(input$map_loaded$Fungi))
leveneTest(input$map_loaded$Fungi ~ input$map_loaded$Environment) # Bad
m1 <- aov(Fungi ~ Environment, data = input$map_loaded)
shapiro.test(m1$residuals) # Bad
hist(m1$residuals)
plot(m1$fitted.values, m1$residuals)
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
pdf("FigsUpdated2/FungalRel.pdf", width = 5, height = 4)
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

# Explore
ggplot(input$map_loaded, aes(reorder(Environment, Fungi, mean), log(Fungi))) +
  geom_violin(outlier.shape = NA) +
  geom_boxplot(outlier.shape = NA, width = 0.1) +
  geom_jitter(size = 1, alpha = 0.2, width = 0.4) +
  geom_text(data = nyi_let1, aes(Environment, y, label = label), 
            size = 4, color = "black") +
  labs(x = NULL, y = "Fungal relative abundance") +
  ylim(-14, 0.1) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        strip.text = element_text(size = 10))

# Try ANOVA with log?
# Tough to deal with zeroes

# Zero-inflated negative binomial?
# Can't use because numeric not integer
z1 <- zeroinfl(Fungi ~ Environment, data = input$map_loaded, dist = "negbin")
summary(z1)

# Zero-inflated Beta regression
d <- input$map_loaded # emmeans can't handle $
b1 <- gamlss(input$map_loaded$Fungi ~ d$Environment,  family = BEZI, trace = F)
summary(b1)
plot(b1)
emmeans(b1, "Environment", type = "response")
pairs(.Last.value)

tuk <- emmeans(object = b1, specs = "Environment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = rep(11, nrow(.))) %>%
  mutate(label = c("a", "a", "abc", "ab", "bc", "abc", "c", "abc", "bc"))
br_plot <- ggplot(input$map_loaded, aes(reorder(Environment, Fungi, median), Fungi*100,
                                        colour = Environment)) +
  geom_jitter(size = 2, alpha = 0.5, width = 0.4, shape = 16) +
  geom_text(data = tuk, aes(Environment, y, label = label), 
            size = 4, color = "black", hjust = 0.5) +
  labs(x = NULL, y = "Fungal % abund.") +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  scale_colour_manual(values = color_mapping) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 10))
br_plot

min(input$map_loaded$Fungi)
max(input$map_loaded$Fungi)
mean(input$map_loaded$Fungi)
se(input$map_loaded$Fungi)



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

nrow(input$taxonomy_loaded) # 5943 total
nrow(input_euk$taxonomy_loaded) # 439 euks
nrow(input_fungi$taxonomy_loaded) # 326 fungi

# Now check reads again (not actually reads, but fungal assigned CDS counts)
sort(colSums(input_fungi$data_loaded))
# Note lots of samples with 0 or very few fungi
# Purposefully not filtering those out those as 0's are interesting in this analysis
# These are extreme environments, some may have few to no fungi
# Further below, however, some analyses will be done with zeroes removed
mean(colSums(input_fungi$data_loaded)) # 556.2304
se(colSums(input_fungi$data_loaded)) # 63.54668
input_fungi$map_loaded$fung_count <- colSums(input_fungi$data_loaded)
input_fungi$map_loaded$present <- ifelse(input_fungi$map_loaded$fung_count > 0,
                                         1,
                                         0)
# Save file
#saveRDS(input_fungi, "input_fungi_updated2.rds")

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

# Check genome size vs fungal genus gene counts
pdf("FigsUpdated2/AssembledMetagenomeSizes_FungalGenusGeneCount.pdf", width = 7, height = 5)
ggplot(input_fungi$map_loaded, aes(GenomeSize, fung_count)) +
  geom_point(size = 1.5, alpha = 0.25) +
  geom_smooth(method = "lm") +
  labs(x = "Assembled genome size", 
       y = "Assigned fungal gene count") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10))
dev.off()
summary(lm(input_fungi$map_loaded$fung_count ~ input_fungi$map_loaded$GenomeSize)) 
# R2 = 0.13, p < 0.001

# Get and plot fraction of samples with 0 versus > 0 fungal reads, and also n
env_prev <- input_fungi$map_loaded %>%
  group_by(Environment) %>%
  summarise(num_present = sum(present),
            num_samples = n(),
            prevalence = round(num_present/num_samples * 100, digits = 2)) %>%
  mutate(num_absent = num_samples - num_present)
sum(env_prev$num_samples)
sum(env_prev$num_present)

# Melt for stacked bar
env_prev_long <- melt(env_prev,
                      id.vars = "Environment",
                      measure.vars = c("num_present", "num_absent"))

pdf("FigsUpdated2/FungalPrevalence.pdf", width = 5, height = 4)
ggplot(env_prev, aes(reorder(Environment, prevalence, mean), prevalence)) +
  geom_bar(stat = "identity") +
  labs(x = NULL, 
       y = "% prevalence of fungi") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

pdf("FigsUpdated2/SampleSize.pdf", width = 5, height = 4.5)
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
  ggtitle("Total sample size = 855\nSamples with Fungi = 732") +
  theme_bw() +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_blank(),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



#### ___Figure 1 ####
# Make multipanel Figure 1
# Update - make map a panel and move genome size to Figure S1
n <- ggplot(env_prev_long, aes(reorder(Environment, value, mean), value, 
                          group = Environment, fill = variable)) +
  geom_bar(stat = "identity") +
  geom_text(data = env_prev, 
            aes(reorder(Environment, num_samples, mean), num_samples+15, 
                label = num_samples), inherit.aes = F) +
  geom_text(aes(x = "Cryosphere - soil", y = 235, label = "total n = 855\nn with fungi = 732"),
            hjust = 0.5, check_overlap = T, size = 3) +
  scale_fill_manual(values = c("#F8766D", "#619CFF"),
                    breaks = c("num_absent", "num_present"),
                    labels = c("Absent", "Present")) +
  labs(x = NULL, 
       y = "Sample size",
       fill = "Fungi") +
  theme_classic() +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_blank(),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank())

prev <- ggplot(env_prev, aes(reorder(Environment, num_samples, mean), prevalence,
                             fill = Environment)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_mapping) +
  labs(x = NULL, 
       y = "Fungal % prev.") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

f1bd <- plot_grid(n, prev, ncol = 1, rel_heights = c(0.4, 0.6), labels = c("B", "D"))
f1bd

f1ac <- plot_grid(fig1a, br_plot, ncol = 1, rel_heights = c(0.4, 0.6), align = "v",
                 labels = c("A", "C"))
f1ac

#f1 <- plot_grid(f1a, br_plot, ncol = 2, labels = "AUTO")
#f1
f1 <- plot_grid(f1ac, f1bd, ncol = 2)
f1

png("FinalFigs/Figure1.png", width = 7, height = 5, units = "in", res = 300)
f1
dev.off()


#### _CPM ####
input_fungi_CPM <- input_fungi
# Replace counts in "data_loaded" with CPM transformed counts
# This is CPM assembled metagenomic base pairs
# Do (count*1000000)/GenomeSize
# i is samples 1 to 1142
# j is taxa 1 to 304
for (i in 1:ncol(input_fungi$data_loaded)) {
  for (j in 1:nrow(input_fungi$data_loaded)) {
    input_fungi_CPM$data_loaded[j, i] <- (input_fungi$data_loaded[j, i]*1000000)/input_fungi$map_loaded$GenomeSize[i]
  }
}

# Make stacked bar plots by taxonomic level
# Re-sort so unassigned and other are on the top
# Use "Paired" palette from RColorBrewer
# Note: "Set2" palette could be another option
# Unassigned gets grey75, Other gets grey90
# Show all phyla; for others show top 15 - use colorRampPalette to expand colors
# Can easily update this later if people want a different number of taxa shown
ntax <- 15
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(ntax)

# Phyla
tax_sum_Phyla <- summarize_taxonomy(input_fungi_CPM, 
                                    level = 2, 
                                    report_higher_tax = F, 
                                    relative = F)
barsPhyla <- plot_taxa_bars(tax_sum_Phyla,
                            input_fungi_CPM$map_loaded,
                            type_header = "Environment",
                            num_taxa = ntax,
                            data_only = T) %>%
  mutate(taxon = fct_rev(taxon))

pdf("FigsUpdated2/CPM_Phyla.pdf", width = 7, height = 5)
fs3a <- ggplot(barsPhyla, aes(reorder(group_by, mean_value, mean), mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Abundance (CPM)", fill = "Phylum") +
  scale_fill_manual(values = brewer.pal(12, "Paired")[7:1]) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.5, unit = "cm"),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_blank())
fs3a
dev.off()
ts_p <- taxa_summary_by_sample_type(tax_sum_Phyla, 
                                    input_fungi_CPM$map_loaded, 
                                    'Environment', 
                                    0.0001, 
                                    'KW')
#write_xlsx(ts_p, "data/TableS2.xlsx", format_headers = F)

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

tax_sum_Class <- summarize_taxonomy(input_fungi_CPM, level = 3, report_higher_tax = T, relative = F)
rownames(tax_sum_Class) <- substring(rownames(tax_sum_Class), 12)
barsClass <- plot_taxa_bars(tax_sum_Class,
               input_fungi_CPM$map_loaded,
               type_header = "Environment",
               num_taxa = ntax,
               data_only = T) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))

pdf("FigsUpdated2/CPM_Class.pdf", width = 7, height = 5)
fs3b <- ggplot(barsClass, aes(reorder(group_by, mean_value, mean), mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Abundance (CPM)", fill = "Class") +
  scale_fill_manual(values = c("grey90", mycolors[15:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.5, unit = "cm"),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_blank())
fs3b
dev.off()
ts_c <- taxa_summary_by_sample_type(tax_sum_Class, 
                                    input_fungi_CPM$map_loaded, 
                                    'Environment', 
                                    0.0001, 
                                    'KW')
# write_xlsx(ts_c, "data/TableS3.xlsx", format_headers = F)

# Multipanel Figure S3
# Horizontal
fs3b <- fs3b +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())
fs3 <- plot_grid(fs3a, fs3b, ncol = 2, labels = "AUTO", rel_widths = c(0.52, 0.48))
fs3

png("FinalFigs/FigureS3.png", width = 8, height = 6, units = "in", res = 300)
fs3
dev.off()

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

pdf("FigsUpdated2/CPM_Order.pdf", width = 7, height = 5)
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

pdf("FigsUpdated2/CPM_Family.pdf", width = 7, height = 5)
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

pdf("FigsUpdated2/CPM_Genus.pdf", width = 7, height = 5)
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
hist(input_fungi_CPM$map_loaded$totalFun)
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
pdf("FigsUpdated2/FungalCPM.pdf", width = 5, height = 4)
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
                                           input_fungi_CPM$map_loaded$`Geographic.Location`,
                                           sep = "_")
tax_sum_Phyla <- summarize_taxonomy(input_fungi_CPM, level = 2, report_higher_tax = F, relative = F)
barsPhyla <- plot_taxa_bars(tax_sum_Phyla,
                            input_fungi_CPM$map_loaded,
                            type_header = c("EnvGeo"),
                            num_taxa = ntax,
                            data_only = T) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Environment", "Location"), sep = "_")

pdf("FigsUpdated2/CPM_Phyla_EnvGeo.pdf", width = 12, height = 8)
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
# Remove zeroes (123 removed, 732 remaining)
countFun <- as.data.frame(sort(colSums(input_fungi_CPM$data_loaded))) %>%
  filter(`sort(colSums(input_fungi_CPM$data_loaded))` == 0)
input_fungi_nz <- filter_data(input_fungi,
                              filter_cat = "sampleID",
                              filter_vals = rownames(countFun))

# Calculate relative abundance of fungi, only for samples with fungi (n = 732)
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

pdf("FigsUpdated2/Rel_Phyla.pdf", width = 7, height = 5)
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

pdf("FigsUpdated2/Rel_Class.pdf", width = 7, height = 5)
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

pdf("FigsUpdated2/Rel_Order.pdf", width = 7, height = 5)
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

pdf("FigsUpdated2/Rel_Family.pdf", width = 7, height = 5)
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

pdf("FigsUpdated2/Rel_Genus.pdf", width = 7, height = 5)
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

# Filter out fungal zeroes from CPM data (filters 123 samples, 732 remaining)
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
adonis2(bc ~ Environment, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.15, p = 0.001
set.seed(1150)
adonis2(bc ~ Habitat, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.40, p = 0.001
set.seed(1150)
adonis2(bc ~ `Ecosystem.Category`, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.02, p = 0.001
set.seed(1150)
adonis2(bc ~ `Ecosystem.Subtype`, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.18, p = 0.001
set.seed(1150)
adonis2(bc ~ `Ecosystem.Type`, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.09, p = 0.001
set.seed(1150)
adonis2(bc ~ `Specific.Ecosystem`, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.12, p = 0.001
anova(betadisper(bc, input_fungi_CPM_nz$map_loaded$Environment)) # Dispersion not homogeneous
pcoa <- cmdscale(bc, k = nrow(input_fungi_CPM_nz$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
input_fungi_CPM_nz$map_loaded$Axis01 <- scores(pcoa)[,1]
input_fungi_CPM_nz$map_loaded$Axis02 <- scores(pcoa)[,2]
input_fungi_CPM_nz$map_loaded$Axis03 <- scores(pcoa)[,3]
micro.hulls <- ddply(input_fungi_CPM_nz$map_loaded, c("Environment"), find_hull)
pdf("FigsUpdated2/PCoA_Genus.pdf", width = 7, height = 5)
g <- ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, Axis02)) +
  #geom_polygon(data = micro.hulls, 
  #             aes(colour = Environment, fill = Environment),
  #             alpha = 0.1, show.legend = F) +
  geom_point(size = 2, alpha = 0.5, shape = 16, aes(colour = Environment),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_colour_manual(values = color_mapping) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        panel.grid = element_blank())
g
dev.off()

# Try ellipse? Bad. Ellipses are huge. Too messy.
ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, Axis02)) +
  stat_ellipse(mapping = aes(x = Axis01, y = Axis02, color = Environment), alpha = 0.8) +
  geom_point(size = 2, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))

# Interactive plot (can hover mouse over points)
ggplotly(g)

# There's basically no clustering of fungal composition by environment at genus level, lots of overlap!
# Extract legend to plot separately as its own panel later
g_leg <- get_legend(g)

# Family
bc_Family <- calc_dm(Family_nz)
set.seed(1150)
adonis2(bc_Family ~ Environment, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.16 , p = 0.001
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
adonis2(bc_Order ~ Environment, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.17 p = 0.001 
anova(betadisper(bc_Order, input_fungi_CPM_nz$map_loaded$Environment)) # Dispersion not homogeneous
pcoa_Order <- cmdscale(bc_Order, k = nrow(input_fungi_CPM_nz$map_loaded) - 1, eig = T)
pcoaA1O <- round((eigenvals(pcoa_Order)/sum(eigenvals(pcoa_Order)))[1]*100, digits = 1)
pcoaA2O <- round((eigenvals(pcoa_Order)/sum(eigenvals(pcoa_Order)))[2]*100, digits = 1)
input_fungi_CPM_nz$map_loaded$Axis01 <- scores(pcoa_Order)[,1]
input_fungi_CPM_nz$map_loaded$Axis02 <- scores(pcoa_Order)[,2]
micro.hulls <- ddply(input_fungi_CPM_nz$map_loaded, c("Environment"), find_hull)
g2 <- ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, -Axis02)) +
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
adonis2(bc_Class ~ Environment, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.18, p = 0.001 
anova(betadisper(bc_Class, input_fungi_CPM_nz$map_loaded$Environment)) # Dispersion not homogeneous
pcoa_Class <- cmdscale(bc_Class, k = nrow(input_fungi_CPM_nz$map_loaded) - 1, eig = T)
pcoaA1C <- round((eigenvals(pcoa_Class)/sum(eigenvals(pcoa_Class)))[1]*100, digits = 1)
pcoaA2C <- round((eigenvals(pcoa_Class)/sum(eigenvals(pcoa_Class)))[2]*100, digits = 1)
input_fungi_CPM_nz$map_loaded$Axis01 <- scores(pcoa_Class)[,1]
input_fungi_CPM_nz$map_loaded$Axis02 <- scores(pcoa_Class)[,2]
micro.hulls <- ddply(input_fungi_CPM_nz$map_loaded, c("Environment"), find_hull)
g3 <- ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, -Axis02)) +
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
adonis2(bc_Phylum ~ Environment, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.22, p = 0.001 
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
pdf("FigsUpdated2/PCoA_AllLevels.pdf", width = 8, height = 5)
plot_grid(g,g1,g2,g3,g4,g_leg, ncol = 3, hjust = "hv")
dev.off()



#### __Prok ####
# Now that we've seen not much clustering of the fungal communities, I'm curious to see how archaeal and bacterial communities look, just for the sake of comparison.

# Archaea
input_arc <- filter_taxa_from_input(input,
                                    taxa_to_keep = "Archaea",
                                    at_spec_level = 1)
nrow(input_arc$data_loaded) # 187
input_arc_CPM <- input_arc
for (i in 1:ncol(input_arc$data_loaded)) {
  for (j in 1:nrow(input_arc$data_loaded)) {
    input_arc_CPM$data_loaded[j, i] <- (input_arc$data_loaded[j, i]*1000000
                                        )/input_arc$map_loaded$GenomeSize[i]
  }
}

# Remove zeroes (42 removed, 813 remaining)
countArc <- as.data.frame(sort(colSums(input_arc_CPM$data_loaded))) %>%
  filter(`sort(colSums(input_arc_CPM$data_loaded))` == 0)
input_arc_CPM_nz <- filter_data(input_arc_CPM,
                                filter_cat = "sampleID",
                                filter_vals = rownames(countArc))

# BC, PERMANOVA, PERMDISP, PCoA
bc_arc <- calc_dm(input_arc_CPM_nz$data_loaded)
set.seed(1150)
adonis2(bc_arc ~ Environment, data = input_arc_CPM_nz$map_loaded) # R2 = 0.22, p = 0.001
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
  geom_point(size = 2, alpha = 0.75, shape = 16, aes(colour = Environment),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1arc, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2arc, "%", sep = "")) +
  scale_fill_manual(values = color_mapping) +
  scale_colour_manual(values = color_mapping) +
  ggtitle("Archaea") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, vjust = -1),
        panel.grid = element_blank())
g_arc

# Bacteria
input_bac <- filter_taxa_from_input(input,
                                    taxa_to_keep = "Bacteria",
                                    at_spec_level = 1)
nrow(input_bac$data_loaded) # 5317 (minus 2 (unclassified, unknown) = 5315 known)
input_bac_CPM <- input_bac
for (i in 1:ncol(input_bac$data_loaded)) {
  for (j in 1:nrow(input_bac$data_loaded)) {
    input_bac_CPM$data_loaded[j, i] <- (input_bac$data_loaded[j, i]*1000000
                                        )/input_bac$map_loaded$GenomeSize[i]
  }
}

# Check zeroes - none, all samples had bacteria
countbac <- as.data.frame(sort(colSums(input_bac_CPM$data_loaded))) %>%
  filter(`sort(colSums(input_bac_CPM$data_loaded))` == 0)

# BC, PERMANOVA, PERMDISP, PCoA
bc_bac <- calc_dm(input_bac_CPM$data_loaded)
set.seed(1150)
adonis2(bc_bac ~ Environment, data = input_bac_CPM$map_loaded) # R2 = 0.21, p = 0.001
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
  geom_point(size = 2, alpha = 0.75, shape = 16, aes(colour = Environment),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1bac, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2bac, "%", sep = "")) +
  scale_fill_manual(values = color_mapping) +
  scale_colour_manual(values = color_mapping) +
  ggtitle("Bacteria") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, vjust = -1),
        panel.grid = element_blank())
g_bac

# Remake first one with fungi title
input_fungi_CPM_nz$map_loaded$Axis01 <- scores(pcoa)[,1]
input_fungi_CPM_nz$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(input_fungi_CPM_nz$map_loaded, c("Environment"), find_hull)
g <- ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F, size = 0.25) +
  geom_point(size = 2, alpha = 0.75, shape = 16, aes(colour = Environment),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_fill_manual(values = color_mapping) +
  scale_colour_manual(values = color_mapping) +
  ggtitle("Fungi") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, vjust = -1),
        panel.grid = element_blank())
g

# Multipanel
pdf("FigsUpdated2/PCoA_ArcBacFun.pdf", width = 8, height = 5)
plot_grid(g_arc,g_bac,g,g_leg, ncol = 2, hjust = "hv", 
          labels = c("A", "B", "C", ""),
          label_x = 0.1)
dev.off()

png("FinalFigs/FigureS2.png", width = 8, height = 6, units = "in", res = 300)
plot_grid(g_arc,g_bac,g,g_leg, ncol = 2, hjust = "hv", 
          labels = c("A", "B", "C", ""),
          label_x = 0.1)
dev.off()



#### _Subset ####
# Lopsided sample sizes may be throwing things off
# Randomly sample 19 samples from each environment and redo
table(input_fungi_CPM_nz$map_loaded$Environment)
table(input_bac_CPM$map_loaded$Environment)
table(input_arc_CPM$map_loaded$Environment)

set.seed(1210)
subset19 <- input_fungi_CPM_nz$map_loaded %>%
  group_by(Environment) %>%
  slice_sample(n = 19)
table(subset19$Environment)
sum(table(subset19$Environment)) # 171

sub_arc <- filter_data(input_arc_CPM_nz,
                       filter_cat = "sampleID",
                       keep_vals = subset19$sampleID) # 165 (6 had zero Archaea)
sub_bac <- filter_data(input_bac_CPM,
                       filter_cat = "sampleID",
                       keep_vals = subset19$sampleID) # 171, good
sub_fun <- filter_data(input_fungi_CPM_nz,
                       filter_cat = "sampleID",
                       keep_vals = subset19$sampleID) # 171, good

bc_arc2 <- calc_dm(sub_arc$data_loaded)
set.seed(1150)
adonis2(bc_arc2 ~ Environment, data = sub_arc$map_loaded) # R2 = 0.35, p = 0.001
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
adonis2(bc_bac2 ~ Environment, data = sub_bac$map_loaded) # R2 = 0.34, p = 0.001
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
adonis2(bc_fun2 ~ Environment, data = sub_fun$map_loaded) # R2 = 0.17, p = 0.001
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

pdf("FigsUpdated2/PCoA_ArcBacFun_n19.pdf", width = 8, height = 5)
plot_grid(g_arc2,g_bac2,g_fun2,g_leg, ncol = 2, hjust = "hv")
dev.off()



#### _Compare BC ####
# Let's dig into the archaeal, bacterial, fungal communities a bit more
# How does Bray-Curtis dissimilarity compare for within or between sample type comparisons?
# Note, the length of the long dataframe should equal (n*(n-1))/2
# Make archaea and bacteria datasets same as fungal for fair comparison
nrow(input_fungi_CPM_nz$map_loaded) # 732
nrow(input_arc_CPM_nz$map_loaded) # 813
nrow(input_bac_CPM$map_loaded) # 855 (all)
sum(input_fungi_CPM_nz$map_loaded$sampleID %in% input_arc_CPM_nz$map_loaded$sampleID)
# 718 samples have all 3 domains. Do for these

bc_comp_fun <- filter_data(input_fungi_CPM_nz,
                           "sampleID",
                           keep_vals = input_arc_CPM_nz$map_loaded$sampleID)
bc_comp_bac <- filter_data(input_bac_CPM,
                           "sampleID",
                           keep_vals = bc_comp_fun$map_loaded$sampleID)
bc_comp_arc <- filter_data(input_arc_CPM_nz,
                           "sampleID",
                           keep_vals = bc_comp_fun$map_loaded$sampleID)

# BC
bc_fun <- calc_dm(bc_comp_fun$data_loaded)
bc_bac <- calc_dm(bc_comp_bac$data_loaded)
bc_arc <- calc_dm(bc_comp_arc$data_loaded)



#### __Genus ####
# Fungi
fun_bray_mat <- as.matrix(bc_fun)
fun_bray_mat[upper.tri(fun_bray_mat, diag = TRUE)] <- NA
fun_bray_df <- as.data.frame(fun_bray_mat)
fun_bray_df$sampleID <- rownames(fun_bray_df)
fun_bray_df_long <- melt(fun_bray_df, id.vars = "sampleID")
fun_bray_df_long <- na.omit(fun_bray_df_long)
fun_bray_df_long$sampleID <- as.factor(fun_bray_df_long$sampleID)
env <- dplyr::select(bc_comp_fun$map_loaded, sampleID, Environment)
fun_bray_df_long <- inner_join(fun_bray_df_long, env, 
                               by = c("sampleID" = "sampleID"))
fun_bray_df_long <- inner_join(fun_bray_df_long, env, 
                               by = c("variable" = "sampleID"))
for (i in 1:nrow(fun_bray_df_long)) {
  ifelse(fun_bray_df_long$Environment.x[i] == fun_bray_df_long$Environment.y[i],
         fun_bray_df_long$comparison[i] <- "within",
         fun_bray_df_long$comparison[i] <- "between")
}
fun_bray_df_long$comparison <- as.factor(fun_bray_df_long$comparison)
table(fun_bray_df_long$comparison)
shapiro.test(fun_bray_df_long$value[1:5000])
t.test(value ~ comparison, data = fun_bray_df_long)
wilcox.test(value ~ comparison, data = fun_bray_df_long)

# Bacteria
bac_bray_mat <- as.matrix(bc_bac)
bac_bray_mat[upper.tri(bac_bray_mat, diag = TRUE)] <- NA
bac_bray_df <- as.data.frame(bac_bray_mat)
bac_bray_df$sampleID <- rownames(bac_bray_df)
bac_bray_df_long <- melt(bac_bray_df, id.vars = "sampleID")
bac_bray_df_long <- na.omit(bac_bray_df_long)
bac_bray_df_long$sampleID <- as.factor(bac_bray_df_long$sampleID)
env <- dplyr::select(bc_comp_bac$map_loaded, sampleID, Environment)
bac_bray_df_long <- inner_join(bac_bray_df_long, env, 
                               by = c("sampleID" = "sampleID"))
bac_bray_df_long <- inner_join(bac_bray_df_long, env, 
                               by = c("variable" = "sampleID"))
for (i in 1:nrow(bac_bray_df_long)) {
  ifelse(bac_bray_df_long$Environment.x[i] == bac_bray_df_long$Environment.y[i],
         bac_bray_df_long$comparison[i] <- "within",
         bac_bray_df_long$comparison[i] <- "between")
}
bac_bray_df_long$comparison <- as.factor(bac_bray_df_long$comparison)
table(bac_bray_df_long$comparison)
shapiro.test(bac_bray_df_long$value[1:5000])
t.test(value ~ comparison, data = bac_bray_df_long)
wilcox.test(value ~ comparison, data = bac_bray_df_long)

# Archaea
arc_bray_mat <- as.matrix(bc_arc)
arc_bray_mat[upper.tri(arc_bray_mat, diag = TRUE)] <- NA
arc_bray_df <- as.data.frame(arc_bray_mat)
arc_bray_df$sampleID <- rownames(arc_bray_df)
arc_bray_df_long <- melt(arc_bray_df, id.vars = "sampleID")
arc_bray_df_long <- na.omit(arc_bray_df_long)
arc_bray_df_long$sampleID <- as.factor(arc_bray_df_long$sampleID)
env <- dplyr::select(bc_comp_arc$map_loaded, sampleID, Environment)
arc_bray_df_long <- inner_join(arc_bray_df_long, env, 
                               by = c("sampleID" = "sampleID"))
arc_bray_df_long <- inner_join(arc_bray_df_long, env, 
                               by = c("variable" = "sampleID"))
for (i in 1:nrow(arc_bray_df_long)) {
  ifelse(arc_bray_df_long$Environment.x[i] == arc_bray_df_long$Environment.y[i],
         arc_bray_df_long$comparison[i] <- "within",
         arc_bray_df_long$comparison[i] <- "between")
}
arc_bray_df_long$comparison <- as.factor(arc_bray_df_long$comparison)
table(arc_bray_df_long$comparison)
shapiro.test(arc_bray_df_long$value[1:5000])
t.test(value ~ comparison, data = arc_bray_df_long)
wilcox.test(value ~ comparison, data = arc_bray_df_long)

# Combine and Graph
fun_bray_df_long$taxon <- "Fungi"
bac_bray_df_long$taxon <- "Bacteria"
arc_bray_df_long$taxon <- "Archaea"
combined <- rbind(fun_bray_df_long, bac_bray_df_long, arc_bray_df_long)
combined$taxcomp <- paste(combined$taxon, combined$comparison, sep = "_")
hist(combined$value)
leveneTest(value ~ taxcomp, data = combined)
m <- aov(value ~ taxcomp, data = combined)
summary(m)
TukeyHSD(m)
shapiro.test(m$residuals[1:5000])
plot(m$fitted.values, m$residuals)
m <- aov(value ~ taxon + comparison, data = combined)
Anova(m, type = "II")
kruskal.test(value ~ taxcomp, data = combined)
kwAllPairsNemenyiTest(combined$value, as.factor(combined$taxcomp))

# Basic barplot
label.df <- data.frame(taxon = c("Fungi", "Fungi",
                                 "Bacteria", "Bacteria",
                                 "Archaea", "Archaea"),
                       comparison = c("between","within",
                                      "between","within",
                                      "between","within"),
                       Value = c(1.05,1.05,1.05,1.05,1.05,1.05),
                       Sig = c("e","f","c","d","a","b"))
label.B <- data.frame(taxon = c("Archaea"),
                       comparison = c("between"),
                       Value = c(1.05),
                       Sig = c("B"))
ggplot(data = combined, aes(comparison, value, colour = taxcomp)) +
  geom_jitter(size = 0.5, alpha = 0.01) +
  geom_boxplot(aes(comparison, value), outlier.shape = NA, color = "black", 
               fill = NA, inherit.aes = F) +
  geom_text(data = label.df, aes(x = comparison, y = Value, label = Sig, group = NULL),
            inherit.aes = F) +
  scale_colour_manual(values = brewer.pal(6, "Paired")[1:6]) +
  labs(x = "Environment comparison",
       y = "Bray-Curtis dissimilarity") +
  facet_wrap(~ taxon, ncol = 3) +
  ylim(0,1.05) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        legend.position = "none")

# With density plot
# Can't do with facet wrap, so plot separately and use plot_grid
p1 <- ggplot(data = arc_bray_df_long, aes(comparison, value, colour = comparison)) +
  geom_jitter(size = 0.5, alpha = 0.01) +
  geom_boxplot(aes(comparison, value), outlier.shape = NA, color = "black", 
               fill = NA, inherit.aes = F) +
  geom_text(data = label.df, aes(x = comparison, y = Value, label = Sig, group = NULL),
            inherit.aes = F) +
  scale_colour_manual(values = brewer.pal(6, "Paired")[1:2]) +
  labs(x = "Environment comparison",
       y = "Bray-Curtis dissimilarity") +
  ylim(0,1.05) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        legend.position = "none")
p2 <- ggdensity(arc_bray_df_long, "value", fill = "comparison", 
                   palette = c(brewer.pal(6, "Paired")[1:2]), size = 0.25) +
  rotate() + 
  clean_theme() + 
  rremove("legend") + 
  xlim(0.55, 1.05)
p_arc <- insert_yaxis_grob(p1, p2, position = "right")

p3 <- ggplot(data = bac_bray_df_long, aes(comparison, value, colour = comparison)) +
  geom_jitter(size = 0.5, alpha = 0.01) +
  geom_boxplot(aes(comparison, value), outlier.shape = NA, color = "black", 
               fill = NA, inherit.aes = F) +
  geom_text(data = label.df, aes(x = comparison, y = Value, label = Sig, group = NULL),
            inherit.aes = F) +
  scale_colour_manual(values = brewer.pal(6, "Paired")[3:4]) +
  labs(x = "Environment comparison",
       y = "Bray-Curtis dissimilarity") +
  ylim(0,1.05) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank())
p4 <- ggdensity(bac_bray_df_long, "value", fill = "comparison", 
                palette = c(brewer.pal(6, "Paired")[3:4]), size = 0.25) +
  rotate() + 
  clean_theme() + 
  rremove("legend") + 
  xlim(0.55, 1.05)
p_bac <- insert_yaxis_grob(p3, p4, position = "right")

p5 <- ggplot(data = fun_bray_df_long, aes(comparison, value, colour = comparison)) +
  geom_jitter(size = 0.5, alpha = 0.01) +
  geom_boxplot(aes(comparison, value), outlier.shape = NA, color = "black", 
               fill = NA, inherit.aes = F) +
  geom_text(data = label.df, aes(x = comparison, y = Value, label = Sig, group = NULL),
            inherit.aes = F) +
  scale_colour_manual(values = brewer.pal(6, "Paired")[5:6]) +
  labs(x = "Environment comparison",
       y = "Bray-Curtis dissimilarity") +
  ylim(0,1.05) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank())
p6 <- ggdensity(fun_bray_df_long, "value", fill = "comparison", 
                palette = c(brewer.pal(6, "Paired")[5:6]), size = 0.25) +
  rotate() + 
  clean_theme() + 
  rremove("legend") + 
  xlim(0.55, 1.05)
p_fun <- insert_yaxis_grob(p5, p6, position = "right")

g1 <- ggdraw(p_arc)
g2 <- ggdraw(p_bac)
g3 <- ggdraw(p_fun)

plot_grid(g1, g2, g3,
          ncol = 3,
          labels = c("a) Archaea", "b) Bacteria", "c) Fungi"),
          rel_widths = c(0.4, 0.3, 0.3))

#### ___Figure 3 ####
# Best option is probably violin plot!
png("FigsUpdated2/CompareBrayCurtis_genus.png", width = 7, height = 3, units = "in", res = 300)
f3b <- ggplot(data = combined, aes(comparison, value)) +
  #geom_jitter(size = 0.5, alpha = 0.01, aes(colour = taxcomp)) + # don't show, too many
  geom_violin() +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2) +
  geom_text(data = label.df, aes(x = comparison, y = Value, label = Sig, group = NULL),
            inherit.aes = F) +
  geom_text(data = label.B,
            aes(x = -Inf, y = Inf, label = Sig, group = NULL, hjust = -0.5, 
                vjust = 2, fontface = "bold"), 
            check_overlap = T, size = 5, inherit.aes = F) +
  #scale_colour_manual(values = brewer.pal(6, "Paired")[1:6]) +
  labs(x = "Environment comparison",
       y = "Bray-Curtis dissimilarity") +
  facet_wrap(~ taxon, ncol = 3) +
  ylim(0,1.05) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"))
f3b
dev.off()

# Make multipanel Figure 3
# Get fungal PCoA, genus level, no title, but A label
# Remake first one with fungi title
input_fungi_CPM_nz$map_loaded$Axis01 <- scores(pcoa)[,1]
input_fungi_CPM_nz$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(input_fungi_CPM_nz$map_loaded, c("Environment"), find_hull)
g <- ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, Axis02)) +
  #geom_polygon(data = micro.hulls, 
  #             aes(colour = Environment, fill = Environment),
  #             alpha = 0.1, show.legend = F, size = 0.25) +
  geom_point(size = 2, alpha = 0.75, shape = 16, aes(colour = Environment),
             show.legend = T) +
  geom_text(aes(x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 2, fontface = "bold"),
            size = 5, check_overlap = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_colour_manual(values = color_mapping) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, vjust = -1),
        panel.grid = element_blank())
g

f3 <- plot_grid(g, f3b, ncol = 1, rel_heights = c(0.5, 0.5))
f3

png("FinalFigs/Figure3.png", width = 7, height = 6, units = "in", res = 300)
f3
dev.off()

combined %>%
  group_by(taxon, comparison) %>%
  summarise(mean = mean(value),
            se = se(value))



#### __Phylum ####
# Wrote function to streamline this some
arc_bray_df_long <- compareBC(input = input_arc_CPM_nz,
                              sampleID = "sampleID",
                              variable = "Environment",
                              level = 2)
bac_bray_df_long <- compareBC(input = input_bac_CPM,
                              sampleID = "sampleID",
                              variable = "Environment",
                              level = 2)
fun_bray_df_long <- compareBC(input = input_fungi_CPM_nz,
                              sampleID = "sampleID",
                              variable = "Environment",
                              level = 2)
# Combine and Graph
fun_bray_df_long$taxon <- "Fungi"
bac_bray_df_long$taxon <- "Bacteria"
arc_bray_df_long$taxon <- "Archaea"
combined <- rbind(fun_bray_df_long, bac_bray_df_long, arc_bray_df_long)
combined$taxcomp <- paste(combined$taxon, combined$comparison, sep = "_")
leveneTest(value ~ taxcomp, data = combined)
m <- aov(value ~ taxcomp, data = combined)
summary(m)
TukeyHSD(m)
shapiro.test(m$residuals[1:5000])
kruskal.test(value ~ taxcomp, data = combined)
kwAllPairsNemenyiTest(combined$value, as.factor(combined$taxcomp))

label.df <- data.frame(taxon = c("Fungi", "Fungi",
                                 "Bacteria", "Bacteria",
                                 "Archaea", "Archaea"),
                       comparison = c("between","within",
                                      "between","within",
                                      "between","within"),
                       Value = c(1.05,1.05,1.05,1.05,1.05,1.05),
                       Sig = c("c","f","c","d","a","b"))

png("FigsUpdated2/CompareBrayCurtis_phylum.png", width = 7, height = 3, units = "in", res = 300)
ggplot(data = combined, aes(comparison, value, colour = taxcomp)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2) +
  geom_text(data = label.df, aes(x = comparison, y = Value, label = Sig, group = NULL),
            inherit.aes = F) +
  scale_colour_manual(values = brewer.pal(6, "Paired")[1:6]) +
  labs(x = "Environment comparison",
       y = "Bray-Curtis dissimilarity") +
  facet_wrap(~ taxon, ncol = 3) +
  ylim(0,1.05) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        legend.position = "none")
dev.off()



#### __Class ####
# Wrote function to streamline this some
arc_bray_df_long <- compareBC(input = input_arc_CPM_nz,
                              sampleID = "sampleID",
                              variable = "Environment",
                              level = 3)
bac_bray_df_long <- compareBC(input = input_bac_CPM,
                              sampleID = "sampleID",
                              variable = "Environment",
                              level = 3)
fun_bray_df_long <- compareBC(input = input_fungi_CPM_nz,
                              sampleID = "sampleID",
                              variable = "Environment",
                              level = 3)
# Combine and Graph
fun_bray_df_long$taxon <- "Fungi"
bac_bray_df_long$taxon <- "Bacteria"
arc_bray_df_long$taxon <- "Archaea"
combined <- rbind(fun_bray_df_long, bac_bray_df_long, arc_bray_df_long)
combined$taxcomp <- paste(combined$taxon, combined$comparison, sep = "_")
leveneTest(value ~ taxcomp, data = combined)
m <- aov(value ~ taxcomp, data = combined)
summary(m)
TukeyHSD(m)
shapiro.test(m$residuals[1:5000])
kruskal.test(value ~ taxcomp, data = combined)
kwAllPairsNemenyiTest(combined$value, as.factor(combined$taxcomp))

png("FigsUpdated2/CompareBrayCurtis_class.png", width = 7, height = 3, units = "in", res = 300)
ggplot(data = combined, aes(comparison, value, colour = taxcomp)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2) +
  geom_text(data = label.df, aes(x = comparison, y = Value, label = Sig, group = NULL),
            inherit.aes = F) +
  scale_colour_manual(values = brewer.pal(6, "Paired")[1:6]) +
  labs(x = "Environment comparison",
       y = "Bray-Curtis dissimilarity") +
  facet_wrap(~ taxon, ncol = 3) +
  ylim(0,1.05) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        legend.position = "none")
dev.off()



#### __Order ####
# Wrote function to streamline this some
arc_bray_df_long <- compareBC(input = input_arc_CPM_nz,
                              sampleID = "sampleID",
                              variable = "Environment",
                              level = 4)
bac_bray_df_long <- compareBC(input = input_bac_CPM,
                              sampleID = "sampleID",
                              variable = "Environment",
                              level = 4)
fun_bray_df_long <- compareBC(input = input_fungi_CPM_nz,
                              sampleID = "sampleID",
                              variable = "Environment",
                              level = 4)
# Combine and Graph
fun_bray_df_long$taxon <- "Fungi"
bac_bray_df_long$taxon <- "Bacteria"
arc_bray_df_long$taxon <- "Archaea"
combined <- rbind(fun_bray_df_long, bac_bray_df_long, arc_bray_df_long)
combined$taxcomp <- paste(combined$taxon, combined$comparison, sep = "_")
leveneTest(value ~ taxcomp, data = combined)
m <- aov(value ~ taxcomp, data = combined)
summary(m)
TukeyHSD(m)
shapiro.test(m$residuals[1:5000])
kruskal.test(value ~ taxcomp, data = combined)
kwAllPairsNemenyiTest(combined$value, as.factor(combined$taxcomp))

png("FigsUpdated2/CompareBrayCurtis_order.png", width = 7, height = 3, units = "in", res = 300)
ggplot(data = combined, aes(comparison, value, colour = taxcomp)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2) +
  geom_text(data = label.df, aes(x = comparison, y = Value, label = Sig, group = NULL),
            inherit.aes = F) +
  scale_colour_manual(values = brewer.pal(6, "Paired")[1:6]) +
  labs(x = "Environment comparison",
       y = "Bray-Curtis dissimilarity") +
  facet_wrap(~ taxon, ncol = 3) +
  ylim(0,1.05) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        legend.position = "none")
dev.off()



#### __Family ####
# Wrote function to streamline this some
arc_bray_df_long <- compareBC(input = input_arc_CPM_nz,
                              sampleID = "sampleID",
                              variable = "Environment",
                              level = 5)
bac_bray_df_long <- compareBC(input = input_bac_CPM,
                              sampleID = "sampleID",
                              variable = "Environment",
                              level = 5)
fun_bray_df_long <- compareBC(input = input_fungi_CPM_nz,
                              sampleID = "sampleID",
                              variable = "Environment",
                              level = 5)
# Combine and Graph
fun_bray_df_long$taxon <- "Fungi"
bac_bray_df_long$taxon <- "Bacteria"
arc_bray_df_long$taxon <- "Archaea"
combined <- rbind(fun_bray_df_long, bac_bray_df_long, arc_bray_df_long)
combined$taxcomp <- paste(combined$taxon, combined$comparison, sep = "_")
leveneTest(value ~ taxcomp, data = combined)
m <- aov(value ~ taxcomp, data = combined)
summary(m)
TukeyHSD(m)
shapiro.test(m$residuals[1:5000])
kruskal.test(value ~ taxcomp, data = combined)
kwAllPairsNemenyiTest(combined$value, as.factor(combined$taxcomp))

png("FigsUpdated2/CompareBrayCurtis_family.png", width = 7, height = 3, units = "in", res = 300)
ggplot(data = combined, aes(comparison, value, colour = taxcomp)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2) +
  geom_text(data = label.df, aes(x = comparison, y = Value, label = Sig, group = NULL),
            inherit.aes = F) +
  scale_colour_manual(values = brewer.pal(6, "Paired")[1:6]) +
  labs(x = "Environment comparison",
       y = "Bray-Curtis dissimilarity") +
  facet_wrap(~ taxon, ncol = 3) +
  ylim(0,1.05) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        legend.position = "none")
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
adonis2(bc ~ Year, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.18, p = 0.001
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
  ggtitle("Year: R2 = 0.18") +
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
table(input_fungi_CPM_nz$map_loaded$`Assembly.Method`)
input_fungi_CPM_nz$map_loaded$Assembler <- dplyr::recode_factor(input_fungi_CPM_nz$map_loaded$`Assembly.Method`,
                                   "AbySS v1.5.0" = "AbySS",
                                   "canu v. 1.7" = "Canu",
                                   "canu v. 1.8" = "Canu",
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
                                   "JGI custom assembly Copeland et. al" = "Custom JGI",
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
                                   "MEGAHIT v. 0.2.0" = "MEGAHIT",
                                   "MegaHit v. 1.02" = "MEGAHIT",
                                   "MEGAHIT v. 1.0.3" = "MEGAHIT",
                                   "MEGAHIT v. 1.1.1" = "MEGAHIT",
                                   "megahit v. 1.1.3" = "MEGAHIT",
                                   "Megahit v. 1.1.3" = "MEGAHIT",
                                   "MegaHIT v. 1.2.9" = "MEGAHIT",
                                   "MEGAHIT v. MEGAHIT v0.2.0" = "MEGAHIT",
                                   "MEGAHIT v. MEGAHIT v1.0.3" = "MEGAHIT",
                                   "MEGAHIT v. MEGAHIT v1.0.6" = "MEGAHIT",
                                   "MEGAHIT v. 1.0.6" = "MEGAHIT",
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
                                   #"" = "Unknown",
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
adonis2(bc ~ Assembler, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.18, p = 0.001
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
  ggtitle("Assembler: R2 = 0.18") +
  guides(colour = guide_legend(override.aes = list(alpha = 1),
                               ncol = 1)) +
  theme_bw() +  
  theme(legend.position = "right",
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))
g_assem

pdf("FigsUpdated2/PCoA_Year_Assembler.pdf", width = 9, height = 4)
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

pdf("FigsUpdated2/SampleSize_Year_Assembler.pdf", width = 9, height = 4)
plot_grid(g_numyear, g_numassem, align = "hv", ncol = 2)
dev.off()



#### _by Ecosystem ####
# Subset the data into each environment and study
# For glacier forefield, it's one study but 4 locations so separate by location
# ^ Now added more GF studies
# Make pie charts of taxa
# Use fungi only input data - get relative abundances of fungal phyla out of just fungi
# Use samples with at least 1 fungal count

# Need to do several for loops
# For subsetting, summarizing, plotting, store dfs in a list to enable for loop/indexing
env <- list()
studyN <- list()
df <- list()
phy <- list()
p_colors <- list()
p <- list()

# Subset by Environment (9 data frames)
for (i in 1:length(levels(input_fungi_nz$map_loaded$Environment))) {
  env[[i]] <- filter_data(input_fungi_nz,
                          filter_cat = "Environment",
                          keep_vals = levels(input_fungi_nz$map_loaded$Environment)[i])
}

# Do one for each, or by individual studies? How many studies in each env?
length(levels(env[[1]]$map_loaded$`Study.Name2`)) # 4 AMD
length(levels(env[[2]]$map_loaded$`Study.Name2`)) # 9 Cryo soil
length(levels(env[[3]]$map_loaded$`Study.Name2`)) # 8 Cryo water
length(levels(env[[4]]$map_loaded$`Study.Name2`)) # 5 Desert
length(levels(env[[5]]$map_loaded$`Study.Name2`)) # 5 Glacial forefield
length(levels(env[[6]]$map_loaded$`Study.Name2`)) # 18 Hot spring
length(levels(env[[7]]$map_loaded$`Study.Name2`)) # 39 Hydro
length(levels(env[[8]]$map_loaded$`Study.Name2`)) # 17 Hypersaline
length(levels(env[[9]]$map_loaded$`Study.Name2`)) # 2 Soda Lake

# Probably best to show some of the variability within environment type
# For example, Shu and Huang 2022 have 4-5 sites for each environment
# Here let's do 1-4, for env. with more than 4 studies, take 4 with greatest sample size
# Could also decide to show sites with most fungi instead of most samples
for (i in 1:length(env)) {
  studyN[[i]] <- as.data.frame(table(env[[i]]$map_loaded$`Study.Name2`)) %>%
    arrange(desc(Freq)) %>%
    slice_head(n = 4)
}

studies <- rbind(studyN[[1]], studyN[[2]], studyN[[3]], studyN[[4]],
                 studyN[[5]], studyN[[6]], studyN[[7]], studyN[[8]],
                 studyN[[9]]) %>%
  mutate(Var1 = as.character(Var1))

# Subset Environments by Studies
counter <- 1
for (i in 1:length(env)) {
  k <- nrow(studyN[[i]])
  for (l in 1:k) {
  df[[counter]] <- filter_data(env[[i]],
                               filter_cat = "Study.Name2",
                               keep_vals = studyN[[i]]$Var1[l])
  counter <- counter + 1
  }
}

# Location and n (for plot titles)
loc_n <- as.data.frame(matrix(NA, nrow = length(df), ncol = 3)) %>%
  set_names(c("Location", "n", "Environment"))
for (i in 1:length(df)) {
  loc_n$Location[i] <- levels(df[[i]]$map_loaded$`Location2`)[1]
  loc_n$n[i] <- nrow(df[[i]]$map_loaded)
  loc_n$Environment[i] <- levels(df[[i]]$map_loaded$Environment)[1]
}
for (i in 1:length(df)) {
  df[[i]]$map_loaded$Location <- paste(loc_n$Location[i],"\n", "(n = ", loc_n$n[i], ")", sep = "")
}
# Manually update some long ones if needed
df[[1]]$map_loaded$Location <- "USA: Iron Mountain, CA\n(n = 17)"
df[[4]]$map_loaded$Location <- "USA: Akron, OH\n(n = 1)"
df[[12]]$map_loaded$Location <- "Antarctica: Canada Glacier\n(n = 3)"
df[[14]]$map_loaded$Location <- "USA: Jornada LTER, NM\n(n = 18)"
df[[16]]$map_loaded$Location <- "USA: Green Butte, UT\n(n = 7)"
df[[17]]$map_loaded$Location <- "Russell Glacier, Greenland\n(n = 24)"
df[[18]]$map_loaded$Location <- "Midre Lovenbreen, Norway\n(n = 12)"
df[[19]]$map_loaded$Location <- "Storglaciaren, Sweden\n(n = 8)"
df[[20]]$map_loaded$Location <- "Rabots glacier, Sweden\n(n = 4)"
df[[28]]$map_loaded$Location <- "Anemone Vent, Axial Seamount\n(n = 12)"
df[[29]]$map_loaded$Location <- "Garden Lake, Australia\n(n = 117)"
df[[30]]$map_loaded$Location <- "Eisfeld solar saltern, Namibia\n(n = 7)"



#### __Phyla ####
# N.B. summarize taxonomy does not work properly with 1 sample
# For those need to transpose the dataframe
for (i in 1:nrow(studies)) {
  phy[[i]] <- summarize_taxonomy(df[[i]], level = 2, report_higher_tax = F, relative = T)
}

for (i in 1:nrow(studies)) {
  if (nrow(phy[[i]]) == 1) {
    phy[[i]] <- as.data.frame(t(phy[[i]]))
    colnames(phy[[i]]) <- rownames(df[[i]]$map_loaded)
  }
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
# Need to do it differently for n > 1 and n = 1
for (i in 1:nrow(studies)) {
  t <- df[[i]]$map_loaded$Location
  p[[i]] <- plot_taxa_bars(phy[[i]], df[[i]]$map_loaded, "Study.Name2", 20) +
    coord_polar("y", start=0) +
    scale_fill_manual(values = p_colors[[i]]) +
    theme_void() +
    ggtitle(t) +
    theme(legend.position = "none",
          plot.title = element_text(size = 5, hjust = 0.5, vjust = -5),
          plot.margin = margin(0, -50, -15, -50, "pt"))
}

# Get legend - use one with all 7 phyla
t <- df[[4]]$map_loaded$Location
p_forleg <- plot_taxa_bars(phy[[4]], df[[4]]$map_loaded, "Study.Name2", 20) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start=0) +
  scale_fill_manual(values = p_colors[[4]]) +
  labs(fill = "Phylum") +
  theme_void() +
  ggtitle(t) +
  theme(legend.position = "right",
        # legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(size = 10, hjust = 0.5, vjust = -15),
        plot.margin = margin(0,0,0,0, "pt"))
p_forleg
pie_leg <- get_legend(p_forleg)

# Troubleshoot
# ptb <- data.frame("group_by" = df[[4]]$map_loaded$Study.Name2,
#                   "taxon" = rownames(phy[[4]]),
#                   "mean_value" = phy[[4]][,1])
# ggplot(ptb, aes(x = group_by, y = mean_value, fill = taxon)) +
#   geom_bar(stat = "identity", width = 1, color = "white") +
#   coord_polar("y", start=0) +
#   scale_fill_manual(values = p_colors[[4]]) +
#   labs(fill = "Phylum") +
#   theme_void() +
#   ggtitle(t) +
#   theme(legend.position = "right",
#         # legend.key.size = unit(0.25, "cm"),
#         plot.title = element_text(size = 10, hjust = 0.5, vjust = -15),
#         plot.margin = margin(0,0,0,0, "pt"))
  

# Make huge multipanel
table(loc_n$Environment)
# 4, 4, 4, 4, 4, 4, 4, 4, 2
pies <- plot_grid(p[[1]], p[[2]], p[[3]], p[[4]],
                  p[[5]], p[[6]], p[[7]], p[[8]],
                  p[[9]], p[[10]], p[[11]], p[[12]],
                  p[[13]], p[[14]], p[[15]], p[[16]],
                  p[[17]], p[[18]], p[[19]], p[[20]],
                  p[[21]], p[[22]], p[[23]], p[[24]],
                  p[[25]], p[[26]], p[[27]], p[[28]],
                  p[[29]], p[[30]], p[[31]], p[[32]],
                  p[[33]], p[[34]], NULL, NULL,
                  ncol = 4)
pies

# Add labels
l1 <- ggplot(data = NULL, aes(x = 1, y = 1)) +
  geom_text(aes(label = "Acid mine\ndrainage", angle = 0, fontface = "bold"), size = 2.5) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "pt"))
l2 <- ggplot(data = NULL, aes(x = 1, y = 1)) +
  geom_text(aes(label = "Cryosphere\nsoil", angle = 0, fontface = "bold"), size = 2.5) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "pt"))
l3 <- ggplot(data = NULL, aes(x = 1, y = 1)) +
  geom_text(aes(label = "Cryosphere\nwater", angle = 0, fontface = "bold"), size = 2.5) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "pt"))
l4 <- ggplot(data = NULL, aes(x = 1, y = 1)) +
  geom_text(aes(label = "Desert", angle = 0, fontface = "bold"), size = 2.5) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "pt"))
l5 <- ggplot(data = NULL, aes(x = 1, y = 1)) +
  geom_text(aes(label = "Glacial\nforefield", angle = 0, fontface = "bold"), size = 2.5) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "pt"))
l6 <- ggplot(data = NULL, aes(x = 1, y = 1)) +
  geom_text(aes(label = "Hot\nspring", angle = 0, fontface = "bold"), size = 2.5) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "pt"))
l7 <- ggplot(data = NULL, aes(x = 1, y = 1)) +
  geom_text(aes(label = "Hydrothermal\nvent", angle = 0, fontface = "bold"), size = 2.5) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "pt"))
l8 <- ggplot(data = NULL, aes(x = 1, y = 1)) +
  geom_text(aes(label = "Hypersaline", angle = 0, fontface = "bold"), size = 2.5) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "pt"))
l9 <- ggplot(data = NULL, aes(x = 1, y = 1)) +
  geom_text(aes(label = "Soda lake", angle = 0, fontface = "bold"), size = 2.5) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "pt"))

# Remake
pies <- plot_grid(l1, p[[1]], p[[2]], p[[3]], p[[4]],
                  l2, p[[5]], p[[6]], p[[7]], p[[8]],
                  l3, p[[9]], p[[10]], p[[11]], p[[12]],
                  l4, p[[13]], p[[14]], p[[15]], p[[16]],
                  l5, p[[17]], p[[18]], p[[19]], p[[20]],
                  l6, p[[21]], p[[22]], p[[23]], p[[24]],
                  l7, p[[25]], p[[26]], p[[27]], p[[28]],
                  l8, p[[29]], p[[30]], p[[31]], p[[32]],
                  l9, p[[33]], p[[34]], NULL, NULL, NULL,
                  ncol = 5,
                  rel_widths = c(0.04, 0.24, 0.24, 0.24, 0.24))
pies

# Add legend
pdf("FigsUpdated2/Pies_Phyla.pdf", width = 8, height = 10)
plot_grid(pies, pie_leg, rel_widths = c(4, 1))
dev.off()



#### __Class ####
cla <- list()
c_colors <- list()
c <- list()

# Get classes with max abundance cutoff
# 37 classes total, need to get this down to 24 max
input_fungi_nz
colMax <- function(data) apply(data, 1, max, na.rm = TRUE)
class_summary <- summarize_taxonomy(input_fungi_nz, level = 3,
                                    report_higher_tax = F, relative = T) %>%
  mutate(maxes = colMax(.))

class_top <- summarize_taxonomy(input_fungi_nz, level = 3,
                                   report_higher_tax = F, relative = T) %>%
  mutate(maxes = colMax(.)) %>%
  filter(maxes > 0.25)
top_classes <- rownames(class_top)
class_bottom <- summarize_taxonomy(input_fungi_nz, level = 3,
                                report_higher_tax = F, relative = T) %>%
  mutate(maxes = colMax(.)) %>%
  filter(maxes <= 0.25) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Other = rowSums(.)) %>%
  t() %>%
  as.data.frame()

class_top_wO <- summarize_taxonomy(input_fungi_nz, level = 3,
                                report_higher_tax = F, relative = T) %>%
  mutate(maxes = colMax(.)) %>%
  filter(maxes > 0.25) %>%
  add_row(class_bottom[18,])
rownames(class_top)[nrow(class_top)] <- "Other"

# Classes
# Studies with 1 sample need to be treated differently
for (i in 1:nrow(studies)) {
  if (nrow(df[[i]]$map_loaded) == 1) {
    csb <- summarize_taxonomy(df[[i]], level = 3, report_higher_tax = F, relative = T) %>%
      t() %>%
      as.data.frame() %>%
      set_names(rownames(df[[i]]$map_loaded)) %>%
      filter(rownames(.) %notin% top_classes) %>%
      t() %>%
      as.data.frame() %>%
      mutate(Other = rowSums(.)) %>%
      t() %>%
      as.data.frame() %>%
      filter(rownames(.) == "Other")
    cst <- summarize_taxonomy(df[[i]], level = 3, report_higher_tax = F, relative = T) %>%
      t() %>%
      as.data.frame() %>%
      set_names(rownames(df[[i]]$map_loaded)) %>%
      filter(rownames(.) %in% top_classes)
    cst <- rbind(cst,csb)
    cla[[i]] <- cst
  }
  
  if (nrow(df[[i]]$map_loaded > 1)) {
   csb <- summarize_taxonomy(df[[i]], level = 3, report_higher_tax = F, relative = T) %>%
            filter(rownames(.) %notin% top_classes) %>%
            t() %>%
            as.data.frame() %>%
            mutate(Other = rowSums(.)) %>%
            t() %>%
            as.data.frame()
    cst <- summarize_taxonomy(df[[i]], level = 3, report_higher_tax = F, relative = T) %>%
      filter(rownames(.) %in% top_classes) %>%
      add_row(csb[nrow(csb),])
    rownames(cst)[nrow(cst)] <- "Other"
    cla[[i]] <- cst
  }
}

# Colors
classes_present <- c(top_classes[1:18],"Ustilaginomycetes","unclassified","Other")
for (i in 1:length(df)) {
  c_colors[[i]] <- data.frame("class" = classes_present,
                "color" = c(colorRampPalette(brewer.pal(12, "Paired"))(length(classes_present)-2),
                                          "grey75", "grey90"))
  c_colors[[i]] <- subset(c_colors[[i]],
                          class %in% rownames(cla[[i]]))
}

# Pies
for (i in 1:nrow(studies)) {
  t <- df[[i]]$map_loaded$Location
  c[[i]] <- plot_taxa_bars(cla[[i]], df[[i]]$map_loaded, "Study.Name2", 37) +
    coord_polar("y", start=0) +
    scale_fill_manual(values = c_colors[[i]]$color) +
    theme_void() +
    ggtitle(t) +
    theme(legend.position = "none",
          plot.title = element_text(size = 5, hjust = 0.5, vjust = -5),
          plot.margin = margin(0, -50, -15, -50, "pt"))
}

# Get legend - top 20 (unclassified last), plus Other
t <- df[[34]]$map_loaded$Location
c_forleg <- plot_taxa_bars(cla[[34]], df[[34]]$map_loaded, "Study.Name2", 37) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start=0) +
  scale_fill_manual(values = c_colors[[34]]$color) +
  labs(fill = "Classes") +
  theme_void() +
  ggtitle(t) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(legend.position = "right",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 10, hjust = 0.5, vjust = -15),
        plot.margin = margin(0,0,0,0, "pt"))
c_forleg
pie_leg_c <- get_legend(c_forleg)

# Make huge multipanel
pies_c <- plot_grid(l1, c[[1]], c[[2]], c[[3]], c[[4]],
                    l2, c[[5]], c[[6]], c[[7]], c[[8]],
                    l3, c[[9]], c[[10]], c[[11]], c[[12]],
                    l4, c[[13]], c[[14]], c[[15]], c[[16]],
                    l5, c[[17]], c[[18]], c[[19]], c[[20]],
                    l6, c[[21]], c[[22]], c[[23]], c[[24]],
                    l7, c[[25]], c[[26]], c[[27]], c[[28]],
                    l8, c[[29]], c[[30]], c[[31]], c[[32]],
                    l9, c[[33]], c[[34]], NULL, NULL, NULL,
                    ncol = 5,
                    rel_widths = c(0.04, 0.24, 0.24, 0.24, 0.24))
pies_c

# No text and no n = 1
pies_c <- plot_grid(c[[1]], c[[2]], NULL, NULL,
                    c[[5]], c[[6]], c[[7]], c[[8]],
                    c[[9]], c[[10]], c[[11]], c[[12]],
                    c[[13]], c[[14]], c[[15]], c[[16]],
                    c[[17]], c[[18]], c[[19]], c[[20]],
                    c[[21]], c[[22]], c[[23]], c[[24]],
                    c[[25]], c[[26]], c[[27]], c[[28]],
                    c[[29]], c[[30]], c[[31]], c[[32]],
                    c[[33]], c[[34]], NULL, NULL, NULL,
                    ncol = 4)
pies_c

# Add legend
pdf("FigsUpdated2/Pies_Classes.pdf", width = 8, height = 10)
plot_grid(pies_c, pie_leg_c, rel_widths = c(3.85, 1.15))
dev.off()

#### ___Figure 4 ####
png("FinalFigs/Figure4forPPT.png", width = 7, height = 9, units = "in", res = 300)
plot_grid(pies_c, pie_leg_c, rel_widths = c(3.73, 1.27))
dev.off()



#### _Top Taxa ####
# Barplots were for top taxa overall, but let's now get the top taxa in each environment
# Use input file for each environment, and summarize with mctoolsr at each taxonomic level
# Combined into data frame and write file so people can see
# Already got 9 input data lists in the "by Ecosystem" section
# These are  CPM, and don't contain any zeroes
# Environments
env <- list()
for (i in 1:length(levels(input_fungi_CPM_nz$map_loaded$Environment))) {
  env[[i]] <- filter_data(input_fungi_CPM_nz,
                          filter_cat = "Environment",
                          keep_vals = levels(input_fungi_CPM_nz$map_loaded$Environment)[i])
}


# Summary taxonomy
phy <- list()
for (i in 1:length(levels(input_fungi_CPM_nz$map_loaded$Environment))) {
  phy[[i]] <- summarize_taxonomy(env[[i]], level = 2, report_higher_tax = F, relative = F)
}
cla <- list()
for (i in 1:length(levels(input_fungi_CPM_nz$map_loaded$Environment))) {
  cla[[i]] <- summarize_taxonomy(env[[i]], level = 3, report_higher_tax = F, relative = F)
}
ord <- list()
for (i in 1:length(levels(input_fungi_CPM_nz$map_loaded$Environment))) {
  ord[[i]] <- summarize_taxonomy(env[[i]], level = 4, report_higher_tax = F, relative = F)
}
fam <- list()
for (i in 1:length(levels(input_fungi_CPM_nz$map_loaded$Environment))) {
  fam[[i]] <- summarize_taxonomy(env[[i]], level = 5, report_higher_tax = F, relative = F)
}
gen <- list()
for (i in 1:length(levels(input_fungi_CPM_nz$map_loaded$Environment))) {
  gen[[i]] <- summarize_taxonomy(env[[i]], level = 6, report_higher_tax = F, relative = F)
}

# Top n taxa
ntax = 10 # set n to 10

top_phy <- list()
for (i in 1:length(levels(input_fungi_CPM_nz$map_loaded$Environment))) {
top_phy[[i]] <- plot_taxa_bars(phy[[i]],
                               env[[i]]$map_loaded,
                               type_header = "Environment",
                               num_taxa = ntax,
                               data_only = T) %>%
  mutate(Level = "Phylum") %>%
  arrange(desc(mean_value))
}

top_cla <- list()
for (i in 1:length(levels(input_fungi_CPM_nz$map_loaded$Environment))) {
  top_cla[[i]] <- plot_taxa_bars(cla[[i]],
                                 env[[i]]$map_loaded,
                                 type_header = "Environment",
                                 num_taxa = ntax,
                                 data_only = T) %>%
    mutate(Level = "Class") %>%
    arrange(desc(mean_value))
}

top_ord <- list()
for (i in 1:length(levels(input_fungi_CPM_nz$map_loaded$Environment))) {
  top_ord[[i]] <- plot_taxa_bars(ord[[i]],
                                 env[[i]]$map_loaded,
                                 type_header = "Environment",
                                 num_taxa = ntax,
                                 data_only = T) %>%
    mutate(Level = "Order") %>%
    arrange(desc(mean_value))
}

top_fam <- list()
for (i in 1:length(levels(input_fungi_CPM_nz$map_loaded$Environment))) {
  top_fam[[i]] <- plot_taxa_bars(fam[[i]],
                                 env[[i]]$map_loaded,
                                 type_header = "Environment",
                                 num_taxa = ntax,
                                 data_only = T) %>%
    mutate(Level = "Family") %>%
    arrange(desc(mean_value))
}

top_gen <- list()
for (i in 1:length(levels(input_fungi_CPM_nz$map_loaded$Environment))) {
  top_gen[[i]] <- plot_taxa_bars(gen[[i]],
                                 env[[i]]$map_loaded,
                                 type_header = "Environment",
                                 num_taxa = ntax,
                                 data_only = T) %>%
    mutate(Level = "Genus") %>%
    arrange(desc(mean_value))
}

#### ___Table S4 ####
# Combine
top_taxa <- rbind(top_phy[[1]], top_phy[[2]], top_phy[[3]],
                  top_phy[[4]], top_phy[[5]], top_phy[[6]],
                  top_phy[[7]], top_phy[[8]], top_phy[[9]],
                  top_cla[[1]], top_cla[[2]], top_cla[[3]],
                  top_cla[[4]], top_cla[[5]], top_cla[[6]],
                  top_cla[[7]], top_cla[[8]], top_cla[[9]],
                  top_ord[[1]], top_ord[[2]], top_ord[[3]],
                  top_ord[[4]], top_ord[[5]], top_ord[[6]],
                  top_ord[[7]], top_ord[[8]], top_ord[[9]],
                  top_fam[[1]], top_fam[[2]], top_fam[[3]],
                  top_fam[[4]], top_fam[[5]], top_fam[[6]],
                  top_fam[[7]], top_fam[[8]], top_fam[[9]],
                  top_gen[[1]], top_gen[[2]], top_gen[[3]],
                  top_gen[[4]], top_gen[[5]], top_gen[[6]],
                  top_gen[[7]], top_gen[[8]], top_gen[[9]]) %>%
  filter(taxon != "Other") %>%
  rename("Environment" = "group_by") %>%
  rename("mean_CPM" = "mean_value") %>%
  mutate(mean_CPM = round(mean_CPM, digits = 2))
#write_xlsx(top_taxa, "data/top_taxa_by_env_final.xlsx", format_headers = F)
#write_xlsx(top_taxa, "data/TableS4.xlsx", format_headers = F)



#### _Alpha ####
# Look at number of different taxa levels by environments
# OTU Richness
input_fungi$map_loaded$rich <- specnumber(input_fungi$data_loaded, 
                                          MARGIN = 2)
max(input_fungi$map_loaded$rich)


# Shannon diversity
input_fungi$map_loaded$shannon <- diversity(input_fungi$data_loaded, 
                                            index = "shannon", 
                                            MARGIN = 2)
max(input_fungi$map_loaded$shannon)

# Stats and graphs
hist(input_fungi$map_loaded$rich)
hist(log(input_fungi$map_loaded$rich + 1))
leveneTest(input_fungi$map_loaded$rich ~ input_fungi$map_loaded$Environment) # Bad
m3 <- aov(rich ~ Environment, data = input_fungi$map_loaded)
summary(m3)
TukeyHSD(m3)
shapiro.test(m3$residuals) # Bad
hist(m3$residuals) # Could argue it approximates normal
plot(m3$fitted.values, m3$residuals)
plot(m3)
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

hist(input_fungi$map_loaded$shannon)
hist(log(input_fungi$map_loaded$shannon + 1))
leveneTest(input_fungi$map_loaded$shannon ~ input_fungi$map_loaded$Environment) # Bad
m4 <- aov(shannon ~ Environment, data = input_fungi$map_loaded)
shapiro.test(m4$residuals) # Bad
summary(m4)
TukeyHSD(m4)
kruskal.test(shannon ~ Environment, data = input_fungi$map_loaded)
nyi4 <- kwAllPairsNemenyiTest(shannon ~ Environment, data = input_fungi$map_loaded)
nyi_table4 <- fullPTable(nyi4$p.value)
nyi_list4 <- multcompLetters(nyi_table4)
nyi_let4 <- as.data.frame(nyi_list4$Letters) %>%
  mutate(label = `nyi_list4$Letters`,
         y = rep(5.5, nrow(.)),
         name = "shannon") %>%
  dplyr::select(-`nyi_list4$Letters`) %>%
  rownames_to_column(var = "Environment")

# Check just richness for ordering
ggplot(input_fungi$map_loaded, aes(reorder(Environment, rich, median), rich)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_violin() +
  geom_jitter(size = 0.5, alpha = 0.2, width = 0.3) +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 3, colour = "red") +
  #geom_text(data = label_df, aes(Environment, y, label = label), 
  #          size = 4, color = "black") +
  labs(x = NULL, y = NULL) +
  #facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        strip.text = element_text(size = 10),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        plot.margin = margin(1, 1, 1, 15))



#### ___Figure 2 ####
label_df <- rbind(nyi_let3, nyi_let4) %>%
  mutate(Environment = factor(Environment,
                              levels = c("Acid mine drainage",
                                         "Hypersaline",
                                         "Glacial forefield",
                                         "Hydrothermal vent",
                                         "Cryosphere - soil",
                                         "Hot spring",
                                         "Soda lake",
                                         "Cryosphere - water",
                                         "Desert")))
facet_df <- c("rich" = "(a) Genus richness",
              "shannon" = "(b) Genus Shannon")
alpha_long <- input_fungi$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon")) %>%
  mutate(Environment = factor(Environment,
                              levels = c("Acid mine drainage",
                                         "Hypersaline",
                                         "Glacial forefield",
                                         "Hydrothermal vent",
                                         "Cryosphere - soil",
                                         "Hot spring",
                                         "Soda lake",
                                         "Cryosphere - water",
                                         "Desert")))

pdf("FigsUpdated2/AlphaDiversity.pdf", width = 6, height = 3)
f2 <- ggplot(alpha_long, aes(Environment, value, colour = Environment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1, alpha = 0.25, shape = 16, width = 0.3) +
  stat_summary(fun.y = mean, geom = "point", shape = 23, 
               size = 4, colour = "black") +
  geom_text(data = label_df, aes(Environment, y, label = label), 
            size = 4, color = "black") +
  labs(x = NULL, y = NULL) +
  scale_colour_manual(values = color_mapping) +
  facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        strip.text = element_text(size = 10),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        plot.margin = margin(1, 1, 1, 15))
f2
dev.off()

png("FinalFigs/Figure2.png", width = 6, height = 3, units = "in", res = 300)
f2
dev.off()



#### _Drivers/Indicators ####
# MULTIPATT (list taxa associated with each group)
# N.B. updated plot_multipatt for use with other projects
# To use original function, call plot_multipatt_fungal

# Phyla
tax_sum_phyla <- summarize_taxonomy(input_fungi_CPM, level = 2, report_higher_tax = F, relative = F)
set.seed(425)
mp_phyla <- multipatt(t(tax_sum_phyla), 
                      input_fungi_CPM$map_loaded$Environment, 
                      func = "r.g", 
                      control = how(nperm=999))
summary(mp_phyla) # 3 soda lake

pdf("FigsUpdated2/mp_Phyla.pdf", width = 5, height = 5)
plot_multipatt_fungal(mp_obj = mp_phyla, 
                      input = input_fungi_CPM,
                      tax_sum = tax_sum_phyla,
                      group = "Environment",
                      qcut = 0.05,
                      rcut = 0)
dev.off()



# Class
tax_sum_class <- summarize_taxonomy(input_fungi_CPM, level = 3, report_higher_tax = F, relative = F)
set.seed(425)
mp_class <- multipatt(t(tax_sum_class), 
                      input_fungi_CPM$map_loaded$Environment, 
                      func = "r.g", 
                      control = how(nperm=999))
summary(mp_class) # 1 AMD, 1 Desert, 22 soda lake
pdf("FigsUpdated2/mp_Class.pdf", width = 5, height = 5)
plot_multipatt_fungal(mp_obj = mp_class, 
                      input = input_fungi_CPM,
                      tax_sum = tax_sum_class,
                      group = "Environment",
                      qcut = 0.05,
                      rcut = 0)
dev.off()

# Order
tax_sum_order <- summarize_taxonomy(input_fungi_CPM, level = 4, report_higher_tax = F, relative = F)
set.seed(425)
mp_order <- multipatt(t(tax_sum_order), 
                      input_fungi_CPM$map_loaded$Environment, 
                      func = "r.g", 
                      control = how(nperm=999))
summary(mp_order) # Acid mine 5, Cryo water 8, Desert 1, Soda lake 37
pdf("FigsUpdated2/mp_Order.pdf", width = 5, height = 7)
plot_multipatt_fungal(mp_obj = mp_order, 
                      input = input_fungi_CPM,
                      tax_sum = tax_sum_order,
                      group = "Environment",
                      qcut = 0.05,
                      rcut = 0)
dev.off()

# Family
tax_sum_family <- summarize_taxonomy(input_fungi_CPM, level = 5, report_higher_tax = F, relative = F)
set.seed(425)
mp_family <- multipatt(t(tax_sum_family), 
                      input_fungi_CPM$map_loaded$Environment, 
                      func = "r.g", 
                      control = how(nperm=999))
summary(mp_family) # 92 to 1 group: Acid mine 5, cryo water 19, desert 4, GF 3, soda lake 61, 
pdf("FigsUpdated2/mp_Family.pdf", width = 5, height = 10)
plot_multipatt_fungal(mp_obj = mp_family, 
                      input = input_fungi_CPM,
                      tax_sum = tax_sum_family,
                      group = "Environment",
                      qcut = 0.05,
                      rcut = 0)
dev.off()

#### ___Figure S4 ####
png("FinalFigs/FigureS4.png", width = 5, height = 10, units = "in", res = 300)
plot_multipatt_fungal(mp_obj = mp_family, 
                      input = input_fungi_CPM,
                      tax_sum = tax_sum_family,
                      group = "Environment",
                      qcut = 0.05,
                      rcut = 0)
dev.off()

# Genus
tax_sum_genus <- summarize_taxonomy(input_fungi_CPM, level = 6, report_higher_tax = F, relative = F)
set.seed(425)
mp_genus <- multipatt(t(tax_sum_genus), 
                      input_fungi_CPM$map_loaded$Environment, 
                      func = "r.g", 
                      control = how(nperm=999))
summary(mp_genus) # 148 to 1 group
# Acid mine 8, cryo water 31, desert 10, glacial forefield 2, hot spring 1, soda lake 96
pdf("FigsUpdated2/mp_Genus.pdf", width = 5, height = 12)
plot_multipatt_fungal(mp_obj = mp_genus, 
                      input = input_fungi_CPM,
                      tax_sum = tax_sum_genus,
                      group = "Environment",
                      qcut = 0.05,
                      rcut = 0)
dev.off()



#### ..........................####
#### 4. Functional ####
# Sent Dongying Wu of IMG staff list of taxonoids and list of fungal phyla
# Dongying ran custom python scripts on JGI super computer to pull out KOs of only fungal phyla
# Folder FungalKOs has a file for each metagenome with the KO hits of the fungal phyla scaffolds
# Already deleted 318 blank files (no fungal KOs); 802 had at least 1 KO
# Note that there is bias in eukaryotic gene calling/KO assignment
# We could also get COG or Pfam profiles if we want
# Code below is how DESeq-transformed abundance table was created
# Need to update with the final 9 cryo and also filter unwanted samples
# 2 of the 9 cryos didn't have any fungal KOs
# 3300038653 and 3300039024
# To reload the KO table and metadata, go to "_Start here" subsection

# Run a for loop to read in the file for each metagenome and combine into 1
setwd("data/FungalKOs/")
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
new_kos <- read.table("data/fungi.ko.txt")
ko_table <- rbind(ko_table, new_kos)

# Add new cryo samples (7 of 9)
setwd("data/FungalKOs_Cryo9/")
list.files()
ko_cry <- list()
ko_input_cry <- data.frame(V1 = "NA",
                           V2 = "NA",
                           V3 = "NA")
ko_table_cry <- ko_input_cry
for (i in 1:length(list.files())) {
  ko_cry[[i]] <- read.delim(list.files()[i], header = F)
  ko_table_cry <- ko_table_cry %>%
    rbind(ko_cry[[i]])
}
setwd("~/Documents/GitHub/Extremophilic_Fungi/")

ko_table <- rbind(ko_table, ko_table_cry)


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
# write.csv(ko_table$KO, file = "KOlist_updated2.csv", row.names = F)

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
# Takes a while so don't want to do everything again!
ko_list1 <- read.csv("data/KOlist_wDefinitions.csv")
ko_list2 <- read.csv("data/KOlist_wDefinitions_new.csv")
old_ko_list <- rbind(ko_list1, ko_list2)

new_ko_list <- ko_list %>%
  filter(KO %notin% old_ko_list$KO)

# for (i in 1:nrow(new_ko_list)) {
#  def <- keggFind(database = "ko", query = new_ko_list$KO[i])
#  if (length(def) != 0) {
#    new_ko_list$Definition[i] <- def
#  }
# }
# write.csv(new_ko_list, file = "KOlist_wDefinitions_new2.csv", row.names = F)

ko_list <- rbind(old_ko_list, new_ko_list)
ko_list$KO_def <- paste(ko_list$KO, ko_list$Definition, sep = " ")

# Make community style table and metadata, match IDs
nrow(input_fungi$map_loaded) # 855
nrow(input_fungi_nz$map_loaded) # 732
ko_comm <- ko_table_MGcount %>%
  t() %>%
  as.data.frame() %>%
  filter(rownames(.) %in% input_fungi$map_loaded$taxon_oid) %>%
  arrange(rownames(.))
nrow(ko_comm) # 676. So 56 with fungi but without KO. 

ko_meta <- input_fungi$map_loaded %>%
  filter(taxon_oid %in% rownames(ko_comm)) %>%
  arrange(taxon_oid) %>%
  left_join(., methods, by = "sampleID") %>%
  mutate(rn = sampleID) %>%
  column_to_rownames(var = "rn")

# Check match (should be zero)
sum(rownames(ko_comm) != ko_meta$taxon_oid)
rownames(ko_comm) <- rownames(ko_meta)

# Check environment sample size
table(ko_meta$Environment) # lowest is soda lake with 18
nrow(ko_meta) # 676

# Get richness and Shannon
ko_meta$richness_KO = specnumber(ko_comm)
ko_meta$shannon_KO = diversity(ko_comm, index = "shannon")
range(ko_meta$richness_KO)
range(ko_meta$shannon_KO)

# Save ko_meta
#saveRDS(ko_meta, "data/ko_meta_final.rds")

pdf("FigsUpdated2/KO_Genus_richness.pdf", width = 6, height = 4)
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
# dds <- DESeqDataSetFromMatrix(countData = t(ko_comm) + 1,
#                              colData = ko_meta,
#                              design = ~ 1)
# dds <- estimateSizeFactors(dds)
# dds <- estimateDispersions(dds)
# ko_comm_DESeq <- as.data.frame(t(counts(dds, normalized = T)))
#Save so you don't have to redo the DESeq (takes a while)
#saveRDS(ko_comm_DESeq, "data/ko_comm_DESeq_updated2.rds")



#### _Start here ####
ko_comm_DESeq <- readRDS("data/ko_comm_DESeq_updated2.rds")
ko_meta <- readRDS("data/ko_meta_final.rds")
table(ko_meta$Environment)



#### _KO Alpha ####
hist(ko_meta$richness_KO)
leveneTest(richness_KO ~ Environment, data = ko_meta) # Bad
m5 <- aov(richness_KO ~ Environment, data = ko_meta)
summary(m5)
shapiro.test(m5$residuals) # Bad
plot(m5$fitted.values, m5$residuals)
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

hist(ko_meta$shannon_KO)
leveneTest(shannon_KO ~ Environment, data = ko_meta) # Bad
m6 <- aov(shannon_KO ~ Environment, data = ko_meta)
summary(m6)
shapiro.test(m6$residuals) # Bad
plot(m6$fitted.values, m6$residuals)
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
pdf("FigsUpdated2/KO_AlphaDiversity.pdf", width = 6, height = 3)
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

min(ko_meta$richness_KO)
max(ko_meta$richness_KO)
mean(ko_meta$richness_KO)
se(ko_meta$richness_KO)
ncol(ko_comm)
ncol(ko_comm_DESeq)



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
sum(ko_richness$richness_KO == 1) # 45 with just 1 KO
sum(ko_richness$richness_KO == 2) # 38 with just 2 KOs
sum(ko_richness$richness_KO == 3) # 39 with just 3 KOs

# Check relationships with genome size and fungal count and genus richness
plot(ko_meta$GenomeSize, ko_meta$richness_KO)
summary(lm(richness_KO ~ GenomeSize, data = ko_meta))
cor.test(ko_meta$GenomeSize, ko_meta$richness_KO, method = "pearson")

plot(ko_meta$fung_count, ko_meta$richness_KO)
pdf("FigsUpdated2/KO_richness_total.pdf", width = 7, height = 5)
ggplot(ko_meta, aes(fung_count, richness_KO)) +
  geom_hline(yintercept = 750, linetype = "dashed") +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  labs(x = "Number of fungal protein coding sequences",
       y = "KO richness") +
  theme_classic()
dev.off()
png("FigsUpdated2/KO_richness_total.png", width = 7, height = 5, units = "in", res = 300)
ggplot(ko_meta, aes(fung_count, richness_KO)) +
  geom_hline(yintercept = 750, linetype = "dashed") +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  labs(x = "Number of fungal protein coding sequences",
       y = "KO richness") +
  theme_classic()
dev.off()
png("FinalFigs/FigureS5.png", width = 7, height = 5, units = "in", res = 300)
ggplot(ko_meta, aes(fung_count, richness_KO)) +
  geom_hline(yintercept = 750, linetype = "dashed") +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  labs(x = "Number of fungal protein coding sequences",
       y = "KO richness") +
  theme_classic()
dev.off()
summary(lm(richness_KO ~ fung_count, data = ko_meta))
summary(lm(richness_KO ~ poly(fung_count, 2, raw = TRUE), data = ko_meta))
cor.test(ko_meta$fung_count, ko_meta$richness_KO, method = "pearson")

plot(ko_meta$rich, ko_meta$richness_KO)
# Very weird! does this have to do with euk gene calling issues? 
# Does the 100-200 genus assigned scaffold range cause issues?

# Need to find good cutoff with some KOs and still high sample size
sum(ko_richness$richness_KO > 10) # 443 with > 10

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
sum(ko_richness$richness_KO > 100) # 157 with > 100
ko_meta_filt <- ko_meta %>%
  filter(richness_KO > 100)
ko_comm_DESeq_filt <- ko_comm_DESeq %>%
  filter(rownames(.) %in% rownames(ko_meta_filt))
sum(rownames(ko_comm_DESeq_filt) != rownames(ko_meta_filt))
table(ko_meta_filt$Environment)
sort(colSums(ko_comm_DESeq_filt))
sort(rowSums(ko_comm_DESeq_filt))
sort(ko_meta_filt$richness_KO)

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
  geom_point(aes(size = richness_KO), alpha = 0.5) +
  scale_size(range = c(1, 10)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Environment",
       size = "Number of KOs",
       title = "KO Bray-Curtis") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        plot.title = element_text(vjust = 0))
g5_ko
ko_pcoa_leg <- get_legend(g5_ko)
g5_ko <- g5_ko +
  theme(legend.position = "none")

jac_ko <- vegdist(ko_comm_DESeq_filt, method = "jaccard")
pcoa1_ko <- cmdscale(jac_ko, k = nrow(ko_meta_filt) - 1, eig = T)
pcoa1A1 <- round((eigenvals(pcoa1_ko)/sum(eigenvals(pcoa1_ko)))[1]*100, digits = 1)
pcoa1A2 <- round((eigenvals(pcoa1_ko)/sum(eigenvals(pcoa1_ko)))[2]*100, digits = 1)
ko_meta_filt$Axis01j <- scores(pcoa1_ko)[,1]
ko_meta_filt$Axis02j <- scores(pcoa1_ko)[,2]
micro.hullsj <- ddply(ko_meta_filt, c("Environment"), find_hullj)
g6_ko <- ggplot(ko_meta_filt, aes(Axis01j, -Axis02j, colour = Environment)) +
  geom_polygon(data = micro.hullsj, aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F) +
  geom_point(aes(size = richness_KO), alpha = 0.5) +
  scale_size(range = c(1, 10)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Environment",
       size = "Number of KOs",
       title = "KO Jaccard") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        plot.title = element_text(vjust = 0))
g6_ko

pdf("FigsUpdated2/KO_PCoA_min100KOs.pdf", width = 8.5, height = 5)
plot_grid(g5_ko, g6_ko, ko_pcoa_leg, ncol = 3, rel_widths = c(2.5,2.5,1.1))
dev.off()

#### __750 ####
# Now increase the cutoff even more - 750!
sum(ko_richness$richness_KO > 750) # 44 with > 750
ko_meta_filt <- ko_meta %>%
  filter(richness_KO > 750)
ko_comm_DESeq_filt <- ko_comm_DESeq %>%
  filter(rownames(.) %in% rownames(ko_meta_filt))
sum(rownames(ko_comm_DESeq_filt) != rownames(ko_meta_filt))
table(ko_meta_filt$Environment) # No hot spring, oh well
sort(colSums(ko_comm_DESeq_filt))
sort(rowSums(ko_comm_DESeq_filt))
sort(ko_meta_filt$richness_KO)

bc_ko <- vegdist(ko_comm_DESeq_filt, method = "bray")
pcoa_ko <- cmdscale(bc_ko, k = nrow(ko_meta_filt) - 1, eig = T)
pcoaA1 <- round((eigenvals(pcoa_ko)/sum(eigenvals(pcoa_ko)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa_ko)/sum(eigenvals(pcoa_ko)))[2]*100, digits = 1)
ko_meta_filt$Axis01 <- scores(pcoa_ko)[,1]
ko_meta_filt$Axis02 <- scores(pcoa_ko)[,2]
range(ko_meta_filt$Axis01)
range(ko_meta_filt$Axis02)
micro.hulls <- ddply(ko_meta_filt, c("Environment"), find_hull)
g5_ko <- ggplot(ko_meta_filt, aes(Axis01, Axis02, colour = Environment)) +
  geom_polygon(data = micro.hulls, aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F) +
  geom_point(aes(size = richness_KO), shape = 16, alpha = 0.75) +
  scale_size(range = c(1, 10)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Environment",
       size = "Number of KOs",
       title = "KO Bray-Curtis") +
  scale_fill_manual(values = color_mapping) +
  scale_colour_manual(values = color_mapping) +
  xlim(-0.2, 0.29) +
  ylim(-0.3, 0.12) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = 0),
        panel.grid = element_blank())
g5_ko
ko_pcoa_leg <- get_legend(g5_ko)
g5_ko <- g5_ko +
  theme(legend.position = "none")

jac_ko <- vegdist(ko_comm_DESeq_filt, method = "jaccard")
pcoa1_ko <- cmdscale(jac_ko, k = nrow(ko_meta_filt) - 1, eig = T)
pcoa1A1 <- round((eigenvals(pcoa1_ko)/sum(eigenvals(pcoa1_ko)))[1]*100, digits = 1)
pcoa1A2 <- round((eigenvals(pcoa1_ko)/sum(eigenvals(pcoa1_ko)))[2]*100, digits = 1)
ko_meta_filt$Axis01j <- scores(pcoa1_ko)[,1]
ko_meta_filt$Axis02j <- scores(pcoa1_ko)[,2]
range(ko_meta_filt$Axis01j)
range(ko_meta_filt$Axis02j)
micro.hullsj <- ddply(ko_meta_filt, c("Environment"), find_hullj)
g6_ko <- ggplot(ko_meta_filt, aes(Axis01j, Axis02j, colour = Environment)) +
  geom_polygon(data = micro.hullsj, aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F) +
  geom_point(aes(size = richness_KO), shape = 16, alpha = 0.75) +
  scale_size(range = c(1, 10)) +
  labs(x = paste("PC1: ", pcoa1A1, "%", sep = ""), 
       y = paste("PC2: ", pcoa1A2, "%", sep = ""),
       colour = "Environment",
       size = "Number of KOs",
       title = "KO Jaccard") +
  scale_fill_manual(values = color_mapping) +
  scale_colour_manual(values = color_mapping) +
  xlim(-0.2, 0.29) +
  ylim(-0.3, 0.12) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = 0),
        panel.grid = element_blank())
g6_ko

pdf("FigsUpdated2/KO_PCoA_min750KOs.pdf", width = 8.5, height = 5)
plot_grid(g5_ko, g6_ko, ko_pcoa_leg, ncol = 3, rel_widths = c(2.5,2.5,1.1),
          labels = c("A", "B"))
dev.off()

#### ___Figure 5 ####
png("FinalFigs/Figure5.png", width = 8.5, height = 5, units = "in", res = 300)
plot_grid(g5_ko, g6_ko, ko_pcoa_leg, ncol = 3, rel_widths = c(2.5,2.5,1.1),
          labels = c("A", "B"))
dev.off()

# Also do with the same subset of the data as taxonomy (19 samples from each env)
# Note numbers are slightly different here because of 0 fungal KOs in some samples
table(ko_meta$Environment)
ko_meta_subset19 <- ko_meta %>%
  filter(sampleID %in% subset19$sampleID)
ko_comm_DESeq_subset19 <- ko_comm_DESeq %>%
  filter(rownames(.) %in% rownames(ko_meta_subset19))
sum(rownames(ko_comm_DESeq_subset19) != rownames(ko_meta_subset19))
table(ko_meta_subset19$Environment)
sort(colSums(ko_comm_DESeq_subset19))
sort(rowSums(ko_comm_DESeq_subset19))
sort(ko_meta_subset19$richness_KO)

bc_ko <- vegdist(ko_comm_DESeq_subset19, method = "bray")
pcoa_ko <- cmdscale(bc_ko, k = nrow(ko_meta_subset19) - 1, eig = T)
pcoaA1 <- round((eigenvals(pcoa_ko)/sum(eigenvals(pcoa_ko)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa_ko)/sum(eigenvals(pcoa_ko)))[2]*100, digits = 1)
ko_meta_subset19$Axis01 <- scores(pcoa_ko)[,1]
ko_meta_subset19$Axis02 <- scores(pcoa_ko)[,2]
micro.hulls <- ddply(ko_meta_subset19, c("Environment"), find_hull)
g7_ko <- ggplot(ko_meta_subset19, aes(Axis01, Axis02, colour = Environment)) +
  geom_polygon(data = micro.hulls, aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F) +
  geom_point(aes(size = richness_KO), alpha = 0.5) +
  scale_size(range = c(1, 10)) +
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

jac_ko <- vegdist(ko_comm_DESeq_subset19, method = "jaccard")
pcoa1_ko <- cmdscale(jac_ko, k = nrow(ko_meta_subset19) - 1, eig = T)
pcoa1A1 <- round((eigenvals(pcoa1_ko)/sum(eigenvals(pcoa1_ko)))[1]*100, digits = 1)
pcoa1A2 <- round((eigenvals(pcoa1_ko)/sum(eigenvals(pcoa1_ko)))[2]*100, digits = 1)
ko_meta_subset19$Axis01j <- scores(pcoa1_ko)[,1]
ko_meta_subset19$Axis02j <- scores(pcoa1_ko)[,2]
micro.hullsj <- ddply(ko_meta_subset19, c("Environment"), find_hullj)
g8_ko <- ggplot(ko_meta_subset19, aes(Axis01j, Axis02j, colour = Environment)) +
  geom_polygon(data = micro.hullsj, aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F) +
  geom_point(aes(size = richness_KO), alpha = 0.5) +
  scale_size(range = c(1, 10)) +
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

pdf("FigsUpdated2/KO_PCoA_n19.pdf", width = 8.5, height = 3.5)
plot_grid(g7_ko, g8_ko, ko_pcoa_leg, ncol = 3, rel_widths = c(2.5,2.5,1))
dev.off()



#### __Stats ####
# Rerun whichever richness filter you want, then run this
# 750
bc_ko <- vegdist(ko_comm_DESeq_filt, method = "bray")
jac_ko <- vegdist(ko_comm_DESeq_filt, method = "jaccard")
set.seed(308)
adonis2(bc_ko ~ Environment, data = ko_meta_filt) # R2 = 0.32, p = 0.001
anova(betadisper(bc_ko, ko_meta_filt$Environment)) # Dispersion homogeneous

set.seed(308)
adonis2(jac_ko ~ ko_meta_filt$Environment) # R2 = 0.30, p = 0.001
anova(betadisper(jac_ko, ko_meta_filt$Environment)) # Dispersion not homogeneous

# Subset 19 (rerun that section first)
bc_ko <- vegdist(ko_comm_DESeq_subset19, method = "bray")
jac_ko <- vegdist(ko_comm_DESeq_subset19, method = "jaccard")
set.seed(308)
adonis2(bc_ko ~ Environment, data = ko_meta_subset19) # R2 = 0.10, p = 0.015
anova(betadisper(bc_ko, ko_meta_subset19$Environment)) # Dispersion not homogeneous

set.seed(308)
adonis2(jac_ko ~ ko_meta_subset19$Environment) # R2 = 0.10, p = 0.000
anova(betadisper(jac_ko, ko_meta_subset19$Environment)) # Dispersion not homogeneous

set.seed(308)
adonis2(bc_ko ~ ko_meta_subset19$Assembler) # R2 = 0.16, p = 0.033
anova(betadisper(bc_ko, ko_meta_subset19$Assembler)) # Dispersion not homogeneous

ko_meta_subset19 <- ko_meta_subset19 %>%
  separate(`Add.Date`, into = c("Day", "Month", "Year"), sep = "/", remove = F)
set.seed(308)
adonis2(bc_ko ~ as.factor(ko_meta_subset19$Year)) # R2 = 0.14, p = 0.02
anova(betadisper(bc_ko, ko_meta_subset19$Year)) # Dispersion not homogeneous



#### __Drivers/Indicators ####
# KO MULTIPATT (list KOs associated with each group)
set.seed(425)
mp <- multipatt(ko_comm_DESeq, 
                ko_meta$Environment, 
                func = "IndVal.g", 
                control = how(nperm=999))
summary(mp) # None!!



#### _KO Top ####
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
                                     "Acid mine drainage" = "Acid mine drainage (n = 23)",
                                     "Cryosphere - soil" = "Cryosphere - soil (n = 46)",
                                     "Cryosphere - water" = "Cryosphere - water (n = 38)",
                                     "Desert" = "Desert (n = 80)",
                                     "Glacial forefield" = "Glacial forefield (n = 49)",
                                     "Hot spring" = "Hot spring (n = 109)",
                                     "Hydrothermal vent" = "Hydrothermal vent (n = 174)",
                                     "Hypersaline" = "Hypersaline (n = 139)",
                                     "Soda lake" = "Soda lake (n = 18)"))

gene_plot_summary <- gene_plot_long %>%
  group_by(Environment, Gene) %>%
  summarise(mean = mean(Abundance),
            se = std.error(Abundance))

pdf("FigsUpdated2/KO_Barplot_Top.pdf", width = 8, height = 5)
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
ann_colors <- list(Environment = c("Acid mine drainage" = hue_pal()(9)[1],
                                   "Cryosphere - soil" = hue_pal()(9)[2],
                                   "Cryosphere - water" = hue_pal()(9)[3],
                                   "Desert" = hue_pal()(9)[4],
                                   "Glacial forefield" = hue_pal()(9)[5],
                                   "Hot spring" = hue_pal()(9)[6],
                                   "Hydrothermal vent" = hue_pal()(9)[7],
                                   "Hypersaline" = hue_pal()(9)[8],
                                   "Soda lake" = hue_pal()(9)[9]))
phm1 <- pheatmap(gene_hm,
         scale = "row",
         show_colnames = F,
         cluster_rows = F,
         annotation_col = ann_cols,
         annotation_colors = ann_colors)
save_pheatmap_pdf(phm1, "FigsUpdated2/KO_heatmap_Top.pdf")



#### _KO Stress ####
# Got list of genes from Quandt Lab
#stress_genes <- read_xlsx("Stress_genes_cleaned.xlsx")

# Add definitions in standardized fashion
#for (i in 1:nrow(stress_genes)) {
#  def <- keggFind(database = "ko", query = stress_genes$KO[i])
#  if (length(def) != 0) {
#    stress_genes$Definition[i] <- def
#  }
#}

#stress_genes <- stress_genes %>%
#  dplyr::select(-Symbol, -Name, -Protein) %>%
#  separate(Definition, into = c("Symbol", "Name"), sep = ";", remove = F) %>%
#  dplyr::select(KO, Definition, Symbol, Name, Stress, Notes)
#write.csv(stress_genes, file = "stress_genes_wDef.csv", row.names = F)
stress_genes <- read.csv("data/stress_genes_wDef.csv") %>%
  filter(KO %in% colnames(ko_comm_DESeq)) %>%
  filter(!duplicated(KO))
# Note 158 of 208 KOs are in the dataset
# Note 2 were duplicates. So there are 157 KOs
not_present <- read.csv("data/stress_genes_wDef.csv") %>%
  filter(KO %notin% colnames(ko_comm_DESeq)) %>%
  filter(!duplicated(KO))
#write.csv(not_present, file = "stress_genes_not_present.csv", row.names = F)

# Updated file annotated by Lara
# From 157, filtered to 56
stress_genes <- read.csv("data/stress_genes_sorted_5.19.23_LV.csv") %>%
  filter(Notes == "ok")

# Data frame
gene_plot <- ko_comm_DESeq %>%
  dplyr::select(stress_genes$KO) %>%
  mutate("Environment" = ko_meta$Environment)



#### __ (I) Stats ####
# Run a loop 
kruskal_results_genes <- as.data.frame(matrix(data = NA, ncol(gene_plot)-1, 3)) %>%
  set_names(c("Gene", "X2", "P"))
for (i in 1:(ncol(gene_plot)-1)) {
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
  pivot_longer(names(gene_plot)[1:nrow(stress_genes)], 
               names_to = "Gene", values_to = "Abundance") %>%
  mutate(Gene = as.factor(Gene)) %>%
  droplevels() %>%
  mutate(Environment = dplyr::recode_factor(Environment,
                                            "Acid mine drainage" = "Acid mine drainage (n = 23)",
                                            "Cryosphere - soil" = "Cryosphere - soil (n = 46)",
                                            "Cryosphere - water" = "Cryosphere - water (n = 38)",
                                            "Desert" = "Desert (n = 80)",
                                            "Glacial forefield" = "Glacial forefield (n = 49)",
                                            "Hot spring" = "Hot spring (n = 109)",
                                            "Hydrothermal vent" = "Hydrothermal vent (n = 174)",
                                            "Hypersaline" = "Hypersaline (n = 139)",
                                            "Soda lake" = "Soda lake (n = 18)"))

gene_plot_summary <- gene_plot_long %>%
  group_by(Environment, Gene) %>%
  summarise(mean = mean(Abundance),
            se = std.error(Abundance))

# Too many to do barplot
ggplot(gene_plot_summary, aes(Environment, mean, fill = Gene, group = Gene)) +
  geom_bar(stat = "identity", position = position_dodge(0.75)) +
  geom_linerange(aes(x = Environment, ymin = mean - se, ymax = mean + se, 
                     group = Gene),
                 position = position_dodge(0.75)) +
  labs(x = NULL, 
       y = "Abundance (DESeq2 normalized)",
       fill = "KO") +
  #scale_fill_manual(values = brewer.pal(10, "Paired"),
  #                  labels = ko_list$KO_def[1:10]) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"),
        legend.position = "right",
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.7), "cm"))

# Heatmap
# Use pheatmap (pretty heatmap) package
# Sort by environment
ko_meta_sorted <- ko_meta %>%
  arrange(Environment)
reorder_idx <- match(rownames(ko_meta_sorted), rownames(ko_comm_DESeq))
ko_comm_DESeq_sorted <- ko_comm_DESeq[reorder_idx,]
sum(rownames(ko_meta_sorted) != rownames(ko_comm_DESeq_sorted))

# Sort by stress
stress_genes_sorted <- stress_genes %>%
  arrange(Stress)
#write.csv(stress_genes_sorted, "stress_genes_sorted_5.19.23.csv")
stress_genes_sorted$KOsymb <- paste(stress_genes_sorted$KO, stress_genes_sorted$Symbol, sep = "; ")

gene_hm <- ko_comm_DESeq_sorted %>%
  dplyr::select(stress_genes_sorted$KO) %>%
  mutate("sampleID" = ko_meta_sorted$sampleID)
rownames(gene_hm) <- gene_hm$sampleID
gene_hm <- gene_hm %>%
  dplyr::select(-sampleID) %>%
  t() %>%
  as.matrix()

ann_cols <- data.frame(row.names = colnames(gene_hm), 
                       Environment = ko_meta_sorted$Environment)
ann_rows <- data.frame(row.names = rownames(gene_hm), 
                       Stress = stress_genes_sorted$Stress)
ann_colors <- list(Environment = c("Acid mine drainage" = hue_pal()(9)[1],
                                   "Cryosphere - soil" = hue_pal()(9)[2],
                                   "Cryosphere - water" = hue_pal()(9)[3],
                                   "Desert" = hue_pal()(9)[4],
                                   "Glacial forefield" = hue_pal()(9)[5],
                                   "Hot spring" = hue_pal()(9)[6],
                                   "Hydrothermal vent" = hue_pal()(9)[7],
                                   "Hypersaline" = hue_pal()(9)[8],
                                   "Soda lake" = hue_pal()(9)[9]),
                   Stress = c("Alkaline pH" = viridis_pal()(17)[1],
                              "Cellular" = viridis_pal()(17)[2],
                              "Chlorine" = viridis_pal()(17)[3],
                              "Cold" = viridis_pal()(17)[4],
                              "Environmental" = viridis_pal()(17)[5],
                              "Environmental-Oxidative" = viridis_pal()(17)[6],
                              "General" = viridis_pal()(17)[7],
                              "Granule" = viridis_pal()(17)[8],
                              "Heat" = viridis_pal()(17)[9],
                              "Heat-Osmotic" = viridis_pal()(17)[10],
                              "High pH" = viridis_pal()(17)[11],
                              "Metal resistance-Zn" = viridis_pal()(17)[12],
                              "Osmotic" = viridis_pal()(17)[13],
                              "Osmotic, Low pH" = viridis_pal()(17)[14],
                              "Oxidative" = viridis_pal()(17)[15],
                              "pH" = viridis_pal()(17)[16],
                              "Starvation" = viridis_pal()(17)[17]))
phm1 <- pheatmap(gene_hm,
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                 border_color = NA,
                 scale = "row",
                 show_colnames = F,
                 angle_col = 315,
                 cluster_rows = F,
                 cluster_cols = F,
                 fontsize_row = 4,
                 annotation_col = ann_cols,
                 annotation_row = ann_rows,
                 annotation_colors = ann_colors)

# Still no good. Summarize by environment so there are only 9 columns
gene_plot_sorted <- ko_comm_DESeq_sorted %>%
  dplyr::select(stress_genes_sorted$KO) %>%
  mutate("Environment" = ko_meta_sorted$Environment)

gene_plot_long_sorted <- gene_plot_sorted %>%
  pivot_longer(names(gene_plot)[1:nrow(stress_genes)], 
               names_to = "Gene", values_to = "Abundance") %>%
  mutate(Environment = dplyr::recode_factor(Environment,
                                            "Acid mine drainage" = "Acid mine drainage (n = 23)",
                                            "Cryosphere - soil" = "Cryosphere - soil (n = 46)",
                                            "Cryosphere - water" = "Cryosphere - water (n = 38)",
                                            "Desert" = "Desert (n = 80)",
                                            "Glacial forefield" = "Glacial forefield (n = 49)",
                                            "Hot spring" = "Hot spring (n = 109)",
                                            "Hydrothermal vent" = "Hydrothermal vent (n = 174)",
                                            "Hypersaline" = "Hypersaline (n = 139)",
                                            "Soda lake" = "Soda lake (n = 18)"))

gene_plot_summary_sorted <- gene_plot_long_sorted %>%
  group_by(Environment, Gene) %>%
  summarise(mean = mean(Abundance),
            se = std.error(Abundance))

gene_hm_summary <- gene_plot_summary_sorted %>%
  dplyr::select(-se) %>%
  pivot_wider(names_from = Environment, values_from = mean) %>%
  column_to_rownames(var = "Gene") %>%
  t() %>%
  as.data.frame() %>%
  dplyr::select(stress_genes_sorted$KO) %>%
  t() %>%
  as.matrix()
rownames(gene_hm_summary) <- stress_genes_sorted$KOsymb

ann_rows <- data.frame(row.names = rownames(gene_hm_summary), 
                       Stress = stress_genes_sorted$Stress)
# Shortened list
ann_colors <- list(Stress = c("Cold" = viridis_pal()(5)[1],
                              "General" = viridis_pal()(5)[2],
                              "Heat" = viridis_pal()(5)[3],
                              "Osmotic" = viridis_pal()(5)[4],
                              "Oxidative" = viridis_pal()(5)[5]))

phm1 <- pheatmap(gene_hm_summary,
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                 legend = T,
                 legend_breaks = c(-2, -1, 0, 1, 2, 2.66),
                 legend_labels = c("-2", "-1", "0", "1", "2", "Abund."),
                 border_color = NA,
                 scale = "row",
                 show_colnames = T,
                 angle_col = 315,
                 cluster_rows = F,
                 cluster_cols = T,
                 method = "ward.D2",
                 fontsize_row = 8,
                 gaps_row = c(12, 25, 29, 38),
                 annotation_row = ann_rows,
                 annotation_colors = ann_colors)
save_pheatmap_pdf <- function(x, filename, width = 7, height = 12) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(phm1, "FigsUpdated2/KO_heatmap_Stress_filtered.pdf")

# Figure 6
pheatmap(gene_hm_summary,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
         legend = T,
         legend_breaks = c(-2, -1, 0, 1, 2, 2.66),
         legend_labels = c("-2", "-1", "0", "1", "2", "Abund."),
         border_color = NA,
         scale = "row",
         show_colnames = T,
         angle_col = 315,
         cluster_rows = F,
         cluster_cols = T,
         clustering_method = "ward.D2",
         fontsize_row = 10,
         gaps_row = c(12, 25, 29, 38),
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         file = "FinalFigs/Figure6.png",
         height = 12,
         width = 7)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# Explore clustering by rows too
phm2 <- pheatmap(gene_hm_summary,
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                 legend = T,
                 legend_breaks = c(-2, -1, 0, 1, 2, 2.65),
                 legend_labels = c("-2", "-1", "0", "1", "2", "Abund."),
                 border_color = NA,
                 scale = "row",
                 show_colnames = T,
                 angle_col = 315,
                 cluster_rows = T,
                 cluster_cols = T,
                 fontsize_row = 4,
                 annotation_row = ann_rows,
                 annotation_colors = ann_colors)




#### ..........................####
#### 5. Other ####
#### _Map ####
# Sample map with ggplot, color by environment
# Need to adjust color scheme here and throughout
# Map will now be Figure 1a
world <- map_data("world")
nrow(input$map_loaded)
input$map_loaded$Latitude <- as.numeric(input$map_loaded$Latitude)
input$map_loaded$Longitude <- as.numeric(input$map_loaded$Longitude)

no_coord <- input$map_loaded %>%
  filter(is.na(Latitude) == T | is.na(Longitude) == T)
table(no_coord$Environment)

# 10 missing
pdf("FigsUpdated2/SampleMap.pdf", width = 8, height = 5)
ggplot() +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "white", fill = "gray90", linewidth = 0.1, alpha = 0.75) +
  geom_point(data = input$map_loaded, 
             aes(x = Longitude, y = Latitude,
                 colour = Environment, shape = Environment),
             size = 2, alpha = 1) +
  scale_shape_manual(values = c(16, 17, 16, 17, 16, 17, 16, 17, 16)) +
  scale_color_manual(values = color_mapping) +
  theme_void() +
  labs(x = NULL,
       y = NULL) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()

fig1a <- ggplot() +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "white", fill = "gray90", linewidth = 0.1, alpha = 0.75) +
  geom_point(data = input$map_loaded, 
             aes(x = Longitude, y = Latitude,
                 colour = Environment),
             size = 2, alpha = 1) +
  #scale_shape_manual(values = c(16, 17, 16, 17, 16, 17, 16, 17, 16)) +
  scale_color_manual(values = color_mapping) +
  theme_void() +
  labs(x = NULL,
       y = NULL) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.margin = margin(c(0,5,5,5)),
        legend.box.margin = margin(c(0,0,0,0)),
        legend.background = element_rect(linewidth = 0.1))
fig1a



#### _Geog. Distance ####
# Calculate geodesic distance among samples
# Test distance decay relationships in community comp.
# Mantel tests between GD and BC
library(geodist)

# Need to filter datasets to samples with coordinates
missing <- input_fungi$map_loaded %>%
  filter(is.na(Longitude) == TRUE)

dat_fun <- filter_data(input_fungi_CPM_nz,
                       filter_cat = "taxon_oid",
                       filter_vals = missing$taxon_oid) # 730

dat_arc <- filter_data(input_arc_CPM_nz,
                       filter_cat = "taxon_oid",
                       filter_vals = missing$taxon_oid) # 804

dat_bac <- filter_data(input_bac_CPM,
                       filter_cat = "taxon_oid",
                       filter_vals = missing$taxon_oid) # 845

coords <- dat_fun$map_loaded %>%
  dplyr::select(Latitude, Longitude)
gd <- geodist(coords, measure = "geodesic")
gd <- dist(gd)
bc <- calc_dm(dat_fun$data_loaded)
set.seed(1130)
mantel(gd, bc, method = "pearson", permutations = 2000, na.rm = T)
# Significant, p = 0.0005, r = 0.12

# Compare to arc and bac
coords_a <- dat_arc$map_loaded %>%
  dplyr::select(Latitude, Longitude)
gd_a <- geodist(coords_a, measure = "geodesic")
gd_a <- dist(gd_a)
bc_a <- calc_dm(dat_arc$data_loaded)
set.seed(1130)
mantel(gd_a, bc_a, method = "pearson", permutations = 2000, na.rm = T)
# Significant, p - 0.0005, r = 0.08

coords_b <- dat_bac$map_loaded %>%
  dplyr::select(Latitude, Longitude)
gd_b <- geodist(coords_b, measure = "geodesic")
gd_b <- dist(gd_b)
bc_b <- calc_dm(dat_bac$data_loaded)
set.seed(1130)
mantel(gd_b, bc_b, method = "pearson", permutations = 2000, na.rm = T)
# Significant, p - 0.04, r = 0.02

# So fungi have a stronger relationship than arc and bac
# Plot correlogs and distances
set.seed(1130)
mc <- mantel.correlog(D.eco = bc, D.geo = gd, r.type = "pearson", nperm = 1000)
plot(mc)
qplot(gd, bc) +
  geom_smooth() +
  theme_bw()

set.seed(1130)
mc_a <- mantel.correlog(D.eco = bc_a, D.geo = gd_a, r.type = "pearson", nperm = 1000)
plot(mc_a)
qplot(gd_a, bc_a) +
  geom_smooth() +
  theme_bw()

set.seed(1130)
mc_b <- mantel.correlog(D.eco = bc_b, D.geo = gd_b, r.type = "pearson", nperm = 1000)
plot(mc_b)
qplot(gd_b, bc_b) +
  geom_smooth() +
  theme_bw()



#### _KO750 Taxa ####
# Get the same 44 samples from the KO750 comparison and make PCoA from fungal_CPM_nz
input_fungi_CPM_nz_44 <- filter_data(input_fungi_CPM_nz,
                                     filter_cat = "taxon_oid",
                                     keep_vals = ko_meta_filt$taxon_oid)
# BC
bc <- calc_dm(input_fungi_CPM_nz_44$data_loaded)
set.seed(1150)
adonis2(bc ~ Environment, data = input_fungi_CPM_nz_44$map_loaded) # R2 = 0.49, p = 0.001
anova(betadisper(bc, input_fungi_CPM_nz_44$map_loaded$Environment)) # Dispersion homogeneous
pcoa <- cmdscale(bc, k = nrow(input_fungi_CPM_nz_44$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
input_fungi_CPM_nz_44$map_loaded$Axis01 <- scores(pcoa)[,1]
input_fungi_CPM_nz_44$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(input_fungi_CPM_nz_44$map_loaded, c("Environment"), find_hull)
bc44 <- ggplot(input_fungi_CPM_nz_44$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.75, shape = 16, aes(colour = Environment), show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       title = "Bray-Curtis",
       colour = "Environment") +
  scale_fill_manual(values = color_mapping) +
  scale_colour_manual(values = color_mapping) +
  xlim(-0.48, 0.35) +
  ylim(-0.31, 0.34) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = 0),
        panel.grid = element_blank())
leg44 <- get_legend(bc44)
bc44 <- bc44 + theme(legend.position = "none")
range(input_fungi_CPM_nz_44$map_loaded$Axis01) # -0.41 to 0.29
range(input_fungi_CPM_nz_44$map_loaded$Axis02) # -0.31 to 0.34

jac <- calc_dm(input_fungi_CPM_nz_44$data_loaded, method = "jaccard")
set.seed(1150)
adonis2(jac ~ Environment, data = input_fungi_CPM_nz_44$map_loaded) # R2 = 0.54, p = 0.001
anova(betadisper(jac, input_fungi_CPM_nz_44$map_loaded$Environment)) # Dispersion homogeneous
pcoa <- cmdscale(jac, k = nrow(input_fungi_CPM_nz_44$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
input_fungi_CPM_nz_44$map_loaded$Axis01 <- scores(pcoa)[,1]
input_fungi_CPM_nz_44$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(input_fungi_CPM_nz_44$map_loaded, c("Environment"), find_hull)
jac44 <- ggplot(input_fungi_CPM_nz_44$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.75, shape = 16, aes(colour = Environment), show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       title = "Jaccard") +
  scale_fill_manual(values = color_mapping) +
  scale_colour_manual(values = color_mapping) +
  xlim(-0.48, 0.35) +
  ylim(-0.31, 0.34) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = 0),
        panel.grid = element_blank())
jac44
range(input_fungi_CPM_nz_44$map_loaded$Axis01) # -0.48 to 0.35
range(input_fungi_CPM_nz_44$map_loaded$Axis02) # -0.16 to 0.18

#### ___Figure S6 ####
png("FinalFigs/FigureS6.png", width = 8.5, height = 5, units = "in", res = 300)
plot_grid(bc44, jac44, leg44, ncol = 3, rel_widths = c(2.5,2.5,1.1),
          labels = c("A", "B"))
dev.off()

#### _Table S1 ####
# Add number of scaffolds
scaff_cryo <- read.delim("data/ScaffoldCount_cryo.txt") %>%
  dplyr::select(taxon_oid, Scaffold.Count.....assembled)
scaff <- read.delim("data/ScaffoldCount.txt") %>%
  dplyr::select(taxon_oid, Scaffold.Count.....assembled) %>%
  rbind(., scaff_cryo)
ts1 <- read_excel("data/TableS1_noscaff.xlsx") %>%
  left_join(., scaff, by = "taxon_oid") %>%
  mutate(ScaffoldCount = Scaffold.Count.....assembled) %>%
  dplyr::select(-Scaffold.Count.....assembled, -IMG.Genome.ID, -Day, -Month, -Year, -sampleID,
                -GenomeSize, -Study.Name2, -Location2, -EnvGeo) %>%
  arrange(Environment, Study.Name, Genome.Name...Sample.Name)
#write_xlsx(ts1, "data/TableS1.xlsx", format_headers = F)

names(ts1)
ts1$ScaffoldCount
sum(is.na(ts1$ScaffoldCount))
check <- ts1 %>%
  filter(is.na(ScaffoldCount) == T)
# Don't have 6 from Schmidt Lab. 



#### End Script ####