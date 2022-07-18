# Extremophilic Fungi Metagenome Analysis
# For Quandt Lab Review Paper
# by Cliff Bueno de Mesquita, JGI, July 2022

#### Retrieving Data ####
# Info on assembling the IMG dataset is in the Data Acquisition.text file
# Ecosystems of interest were searched by ecosystem or name
# Domain = *Microbiome
# GOLD Analysis Project Type = Metagenome Analysis
# Preprocessing yielded 1560 samples ("Extreme_Combined_Prefilt.txt")
# Postprocessing:
# -Filter Sequence center != DOE Joint Genome Institute (JGI) and Is Public = Yes
# -Filter Sequence center = DOE Joint Genome Institute (JGI) and JGI Data Utilization Status = Unrestricted
# -Filter out duplicates that were reassembled with SPAdes and contain (SPAdes) at the end of the Genome Name/ Sample Name
# -Filter out engineered ecosystems (e.g. bioreactors)
# -Check that there are no duplicate taxonoids
# -Make comma separated list of taxonoids to search for in Genome Search

#### Setup ####
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
suppressWarnings(suppressMessages(library(FSA))) # For se
suppressWarnings(suppressMessages(library(mctoolsr))) # For taxonomic analysis
suppressWarnings(suppressMessages(library(cowplot))) # For multipanel
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

setwd("~/Documents/GitHub/Extremophilic_Fungi/")
options(max.print = 20000000)
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
find_hullj <- function(df) df[chull(df$Axis01j, df$Axis02j),]
`%notin%` <- Negate(`%in%`)

#### Postprocessing ####
prefilt <- read.delim2("Extreme_Combined_Prefilt.txt", header = T) %>%
  mutate(Use = ifelse(Sequencing.Center != "DOE Joint Genome Institute (JGI)" &
                        Is.Public == "Yes",
                      "Use",
                      ifelse(Sequencing.Center != "DOE Joint Genome Institute (JGI)" &
                               Is.Public == "No",
                             "Don't use",
                             ifelse(Sequencing.Center == "DOE Joint Genome Institute (JGI)" &
                                      JGI.Data.Utilization.Status == "Unrestricted",
                                    "Use",
                                    "Don't use"))))
table(prefilt$Use)
postfilt <- prefilt %>%
  filter(Phylum != "Engineered") %>%
  filter(Use == "Use") %>%
  filter(!grepl("SPAdes", Genome.Name...Sample.Name)) %>%
  mutate(sampleID = paste("X", taxon_oid, sep = "")) %>%
  mutate(Environment = 
           ifelse(Family == "Hydrothermal vents",
                  "Hydrothermal vent",
                  ifelse(Family == "Hot (42-90C)",
                         "Hot spring",
                         ifelse(Family == "Hypersaline" |
                                  Family == "Salt crystallizer ponds" |
                                  Study.Name ==  "Salt pond water, soil and salt crust microbial communities from South San Francisco under conditions of wetland restoration.",
                                "Hypersaline",
                                ifelse(Family == "Alkaline",
                                       "Soda lake",
                                       ifelse(Family == "Ice",
                                              "Cryosphere",
                                              ifelse(Family == "Glacier",
                                                     "Glacial forefield",
                                                     ifelse(Family == "Desert",
                                                            "Desert",
          ifelse(grep("Acid mine drainage", Study.Name, ignore.case = T),
                 "Acid mine drainage",
                 "Heavy metal"))))))))) %>%
  select(sampleID, everything())

# Tried making own heavy metal category, but later decided to just merge it with Acid mine drainage
# for (i in 1:nrow(postfilt)) {
#   if (grepl("Heavy metal|arsenic|copper|mercury", postfilt$Study.Name[i])) {
#     postfilt$Environment[i] <- "Heavy metal"
#   }
# }

table(postfilt$Environment)
checkn <- as.data.frame(table(postfilt$Environment))
sum(checkn$Freq) == nrow(postfilt) # Good
sum(duplicated(postfilt$taxon_oid)) # Good

# Save as metadata file
write.table(postfilt, "metadata.txt", sep = "\t" , row.names = F)

# Get list, note that you can only search for 1000 at a time so split up
# We need to split into 2 anyways for the Statistical Analysis tool
# Note: this was done before removing engineered samples, but doesn't matter because we filter input by the metadata.txt file
taxonoid_list1 <- postfilt %>%
  select(taxon_oid) %>%
  mutate(comma = ",") %>%
  mutate(list = paste(taxon_oid, comma, sep = "")) %>%
  select(list) %>%
  slice_head(n = 1000)
taxonoid_list1$list[nrow(taxonoid_list1)] <- gsub(",", "", taxonoid_list1$list[nrow(taxonoid_list1)])
write.csv(taxonoid_list1, "taxonoid_list1.csv", row.names = F)

taxonoid_list2 <- postfilt %>%
  select(taxon_oid) %>%
  mutate(comma = ",") %>%
  mutate(list = paste(taxon_oid, comma, sep = "")) %>%
  select(list) %>%
  slice_tail(n = nrow(postfilt) - 1000)
taxonoid_list2$list[nrow(taxonoid_list2)] <- gsub(",", "", taxonoid_list2$list[nrow(taxonoid_list2)])
write.csv(taxonoid_list2, "taxonoid_list2.csv", row.names = F)

# Now use Genome Search by taxonoid to get a set of those and download genus taxonomy table
# Send set to Dongying Wu for fungal function

#### ..........................####
#### Taxonomic ####
# Metadata ("mapping file") downloaded from IMG and processed above

# Tax table for mctoolsr (only need to do once, skip to Import)
# t <- read.delim("Extreme_1155/UI_data_output.txt") %>%
#  dplyr::select(-c(3:9)) %>%
#  dplyr::rename(taxonomy = FeatureName) %>%
#  dplyr::rename(ASV_ID = Feature) %>%
#  select(ASV_ID, 3:ncol(.), taxonomy)
# names(t) <- abbreviate(names(t), minlength = 11)
# table.fp <- "~/Documents/GitHub/Extremophilic_Fungi"
# out_fp <- paste0(table.fp, "/genus_table_mctoolsr.txt")
# names(t)[1] = "#ASV_ID"
# write("#Exported for mctoolsr", out_fp)
# suppressWarnings(write.table(t, out_fp, sep = "\t", row.names = FALSE, append = TRUE))

#### _Setup ####
# Import with mctoolsr (matches sampleIDs, 1141 samples)
tax_table_fp <- file.path("genus_table_mctoolsr.txt")
map_fp <- file.path("metadata.txt")
input = load_taxa_table(tax_table_fp, map_fp)

# Update map_loaded - sampleID, GenomeSize, Environment
input$map_loaded <- input$map_loaded %>%
  mutate(sampleID = paste("X", taxon_oid, sep = ""))
input$map_loaded$GenomeSize <- input$map_loaded$Genome.Size.....assembled
input$map_loaded$Environment <- as.factor(input$map_loaded$Environment)

# Check sequencing depth 
sort(colSums(input$data_loaded))
mean(colSums(input$data_loaded)) # 308597.4
se(colSums(input$data_loaded)) # 15961.91
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

ggplot(input$map_loaded, aes(Genome.Size.....assembled, count)) +
  geom_point(size = 1.5, alpha = 0.25) +
  geom_smooth(method = "lm") +
  labs(x = "Assembled genome size", 
       y = "Assigned genus reads") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10))

# Filter out samples with no genus level reads (removes 10 samples, 1131 remaining)
count <- as.data.frame(sort(colSums(input$data_loaded))) %>%
  filter(`sort(colSums(input$data_loaded))` == 0)
input <- filter_data(input,
                     filter_cat = "sampleID",
                     filter_vals = rownames(count))

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
range(tax_sum_Domain[3,]) # 0 to 0.63, so can be high

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

# Too many, filter top ones (> 0.5%), want to see which environments the top samples are from
eukabund <- data.frame("Euks" = input$map_loaded$Euks,
                       "sampleID" = input$map_loaded$sampleID) %>%
  filter(Euks > 0.05) %>%
  column_to_rownames(var = "sampleID")
topeuk <- filter_data(input,
                      filter_cat = "sampleID",
                      keep_vals = rownames(eukabund))
pdf("Figs/EukTopSamples.pdf", width = 7, height = 5)
ggplot(topeuk$map_loaded, aes(reorder(sampleID, Euks, mean), Euks, fill = Environment)) +
  geom_bar(stat = "identity", color = NA) +
  labs(x = NULL, 
       y = "Relative abundance") +
  ggtitle("Samples with Eukaryota > 0.5%") +
  theme_classic() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
dev.off()

#### __Check Fungi ####
# First do same as done above for euk but for extracted fungal phyla
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

# Get sum of fungi
fungi <- as.data.frame(t(colSums(fungal_phyla)))
fungi_t <- as.data.frame(t(fungi))

input$map_loaded$Fungi <- fungi_t$V1
ggplot(input$map_loaded, aes(reorder(sampleID, Fungi, mean), Fungi, fill = Environment)) +
  geom_bar(stat = "identity", color = NA) +
  labs(x = NULL, 
       y = "Relative abundance") +
  theme_classic() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Too many, filter top ones (> 0.1%), want to see which environments the top samples are from
funabund <- data.frame("Fungi" = input$map_loaded$Fungi,
                       "sampleID" = input$map_loaded$sampleID) %>%
  filter(Fungi > 0.01) %>%
  column_to_rownames(var = "sampleID")
topfun <- filter_data(input,
                      filter_cat = "sampleID",
                      keep_vals = rownames(funabund))
pdf("Figs/FungalTopSamples.pdf", width = 7, height = 5)
ggplot(topfun$map_loaded, aes(reorder(sampleID, Fungi, mean), Fungi, fill = Environment)) +
  geom_bar(stat = "identity", color = NA) +
  labs(x = NULL, 
       y = "Relative abundance") +
  ggtitle("Samples with Fungi > 0.1%") +
  theme_classic() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
dev.off()

# Also plot total fungal relative abundance by environment
input_fungi_CPM$map_loaded$totalFun <- colSums(input_fungi_CPM$data_loaded)
View(input_fungi_CPM$map_loaded$totalFun)
leveneTest(input_fungi_CPM$map_loaded$totalFun ~ input_fungi_CPM$map_loaded$Environment)
m2 <- aov(totalFun ~ Environment, data = input_fungi_CPM$map_loaded)
shapiro.test(m2$residuals)
summary(m2)
tuk2 <- emmeans(object = m2, specs = "Environment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "totalFun",
         y = max(input_fungi_CPM$map_loaded$totalFun)+(max(input_fungi_CPM$map_loaded$totalFun)-min(input_fungi_CPM$map_loaded$totalFun))/20)
pdf("Figs/FungalRel.pdf", width = 5, height = 4)
ggplot(input_fungi_CPM$map_loaded, aes(reorder(Environment, totalFun, mean), totalFun)) +
  # geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.5, alpha = 0.2, width = 0.4) +
  geom_text(data = tuk2, aes(Environment, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = NULL, y = "Fungal abundance (CPM)") +
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

nrow(input$taxonomy_loaded) # 4031 total
nrow(input_euk$taxonomy_loaded) # 399 euks
nrow(input_fungi$taxonomy_loaded) # 303 fungi

# Now check reads again
sort(colSums(input_fungi$data_loaded))
# Note lots of samples with 0 or very few fungi
# Purposefully not filtering those out those as 0's are interesting in this analysis
# These are extreme envrionments, some may have few to no fungi
mean(colSums(input_fungi$data_loaded)) # 408
se(colSums(input_fungi$data_loaded)) # 51
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

# Check genome size vs fungi
ggplot(input_fungi$map_loaded, aes(Genome.Size.....assembled, fung_count)) +
  geom_point(size = 1.5, alpha = 0.25) +
  geom_smooth(method = "lm") +
  labs(x = "Assembled genome size", 
       y = "Assigned fungal genus reads") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10))

# Get and plot fraction of samples with 0 versus > 0 fungal reads, and also n
env_prev <- input_fungi$map_loaded %>%
  group_by(Environment) %>%
  summarise(num_present = sum(present),
            num_samples = n(),
            prevalence = round(num_present/num_samples * 100, digits = 2)) %>%
  mutate(num_absent = num_samples - num_present)

# Melt for stacked bar
env_prev_long <- melt(env_prev,
                      id.vars = "Environment",
                      measure.vars = c("num_present", "num_absent"))

pdf("Figs/Prevalence.pdf", width = 5, height = 4)
ggplot(env_prev, aes(reorder(Environment, prevalence, mean), prevalence)) +
  geom_bar(stat = "identity") +
  labs(x = NULL, 
       y = "% prevalence of fungi") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

ggplot(env_prev, aes(reorder(Environment, num_samples, mean), num_samples)) +
  geom_bar(stat = "identity") +
  geom_text(data = env_prev, 
            aes(reorder(Environment, num_samples, mean), num_samples+10, 
                          label = num_samples), inherit.aes = F) +
  labs(x = NULL, 
       y = "Sample size") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

pdf("Figs/SampleSize.pdf", width = 5, height = 4.5)
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
  ggtitle("Total sample size = 1131\nSamples with Fungi = 925") +
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
# i is samples 1 to 1145
# j is taxa 1 to 303
for (i in 1:ncol(input_fungi$data_loaded)) {
  for (j in 1:nrow(input_fungi$data_loaded)) {
    input_fungi_CPM$data_loaded[j, i] <- (input_fungi$data_loaded[j, i]*1000000)/input_fungi$map_loaded$GenomeSize[i]
  }
}

# Make stacked bar plots by taxonomic level
# Resort so unassigned and other are on the top
# Use "Paired" palette from RColorBrewer
# Unassigned gets grey75, Other gets grey90
# Show all phyla; for others show top 15
ntax <- 15
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(ntax)
mycolors2 <- colorRampPalette(brewer.pal(8, "Set2"))(ntax)

# Phyla
tax_sum_Phyla <- summarize_taxonomy(input_fungi_CPM, level = 2, report_higher_tax = F, relative = F)
pdf("Figs/CPM_Phyla.pdf", width = 7, height = 5)
plot_taxa_bars(tax_sum_Phyla,
               input_fungi_CPM$map_loaded,
               type_header = "Environment",
               num_taxa = nrow(tax_sum_Phyla),
               data_only = F) +
  scale_fill_brewer(palette = "Paired") +
  theme_classic() +
  labs(x = "Environment",
       y = "Abundance (CPM)",
       fill = "Phylum") +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

# Class
tax_sum_Class <- summarize_taxonomy(input_fungi_CPM, level = 3, report_higher_tax = T, relative = F)
rownames(tax_sum_Class) <- substring(rownames(tax_sum_Class), 12)
plot_taxa_bars(tax_sum_Class,
               input_fungi_CPM$map_loaded,
               type_header = "Environment",
               num_taxa = nrow(tax_sum_Class),
               data_only = F) +
  theme_classic() +
  labs(x = "Environment",
       y = "Abundance (CPM)",
       fill = "Class") +
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

pdf("Figs/CPM_Class.pdf", width = 7, height = 5)
ggplot(barsClass, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Abundance (CPM)", fill = "Class") +
  scale_fill_manual(values = c("grey90", mycolors[15:1])) +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))
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

pdf("Figs/CPM_Order.pdf", width = 7, height = 5)
ggplot(barsOrder, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Abundance (CPM)", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", mycolors[14:1])) +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))
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

pdf("Figs/CPM_Family.pdf", width = 7, height = 5)
ggplot(barsFamily, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Abundance (CPM)", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", mycolors[14:1])) +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))
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

pdf("Figs/CPM_Genus.pdf", width = 7, height = 5)
ggplot(barsGenus, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Abundance (CPM)", fill = "Genus") +
  scale_fill_manual(values = c("grey90", mycolors[15:1])) +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))
dev.off()

taxa_summary_by_sample_type(tax_sum_Genus, input_fungi_CPM$map_loaded, 'Environment', 0.0001, 'KW')

# Also plot total fungal CPM by environment
input_fungi_CPM$map_loaded$totalFun <- colSums(input_fungi_CPM$data_loaded)
View(input_fungi_CPM$map_loaded$totalFun)
leveneTest(input_fungi_CPM$map_loaded$totalFun ~ input_fungi_CPM$map_loaded$Environment)
m2 <- aov(totalFun ~ Environment, data = input_fungi_CPM$map_loaded)
shapiro.test(m2$residuals)
summary(m2)
tuk2 <- emmeans(object = m2, specs = "Environment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "totalFun",
         y = max(input_fungi_CPM$map_loaded$totalFun)+(max(input_fungi_CPM$map_loaded$totalFun)-min(input_fungi_CPM$map_loaded$totalFun))/20)
pdf("Figs/FungalCPM.pdf", width = 5, height = 4)
ggplot(input_fungi_CPM$map_loaded, aes(reorder(Environment, totalFun, mean), totalFun)) +
  # geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.5, alpha = 0.2, width = 0.4) +
  geom_text(data = tuk2, aes(Environment, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = NULL, y = "Fungal abundance (CPM)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        strip.text = element_text(size = 10))
dev.off()

sort(colSums(input_fungi_CPM$data_loaded))

# Could also try plot by location with environment facets



#### _Relative ####
# Calculate relative abundance of fungi, for only samples with fungi (n = 925)
# Make relative abundance stacked bar plots by taxonomic level
# Re-sort so unassigned and other are on the top
# Use "Paired" palette from RColorBrewer
# Unassigned gets grey75, Other gets grey90
# Show all phyla; for others show top 15
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

pdf("Figs/Rel_Phyla.pdf", width = 7, height = 5)
ggplot(barsPhyla, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = brewer.pal(12, "Paired")[7:1]) +
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

pdf("Figs/Rel_Class.pdf", width = 7, height = 5)
ggplot(barsClass, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Relative abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", mycolors[14:1])) +
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

pdf("Figs/Rel_Order.pdf", width = 7, height = 5)
ggplot(barsOrder, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Relative abundance", fill = "Order") +
  scale_fill_manual(values = c("grey90", mycolors[15:1])) +
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

pdf("Figs/Rel_Family.pdf", width = 7, height = 5)
ggplot(barsFamily, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Relative abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", mycolors[14:1])) +
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

pdf("Figs/Rel_Genus.pdf", width = 7, height = 5)
ggplot(barsGenus, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Environment", y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey90", mycolors[15:1])) +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))
dev.off()

taxa_summary_by_sample_type(tax_sum_Genus, input_fungi_nz$map_loaded, 'Environment', 0.0001, 'KW')

# Also plot total fungal relative abundance by environment
input_fungi_nz$map_loaded$totalFun <- colSums(input_fungi_nz$data_loaded)
View(input_fungi_nz$map_loaded$totalFun)
leveneTest(input_fungi_nz$map_loaded$totalFun ~ input_fungi_nz$map_loaded$Environment)
m2 <- aov(totalFun ~ Environment, data = input_fungi_nz$map_loaded)
shapiro.test(m2$residuals)
summary(m2)
tuk2 <- emmeans(object = m2, specs = "Environment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "totalFun",
         y = max(input_fungi_nz$map_loaded$totalFun)+(max(input_fungi_nz$map_loaded$totalFun)-min(input_fungi_nz$map_loaded$totalFun))/20)
pdf("Figs/FungalCPM.pdf", width = 5, height = 4)
ggplot(input_fungi$map_loaded, aes(reorder(Environment, totalFun, mean), totalFun)) +
  # geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.5, alpha = 0.2, width = 0.4) +
  geom_text(data = tuk2, aes(Environment, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = NULL, y = "Fungal Relative abundance") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        strip.text = element_text(size = 10))
dev.off()

sort(colSums(input_fungi$data_loaded))



#### _PCoAs ####
# PCoA and PERMANOVA by Environment
# Also check Habitat, Ecosystem.Category, Ecosystem.Subtype, Ecosystem.Type, Specific.Ecosystem
# Filter out fungal zeroes (filters 206 samples, 925 remaining)
countFun <- as.data.frame(sort(colSums(input_fungi_CPM$data_loaded))) %>%
  filter(`sort(colSums(input_fungi_CPM$data_loaded))` == 0)
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
set.seed(1150)
adonis2(bc ~ Environment, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.074, p = 0.001
adonis2(bc ~ Habitat, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.338, p = 0.001
adonis2(bc ~ Ecosystem.Category, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.077
adonis2(bc ~ Ecosystem.Subtype, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.077
adonis2(bc ~ Ecosystem.Type, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.077
adonis2(bc ~ Specific.Ecosystem, data = input_fungi_CPM_nz$map_loaded) # R2 = 0.077
anova(betadisper(bc, input_fungi_CPM_nz$map_loaded$Environment)) # Dispersion not homogeneous
pcoa <- cmdscale(bc, k = nrow(input_fungi_CPM_nz$map_loaded) - 1, eig = T)
eigenvals(pcoa)/sum(eigenvals(pcoa)) # 12.1, 9.2 % variation explained
input_fungi_CPM_nz$map_loaded$Axis01 <- scores(pcoa)[,1]
input_fungi_CPM_nz$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(input_fungi_CPM_nz$map_loaded, c("Environment"), find_hull)
pdf("Figs/PCoA_Genus.pdf", width = 7, height = 5)
g <- ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = "PC1: 12.1%", 
       y = "PC2: 9.2%") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
g
dev.off()
ggplotly(g)
# There's basically no clustering of fungal composition by environment at genus level, lots of overlap!
# Extract legend to plot separately as its own panel later
g_leg <- get_legend(g)

# Family
bc_Family <- calc_dm(Family_nz)
set.seed(1150)
adonis2(bc_Family ~ Environment, data = input_fungi_CPM_nz$map_loaded) # R2 = , p = 
anova(betadisper(bc_Family, input_fungi_CPM_nz$map_loaded$Environment)) # Dispersion not homogeneous
pcoa_Family <- cmdscale(bc_Family, k = nrow(input_fungi_CPM_nz$map_loaded) - 1, eig = T)
eigenvals(pcoa_Family)/sum(eigenvals(pcoa_Family)) # 12.1, 11.3 % variation explained
input_fungi_CPM_nz$map_loaded$Axis01 <- scores(pcoa_Family)[,1]
input_fungi_CPM_nz$map_loaded$Axis02 <- scores(pcoa_Family)[,2]
micro.hulls <- ddply(input_fungi_CPM_nz$map_loaded, c("Environment"), find_hull)
g1 <- ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F, size = 0.25) +
  geom_point(size = 1, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = "PC1: 12.1%", 
       y = "PC2: 11.3%") +
  ggtitle("Family") +
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
adonis2(bc_Order ~ Environment, data = input_fungi_CPM_nz$map_loaded) # R2 = , p = 
anova(betadisper(bc_Order, input_fungi_CPM_nz$map_loaded$Environment)) # Dispersion not homogeneous
pcoa_Order <- cmdscale(bc_Order, k = nrow(input_fungi_CPM_nz$map_loaded) - 1, eig = T)
eigenvals(pcoa_Order)/sum(eigenvals(pcoa_Order)) # 16.8, 11.7 % variation explained
input_fungi_CPM_nz$map_loaded$Axis01 <- scores(pcoa_Order)[,1]
input_fungi_CPM_nz$map_loaded$Axis02 <- scores(pcoa_Order)[,2]
micro.hulls <- ddply(input_fungi_CPM_nz$map_loaded, c("Environment"), find_hull)
g2 <- ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F, size = 0.25) +
  geom_point(size = 1, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = "PC1: 16.8%", 
       y = "PC2: 11.7%") +
  ggtitle("Order") +
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
adonis2(bc_Class ~ Environment, data = input_fungi_CPM_nz$map_loaded) # R2 = , p = 
anova(betadisper(bc_Class, input_fungi_CPM_nz$map_loaded$Environment)) # Dispersion not homogeneous
pcoa_Class <- cmdscale(bc_Class, k = nrow(input_fungi_CPM_nz$map_loaded) - 1, eig = T)
eigenvals(pcoa_Class)/sum(eigenvals(pcoa_Class)) # 21.4, 11.9 % variation explained
input_fungi_CPM_nz$map_loaded$Axis01 <- scores(pcoa_Class)[,1]
input_fungi_CPM_nz$map_loaded$Axis02 <- scores(pcoa_Class)[,2]
micro.hulls <- ddply(input_fungi_CPM_nz$map_loaded, c("Environment"), find_hull)
g3 <- ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F, size = 0.25) +
  geom_point(size = 1, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = "PC1: 21.4%", 
       y = "PC2: 11.9%") +
  ggtitle("Class") +
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
adonis2(bc_Phylum ~ Environment, data = input_fungi_CPM_nz$map_loaded) # R2 = , p = 
anova(betadisper(bc_Phylum, input_fungi_CPM_nz$map_loaded$Environment)) # Dispersion not homogeneous
pcoa_Phylum <- cmdscale(bc_Phylum, k = nrow(input_fungi_CPM_nz$map_loaded) - 1, eig = T)
eigenvals(pcoa_Phylum)/sum(eigenvals(pcoa_Phylum)) # 41.3, 18.2 % variation explained
input_fungi_CPM_nz$map_loaded$Axis01 <- scores(pcoa_Phylum)[,1]
input_fungi_CPM_nz$map_loaded$Axis02 <- scores(pcoa_Phylum)[,2]
micro.hulls <- ddply(input_fungi_CPM_nz$map_loaded, c("Environment"), find_hull)
g4 <- ggplot(input_fungi_CPM_nz$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F, size = 0.25) +
  geom_point(size = 1, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = "PC1: 41.3%", 
       y = "PC2: 18.2%") +
  ggtitle("Phylum") +
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
  labs(x = "PC1: 12.1%", 
       y = "PC2: 9.2%") +
  ggtitle("Genus") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, vjust = -1))
g

# Multipanel
pdf("Figs/PCoA_AllLevels.pdf", width = 8, height = 5)
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

# Remove zeroes
countArc <- as.data.frame(sort(colSums(input_arc$data_loaded))) %>%
  filter(`sort(colSums(input_arc$data_loaded))` == 0)
input_arc_nz <- filter_data(input_arc,
                            filter_cat = "sampleID",
                            filter_vals = rownames(countArc))

# BC, PERMANOVA, PERMDISP, PCoA
bc_arc <- calc_dm(input_arc_nz$data_loaded)
set.seed(1150)
adonis2(bc_arc ~ Environment, data = input_arc_nz$map_loaded) # R2 = 0.16, p = 0.001
anova(betadisper(bc_arc, input_arc_nz$map_loaded$Environment)) # Dispersion not homogeneous
pcoa_arc <- cmdscale(bc_arc, k = nrow(input_arc_nz$map_loaded) - 1, eig = T)
eigenvals(pcoa_arc)/sum(eigenvals(pcoa_arc)) # 21.6, 13.7 % variation explained
input_arc_nz$map_loaded$Axis01 <- scores(pcoa_arc)[,1]
input_arc_nz$map_loaded$Axis02 <- scores(pcoa_arc)[,2]
micro.hulls <- ddply(input_arc_nz$map_loaded, c("Environment"), find_hull)
g_arc <- ggplot(input_arc_nz$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F, size = 0.25) +
  geom_point(size = 1, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = "PC1: 21.6%", 
       y = "PC2: 13.7%") +
  ggtitle("Archaea") +
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
nrow(input_bac$data_loaded) # 2679
input_bac_CPM <- input_bac
for (i in 1:ncol(input_bac$data_loaded)) {
  for (j in 1:nrow(input_bac$data_loaded)) {
    input_bac_CPM$data_loaded[j, i] <- (input_bac$data_loaded[j, i]*1000000)/input_bac$map_loaded$GenomeSize[i]
  }
}

# Check zeroes - none, all samples had bacteria
countbac <- as.data.frame(sort(colSums(input_bac$data_loaded))) %>%
  filter(`sort(colSums(input_bac$data_loaded))` == 0)

# BC, PERMANOVA, PERMDISP, PCoA
bc_bac <- calc_dm(input_bac$data_loaded)
set.seed(1150)
adonis2(bc_bac ~ Environment, data = input_bac$map_loaded) # R2 = 0.14, p = 0.001
anova(betadisper(bc_bac, input_bac$map_loaded$Environment)) # Dispersion not homogeneous
pcoa_bac <- cmdscale(bc_bac, k = nrow(input_bac$map_loaded) - 1, eig = T)
eigenvals(pcoa_bac)/sum(eigenvals(pcoa_bac)) # 25.6, 14.7 % variation explained
input_bac$map_loaded$Axis01 <- scores(pcoa_bac)[,1]
input_bac$map_loaded$Axis02 <- scores(pcoa_bac)[,2]
micro.hulls <- ddply(input_bac$map_loaded, c("Environment"), find_hull)
g_bac <- ggplot(input_bac$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Environment, fill = Environment),
               alpha = 0.1, show.legend = F, size = 0.25) +
  geom_point(size = 1, alpha = 0.5, aes(colour = Environment),
             show.legend = T) +
  labs(x = "PC1: 25.6%", 
       y = "PC2: 14.7%") +
  ggtitle("Bacteria") +
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
  labs(x = "PC1: 12.2%", 
       y = "PC2: 9.5%") +
  ggtitle("Fungi") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, vjust = -1))
g

# Multipanel
pdf("Figs/PCoA_ArcBacFun.pdf", width = 8, height = 5)
plot_grid(g_arc,g_bac,g,g_leg, ncol = 2, hjust = "hv")
dev.off()


#### _by Ecosystem ####
# Subset the data into each environment
# Make pie charts of taxa
# Use fungi only input data - get relative abundances of fungal phyla out of just fungi
# Use samples with at least 1 fungal count
input_fungi_nz <- filter_data(input_fungi,
                              filter_cat = "sampleID",
                              filter_vals = rownames(countFun))
input_fungi_nz$map_loaded$Study.Name <- as.factor(input_fungi_nz$map_loaded$Study.Name)
levels(input_fungi_nz$map_loaded$Environment)

# Need to do several for loops
# For subseting, summarizing, plotting, store dfs in a list to enable for loop/indexing
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
length(levels(env[[1]]$map_loaded$Study.Name))
length(levels(env[[2]]$map_loaded$Study.Name))
length(levels(env[[3]]$map_loaded$Study.Name))
length(levels(env[[4]]$map_loaded$Study.Name))
length(levels(env[[5]]$map_loaded$Study.Name))
length(levels(env[[6]]$map_loaded$Study.Name))
length(levels(env[[7]]$map_loaded$Study.Name))
length(levels(env[[8]]$map_loaded$Study.Name))

# Probably best to show some of the variability within environment type
# For example, Shu and Huang 2022 have 4-5 sites for each environment
# Here let's do 1-6, for env. with more than 6 studies, take 6 with greatest sample size
for (i in 1:length(env)) {
  studyN[[i]] <- as.data.frame(table(env[[i]]$map_loaded$Study.Name)) %>%
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
  df[[c]] <- filter_data(env[[i]],
                         filter_cat = "Study.Name",
                         keep_vals = studyN[[i]]$Var1[l])
  counter <- counter + 1
  }
}

# Location and n (for plot titles), do manually
df[[1]]$map_loaded$Location <- "Richmond Mine\nUSA\n(n = 17)"
df[[2]]$map_loaded$Location <- "Malanjkhand copper mine\nIndia\n(n = 10)"
df[[3]]$map_loaded$Location <- "Los Rueldos mercury mine\nSpain\n(n = 3)"
df[[4]]$map_loaded$Location <- "Various Mines\nChina\n(n = 3)"
df[[5]]$map_loaded$Location <- "Rothamsted Park\nUK\n(n = 3)"
df[[6]]$map_loaded$Location <- "Richmond Mine\nUSA\n(n = 1)"
df[[7]]$map_loaded$Location <- "Glacial meltwater/mats\nAntarctica\n(n = 17)"
df[[8]]$map_loaded$Location <- "Cryoconites\nGreenland\n(n = 12)"
df[[9]]$map_loaded$Location <- "Glacier sediment\nAntarctica\n(n = 6)"
df[[10]]$map_loaded$Location <- "Glacial ice\nCanada\n(n = 4)"
df[[11]]$map_loaded$Location <- "Glacial sediment\nCanada/Iceland\n(n = 3)"
df[[12]]$map_loaded$Location <- "Cryoconite\nItaly\n(n = 2)"
df[[13]]$map_loaded$Location <- "Polar Desert soil\nAntarctica\n(n = 33)"
df[[14]]$map_loaded$Location <- "Glacial forefield soil\nSweden/Norway/Greenland\n(n = 64)"
df[[15]]$map_loaded$Location <- "Yellowstone hot springs\nUSA\n(n = 42)"
df[[16]]$map_loaded$Location <- "Various hot springs\nUSA/Canada/China/South Africa\n(n = 39)"
df[[17]]$map_loaded$Location <- "Yellowstone hot spring mats\nUSA\n(n = 34)"
df[[18]]$map_loaded$Location <- "Waikite Valley hot spring mats\nNew Zealand\n(n = 18)"
df[[19]]$map_loaded$Location <- "Yellowstone hot spring sediment\nUSA\n(n = 11)"
df[[20]]$map_loaded$Location <- "Great Boiling Spring sediment\nUSA\n(n = 9)"
df[[21]]$map_loaded$Location <- "Guaymas Basin sediment/mats\nMexico\n(n = 47)"
df[[22]]$map_loaded$Location <- "Mid Cayman Rise\n(n = 27)"
df[[23]]$map_loaded$Location <- "Various vents\nPacific/Atlantic\n(n = 21)"
df[[24]]$map_loaded$Location <- "Mid Cayman Rise plume\n(n = 14)"
df[[25]]$map_loaded$Location <- "Various vents\nPacific\n(n = 13)"
df[[26]]$map_loaded$Location <- "Axial seamount\n(n = 12)"
df[[27]]$map_loaded$Location <- "Various lakes\nAustralia\n(n = 117)"
df[[28]]$map_loaded$Location <- "Organic Lake\nAntarctica\n(n = 34)"
df[[29]]$map_loaded$Location <- "Unrestored salterns\nUSA\n(n = 12)"
df[[30]]$map_loaded$Location <- "Salton Sea\nUSA\n(n = 10)"
df[[31]]$map_loaded$Location <- "Salterns\nNamibia\n(n = 7)"
df[[32]]$map_loaded$Location <- "Bras del Port saltern\nSpain\n(n = 6)"
df[[33]]$map_loaded$Location <- "Alkaline sediment\nRussia/Germany\n(n = 11)"
df[[34]]$map_loaded$Location <- "Alkaline water\nItaly/Philippines/Costa Rica\n(n = 8)"
df[[35]]$map_loaded$Location <- "Mine pit pond\nUSA\n(n = 4)"
df[[36]]$map_loaded$Location <- "Bras del Port saltern\nSpain\n(n = 1)"
df[[37]]$map_loaded$Location <- "Tuz Lake\nTurkey\n(n = 1)"
df[[38]]$map_loaded$Location <- "Diamante Lake biofilm\nAgrgentina\n(n = 1)"

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
  p[[i]] <- plot_taxa_bars(phy[[i]], df[[i]]$map_loaded, "Study.Name", 20) +
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
p_forleg <- plot_taxa_bars(phy[[7]], df[[7]]$map_loaded, "Study.Name", 20) +
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
pies <- plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]],
                  p[[7]], p[[8]], p[[9]], p[[10]], p[[11]], p[[12]],
                  p[[13]], NULL, NULL, NULL, NULL, NULL,
                  p[[14]], NULL, NULL, NULL, NULL, NULL,
                  p[[15]], p[[16]], p[[17]], p[[18]], p[[19]], p[[20]],
                  p[[21]], p[[22]], p[[23]], p[[24]], p[[25]], p[[26]],
                  p[[27]], p[[28]], p[[29]], p[[30]], p[[31]], p[[32]],
                  p[[33]], p[[34]], p[[35]], p[[36]], p[[37]], p[[38]],
                  ncol = 6)
pies

# Add legend
pdf("Figs/Pies_Phyla.pdf", width = 8.5, height = 6.5)
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
  c[[i]] <- plot_taxa_bars(cla[[i]], df[[i]]$map_loaded, "Study.Name", 20) +
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
c_forleg <- plot_taxa_bars(cla[[7]], df[[7]]$map_loaded, "Study.Name", 37) +
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
pies_c <- plot_grid(c[[1]], c[[2]], c[[3]], c[[4]], c[[5]], c[[6]],
                    c[[7]], c[[8]], c[[9]], c[[10]], c[[11]], c[[12]],
                    c[[13]], NULL, NULL, NULL, NULL, NULL,
                    c[[14]], NULL, NULL, NULL, NULL, NULL,
                    c[[15]], c[[16]], c[[17]], c[[18]], c[[19]], c[[20]],
                    c[[21]], c[[22]], c[[23]], c[[24]], c[[25]], c[[26]],
                    c[[27]], c[[28]], c[[29]], c[[30]], c[[31]], c[[32]],
                    c[[33]], c[[34]], c[[35]], c[[36]], c[[37]], c[[38]],
                    ncol = 6)
pies_c

# Add legend
pdf("Figs/Pies_Classes.pdf", width = 8.5, height = 6.5)
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
leveneTest(input_fungi$map_loaded$rich ~ input_fungi$map_loaded$Environment)
m <- aov(rich ~ Environment, data = input_fungi$map_loaded)
summary(m)
shapiro.test(m$residuals)
tuk <- emmeans(object = m, specs = "Environment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(input_fungi$map_loaded$rich)+(max(input_fungi$map_loaded$rich)-min(input_fungi$map_loaded$rich))/20)

leveneTest(input_fungi$map_loaded$shannon ~ input_fungi$map_loaded$Environment)
m1 <- aov(shannon ~ Environment, data = input_fungi$map_loaded)
shapiro.test(m1$residuals)
summary(m1)
tuk1 <- emmeans(object = m1, specs = "Environment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         y = max(input_fungi$map_loaded$shannon)+(max(input_fungi$map_loaded$shannon)-min(input_fungi$map_loaded$shannon))/20)
label_df <- rbind(tuk, tuk1)
facet_df <- c("rich" = "(a) Genus richness",
              "shannon" = "(b) Genus Shannon")
alpha_long <- input_fungi$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
pdf("Figs/Alpha.pdf", width = 6, height = 3)
ggplot(alpha_long, aes(reorder(Environment, value, mean), value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.5, alpha = 0.2, width = 0.3) +
  geom_text(data = label_df, aes(Environment, y, label = str_trim(.group)), 
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



#### ..........................####
#### Functional ####
# Sent Dongying Wu of IMG staff list of taxonoids and list of fungal phyla
# Dongying ran custom python scripts on JGI super computer to pull out KOs of only fungal phyla
# Folder FungalKOs has a file for each metagenome with the KO hits of the fungal phyla scaffolds
# Already deleted 318 blank files (no fungal KOs); 837 had at least 1 KO
# Note that there is bias in eukaryote gene calling/KO assignment
# Run a for loop to read in the file for each metagenome and combine into 1
setwd("FungalKOs/")
ko <- list()
ko_input <- data.frame(V1 = "NA",
                       V2 = "NA",
                       V3 = "NA")
for (i in 1:length(list.files())) {
  ko[[i]] <- read.delim(list.files()[i], header = F)
  ko_table <- ko_table %>%
    rbind(ko[[i]])
}

setwd("~/Documents/GitHub/Extremophilic_Fungi/")

# Clean up table
ko_table_wTax <- ko_input %>%
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
write.csv(ko_table$KO, file = "KOlist.csv", row.names = F)

# List of KOs sorted by overall abundance
ko_list <- ko_input %>%
  filter(V1 != "NA") %>%
  dplyr::select(V1) %>%
  group_by(V1) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(Definition = "NA") %>%
  separate(V1, into = c("Junk", "KO"), sep = ":", remove = F)

# Add definitions to list
for (i in 1:nrow(ko_list)) {
  def <- keggFind(database = "ko", query = ko_list$KO[i])
  if (length(def) != 0) {
    ko_list$Definition[i] <- def
  }
}

# Make community style table and metadata, match IDs
ko_comm <- ko_table_MGcount %>%
  t() %>%
  as.data.frame() %>%
  filter(rownames(.) %in% input_fungi$map_loaded$taxon_oid) %>%
  arrange(rownames(.))

ko_meta <- input_fungi$map_loaded %>%
  filter(taxon_oid %in% rownames(ko_comm)) %>%
  arrange(taxon_oid)

# Check match (should be zero)
sum(rownames(ko_comm) != ko_meta$taxon_oid)

# Check environment sample size
table(ko_meta$Environment)

# Get richness
ko_meta$richness_KO = specnumber(ko_comm)
range(ko_meta$richness_KO)
pdf("Figs/KO_Genus_richness.pdf", width = 6, height = 4)
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
#                              design = ~ Study.Name)
# dds <- estimateSizeFactors(dds)
# dds <- estimateDispersions(dds)
# ko_comm_DESeq <- as.data.frame(t(counts(dds, normalized = T)))
# Save so you don't have to redo the DESeq (takes a while)
# saveRDS(ko_comm_DESeq, "ko_comm_DESeq.rds")
ko_comm_DESeq <- readRDS("ko_comm_DESeq.rds")

#### _KO Richness ####
leveneTest(richness_KO ~ Environment, data = ko_meta)
m <- aov(richness_KO ~ Environment, data = ko_meta)
summary(m)
shapiro.test(m$residuals)
TukeyHSD(m)
kruskal.test(richness_KO ~ Environment, data = ko_meta)
kwAllPairsNemenyiTest(richness_KO ~ Environment, data = ko_meta)
pdf("Figs/KO_richness.pdf", width = 6, height = 3)
ggplot(ko_meta, aes(reorder(Environment, richness_KO, mean), richness_KO)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1, alpha = 0.2, width = 0.3) +
  labs(x = "Environment", 
       y = "# KOs") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        strip.text = element_text(size = 10))
dev.off()



#### _KO Community ####
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

plot_grid(g1_ko, g2_ko, align = "hv", ncol = 2, rel_widths = c(1,1.515))

# These are incredibly similar and also show a dense cluster with barely any dissimilarity
# I suspect this is driven by KO count, metagenomes with only 1 KO can't be too dissimilar
# Check number of ones
ko_richness <- ko_meta %>%
  dplyr::select(richness_KO) %>%
  arrange(desc(richness_KO)) %>%
  mutate(index = seq(1:nrow(ko_richness)))
plot(rownames(ko_list), ko_list$n)
plot(ko_richness$index, ko_richness$richness_KO)
sum(ko_richness$richness_KO == 1) # 77 with just 1 KO
sum(ko_richness$richness_KO == 2) # 51 with just 2 KOs
sum(ko_richness$richness_KO == 3) # 55 with just 3 KOs

# Need to find good cutoff with some KOs and still high sample size
sum(ko_richness$richness_KO > 10) # 492 with > 10

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
sum(ko_richness$richness_KO > 100) # 161 with > 100
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

pdf("Figs/KO_PCoA.pdf", width = 8.5, height = 3.5)
plot_grid(g5_ko, g6_ko, ko_pcoa_leg, align = "hv", ncol = 3, rel_widths = c(2.5,2.5,1))
dev.off()


#### __Stats ####
set.seed(308)
adonis2(bc_ko ~ Environment, data = ko_meta_filt) # R2 = 0.18, p = 0.001
anova(betadisper(bc_ko, ko_meta_filt$Environment)) # Dispersion not homogeneous

set.seed(308)
adonis2(jac_ko ~ ko_meta_filt$Environment) # R2 = 0.17, p = 0.001
anova(betadisper(jac_ko, ko_meta_filt$Environment)) # Dispersion not homogeneous



#### __Drivers/Indicators ####
# KO MULTIPATT (list KOs associated with each group)
set.seed(425)
mp <- multipatt(ko_comm_DESeq, 
                ko_meta$Environment, 
                func = "IndVal.g", 
                control = how(nperm=999))
summary(mp) # None!!

#### _ Specific KOs ####
# Extract and analyze list of KOs of interest
# Need to get list of KOs from Lara
# For now make barplot and heatmap of top KOs
 

#### ..........................####
#### Other ####
#### _Map ####
# Sample map with ggplot
world <- map_data("world")
input$map_loaded$Latitude <- as.numeric(input$map_loaded$Latitude)
input$map_loaded$Longitude <- as.numeric(input$map_loaded$Longitude)

# 21 missing
pdf("Figs/SampleMap.pdf", width = 8, height = 5)
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


