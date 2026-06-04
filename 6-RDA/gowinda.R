## revision 1: gowinda

setwd("/Users/lilydurkee/colias2025/")

library(ggplot2)
library(dplyr)
library(tidyverse)

# 1. Load the Gowinda Results
gowinda <- read.table("gowinda_results.txt", sep="\t", quote="", fill=TRUE, stringsAsFactors=FALSE)
colnames(gowinda) <- c("GO_ID", "Expectation", "Observed", "P_value", "FDR", 
                  "n_genes_sig", "n_genes_total", "n_genes_univ", "extra", "Merged_Col")

# 2. Load the Name Key
gene_key <- read.table("gene_names_key.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
colnames(gene_key) <- c("Gene_ID", "Description")

# 3. Clean the 'smashed' column to isolate Gene IDs
# This separates the 'Biological_Process' prefix and splits genes into long format
sig_genes_mapped <- gowinda %>%
  # Remove the prefixes (Biological_Process, etc.)
  mutate(Gene_List = gsub("Biological_Process|Molecular_Function|Cellular_Component", "", Merged_Col)) %>%
  # Split the comma-separated string into individual rows
  separate_rows(Gene_List, sep = ",") %>%
  # IMPORTANT: Force the IDs to UPPERCASE to match your gene_key
  mutate(Gene_List = toupper(Gene_List)) %>%
  # Now join with the gene_key (ensure gene_key$Gene_ID is also uppercase)
  left_join(gene_key %>% mutate(Gene_ID = toupper(Gene_ID)), 
            by = c("Gene_List" = "Gene_ID"))

# This should now show names instead of NAs
head(sig_genes_mapped, 10)

#### map GO to function ####

#install.packages("BiocManager")
#BiocManager::install("GO.db")
library(GO.db)

# Create a mapping of IDs to Terms
goterms <- Term(GOTERM)

# Add the terms to your data
sig_genes_mapped <- sig_genes_mapped %>%
  mutate(GO_Term_Name = goterms[GO_ID])


# 4. Create the Top Hits Summary for your Table
# This groups by GO term and lists the top 3 protein names for each
summary_table <- sig_genes_mapped %>%
  group_by(GO_ID, GO_Term_Name, P_value, FDR, Expectation, Observed) %>%
  summarize(Top_Genes = paste(unique(na.omit(Description))[1:3], collapse = "; "), .groups = 'drop') %>%
  arrange(P_value)

#write.csv(summary_table, "summary_table-gowinda-update.csv")
# Now when you plot, use 'GO_Term_Name' for the Y-axis

# 5. Create the Dot Plot (using GO IDs or names)
# if you want GO names, you can add the GO.db mapping here too

# wrap the long GO names (manually inserting a newline every 45 characters)
# This prevents the Y-axis from getting too wide
summary_table$GO_Wrapped <- gsub("(.{1,45})(\\s|$)", "\\1\n", summary_table$GO_Term_Name)

# create the plot using standard R expressions for the title
ggplot(head(summary_table, 15), aes(x = Observed, y = reorder(GO_Wrapped, Observed))) +
  geom_point(aes(size = Observed, color = P_value)) +
  
  # Color scale: Red for significance
  scale_color_gradient(low = "firebrick", high = "dodgerblue") +
  
  theme_bw() +
  
  # Use expression() to italicize the species name without extra packages
  labs(
    title = expression(paste("Top enriched pathways in ", italic("C. p. eriphyle"))),
    subtitle = "Significant genes (p-value < 0.05)",
    x = "Number of significant genes hit",
    y = NULL, # Removed 'Gene Ontology ID' label to save space
    size = "count",
    color = "p-value"
  ) +
  
  # Adjust text size and spacing
  theme(
    axis.text.y = element_text(size = 9, lineheight = 0.8),
    plot.title = element_text(size = 14, face = "bold")
  )

#### g:Profiler prep ####
library(dplyr)
library(tidyr)
library(GO.db)

# 1. Load your 3-column map
# Use 'fill=TRUE' in case some lines are slightly different
go_map <- read.table("colias_go_map-FINAL7.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE, quote="")
colnames(go_map) <- c("GO_ID", "Category", "Gene_ID")

# 2. Fetch the REAL Biological Process names (instead of just 'Biological_Process')
goterms <- Term(GOTERM)
go_map$Term_Name <- goterms[go_map$GO_ID]

# Handle any missing names
go_map$Term_Name[is.na(go_map$Term_Name)] <- "Unknown_Process"

# 3. Collapse the data into the "Wide" GMT format
# This puts all genes for a GO ID into one single tab-separated row
gmt_ready <- go_map %>%
  group_by(GO_ID, Term_Name) %>%
  summarize(Gene_List = paste(unique(Gene_ID), collapse = "\t"), .groups = 'drop')


# 4. Save the file
# This format: ID <tab> Name <tab> Gene1 <tab> Gene2...

# write.table(gmt_ready, "Colias_eriphyle_custom.gmt", 
#             sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# 

#### g:Profiler results ####

# Load libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

gProf <- read.csv("gProfiler_custom_2-25-2026.csv") 

#### Version 1 ####

# # 1. Create the plot data
# gProf.plot_data <- gProf %>%
#   # Use .data$ to tell R to look for columns, not functions
#   arrange(source, adjusted_p_value) %>%
#   mutate(
#     term_index = row_number(),
#     # Create a cleaner label for the plot
#     short_name = ifelse(nchar(term_name) > 30, 
#                         paste0(substr(term_name, 1, 27), "..."), 
#                         term_name)
#   )
# 
# # 2. Identify top terms for labeling (Top 2 per source)
# # This ensures we label the most significant points in each category
# top_labels <- gProf.plot_data %>%
#   group_by(source) %>%
#   slice_max(order_by = negative_log10_of_adjusted_p_value, n = 2) %>%
#   ungroup()
# 
# # 3. Calculate X-axis midpoints for the source labels
# source_midpoints <- gProf.plot_data %>%
#   group_by(source) %>%
#   summarize(mid = mean(term_index))

##### Version 2 ####

# 
# # 1. Process data with smarter grouping
# gProf.plot_data <- gProf %>%
#   mutate(
#     # Identify sub-groups within your custom source
#     group_source = case_when(
#       grepl("GO:0008150|GO:0003674|GO:0005575", term_id) ~ "GO:Root", # Catch-all
#       grepl("GO:0005", term_id) ~ "GO:CC", # Cellular Component
#       grepl("GO:0008", term_id) ~ "GO:MF", # Molecular Function
#       grepl("GO:000", term_id)  ~ "GO:BP", # Biological Process
#       TRUE ~ gsub(":.*", "", term_id)      # Fallback: extract prefix before ":"
#     )
#   ) %>%
#   # Sort by our new group_source
#   arrange(group_source, adjusted_p_value) %>%
#   mutate(
#     term_index = row_number(),
#     short_name = ifelse(nchar(term_name) > 30, 
#                         paste0(substr(term_name, 1, 27), "..."), 
#                         term_name)
#   )
# 
# # 2. Pick more labels (Top 5 per group)
# top_labels <- gProf.plot_data %>%
#   group_by(group_source) %>%
#   slice_max(order_by = negative_log10_of_adjusted_p_value, n = 5) %>%
#   ungroup()
# 
# # 3. Calculate midpoints for the new groups
# group_midpoints <- gProf.plot_data %>%
#   group_by(group_source) %>%
#   summarize(mid = mean(term_index))
# 
# # 4. Plot
# ggplot(gProf.plot_data, aes(x = term_index, y = negative_log10_of_adjusted_p_value, color = group_source)) +
#   geom_point(aes(size = intersection_size), alpha = 0.5) +
#   geom_text_repel(
#     data = top_labels,
#     aes(label = short_name),
#     size = 2.8,
#     fontface = "italic",
#     box.padding = 0.5,
#     max.overlaps = 30, # Increased to allow more labels
#     segment.color = 'grey50'
#   ) +
#   scale_x_continuous(breaks = group_midpoints$mid, labels = group_midpoints$group_source) +
#   scale_color_brewer(palette = "Dark2") + 
#   labs(
#     title = "Functional Enrichment: Colias p. eriphyle",
#     x = "Functional Category (Derived from Term ID)",
#     y = "-log10(Adj. P-value)",
#     size = "Gene Count",
#     color = "Category"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 0, hjust = 0.5), # Standard orientation often looks cleaner here
#     panel.grid.minor = element_blank()
#   )

#### Version 3 ####
library(ggplot2)
library(dplyr)
library(ggrepel)
library(stringr)

# 1. Process and Filter Data
gProf.plot_data <- gProf %>%
  mutate(
    group_source = case_when(
      grepl("GO:0005", term_id) ~ "Cellular Component",
      grepl("GO:0008", term_id) ~ "Molecular Function",
      grepl("GO:000", term_id)  ~ "Biological Process",
      TRUE ~ "Other GO"
    )
  ) %>%
  # Filter out very broad terms (e.g., > 1000 genes) to keep it specific
  #filter(term_size < 1000) %>% 
  arrange(group_source, adjusted_p_value) %>%
  mutate(
    term_index = row_number(),
    # Wrap text to 20 chars for vertical stacking
    wrapped_name = str_wrap(term_name, width = 20)
  )

# 2. SET THRESHOLD HERE
# Try 8 or 10. If too many show up, raise it to 12.
sig_threshold <- 5

top_labels <- gProf.plot_data %>%
  filter(negative_log10_of_adjusted_p_value > sig_threshold)

# Quick check: How many labels will show up?
print(paste("Number of labels to be plotted:", nrow(top_labels)))

# 3. Calculate X-axis midpoints
group_midpoints <- gProf.plot_data %>%
  group_by(group_source) %>%
  summarize(mid = mean(term_index))

# 4. Final Plot
ggplot(gProf.plot_data, aes(x = term_index, y = negative_log10_of_adjusted_p_value, color = group_source)) +
  # Use a smaller alpha (0.3) so the background dots don't distract from labels
  geom_point(aes(size = intersection_size), alpha = 0.3) +
  geom_text_repel(
    data = top_labels,
    aes(label = wrapped_name),
    size = 2.8,
    lineheight = 0.8,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = 50, 
    force = 5,          # Higher force pushes labels away from crowded dots
    segment.alpha = 0.5 # Makes the connector lines subtle
  ) +
  scale_x_continuous(breaks = group_midpoints$mid, labels = group_midpoints$group_source) +
  scale_color_brewer(palette = "Dark2") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + # Extra room at the top
  labs(
    title = "gProfiler: Functional Gene Enrichment",
    #subtitle = paste("Colias p. eriphyle (-log10 p-value >", sig_threshold, ")"),
    x = "Gene Ontology Category",
    y = "-log10(Adj. p-value)",
    size = "Genes",
    color = "Category"
  ) +
  theme_classic() +
  theme(
    # --- REMOVE X-AXIS LABELS AND TICKS ---
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # --------------------------------------
    #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold", size = 10),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

#### COMPARISON OF PCAs: ANGSD vs GATK vs PSEUDO-HAP ####

## angsd

# 1. Load the Covariance Matrix
cov_mat <- as.matrix(read.table("colias-thinned-1x-PCA.cov"))

# 2. Extract Eigenvalues/Vectors
eig <- eigen(cov_mat)

# 3. Create a Dataframe for Plotting
# Note: Ensure bams.filelist is in the same folder
raw_paths <- readLines("bams.filelist")
clean_filenames <- basename(raw_paths)
sample_ids <- gsub("[_\\.].*", "", clean_filenames)

pca_angsd <- data.frame(
  ID = sample_ids,
  PC1 = eig$vectors[,1],
  PC2 = eig$vectors[,2],
  site = substr(sample_ids, 1, 2)
)

site_info <- read.csv("env_by_site.csv")

pca_angsd <- merge(pca_angsd, site_info, by="site")

# 4. Simple Plot
library(ggplot2)

p1 <- ggplot(pca_angsd, aes(x=PC1, y=PC2, fill=pair, shape=cat)) +
  geom_point(size=3) +
  theme_bw() +
  labs(x="PC1", y="PC2", title="Single-allele sampling (ANGSD)") +
  scale_fill_manual(values=c("#031A6B", "dodgerblue", "lightblue2",  
                             "#750D37", "firebrick1", "#ff9896",
                             "darkgray", 
                             "#AA8EF4"), guide="none")+
  scale_shape_manual(values = c(24, 21))

## gatk (publbished PCA)

pca_gatk <- read.csv("colias_all_PCA.csv")
pca_gatk$site <- substr(pca_gatk$ID, 1, 2)

pca_gatk <- merge(pca_gatk, site_info, by="site")

p2 <- ggplot(pca_gatk, aes(x=-PC1, y=PC2, fill=pair, shape=cat)) +
  geom_point(size=3) +
  theme_bw() +
  labs(x="PC1", y="PC2", title = "Called genotypes (GATK)") +
  scale_fill_manual(values=c("#031A6B", "dodgerblue", "lightblue2",  
                             "#750D37", "firebrick1", "#ff9896",
                             "darkgray",
                             "#AA8EF4"), guide="none")+
  scale_shape_manual(values = c(24, 21))

library(ggpubr)

ggarrange(p1, p2, #p3,
          nrow = 1,
          common.legend = TRUE, # Use a common legend
          legend = "bottom",    # Position the common legend
          labels = c("A)", "B)"), # Add labels
          align = "hv"          # Align plots horizontally and vertically
)


## stats
library(vegan)

# Reorder GATK to match the ANGSD order exactly
pca_gatk_ordered <- pca_gatk[match(sample_ids, pca_gatk$ID), ]
pca_angsd_ordered <- pca_angsd[match(sample_ids, pca_angsd$ID), ]

# Statistical Test 1: GATK vs. Likelihoods
pro_beagle <- protest(pca_gatk_ordered[, c("PC1", "PC2")], pca_angsd_ordered[, c("PC1", "PC2")], permutations = 999)
summary(pro_beagle)
pro_beagle

cat("Correlation (GATK vs Likelihoods):", pro_beagle$t0, "\n")

#### top RDA SNPs/genes table ####

library(dplyr)
library(tidyverse)

# Load your full dataset
genes <- as.data.frame(read.csv("sig.genes-all.csv"))
genes$abs_loading <- abs(genes$loading)

# 1. Calculate absolute loading and select top 50
top_50_genes <- genes %>%
  filter(d < 500) %>%
  mutate(abs_loading = abs(loading)) %>%
  arrange(desc(abs_loading)) %>%
  dplyr::slice(1:50) %>%
  # 2. Select and rename columns for a professional look
  dplyr::select(
    GeneName = GeneName,
    Product = product,
    Chrom = CHROM,
    Position = POS,
    Effect = type,
    loc = loc,
    RDA_axis = axis,
    RDA_loading = loading,
    TopPredictor = predictor,
    corr = correlation,
    dist = d
  )

# 3. Clean up the 'Product' names (removing %2C etc.)
top_50_genes$Product <- gsub("%2C", ",", top_50_genes$Product)

# Export for Supplementary Materials
#write.csv(top_50_genes, "SupplementaryTable2_Top50_SNP-RDA.csv", row.names = FALSE)

#### load correlations ####
env_corr <- read.csv("RDA.freq_env.results.csv")

top_50_genes.env <- merge(top_50_genes, env_corr, by="loc", all=F)
write.csv(top_50_genes.env, "SupplementaryTable2_Top50SNP_RDA_env.csv", row.names = FALSE)

#### MAF FDR correction #####
# 1. Load your data (replace 'your_file.csv' with your actual filename)
df <- read.csv("RDA.freq_env.results.csv")

# 2. Calculate q-values (FDR) for each environmental variable
# We use the Benjamini-Hochberg ("BH") method
df$elev_qval   <- p.adjust(df$elev_p.val, method = "BH")
df$precip_qval <- p.adjust(df$precip_p.val, method = "BH")
df$srad_qval   <- p.adjust(df$srad_p.val, method = "BH")

# 3. Create a logical column for Criterion #3 
# (Did the SNP pass FDR < 0.05 for at least one variable?)
df$is_significant <- df$elev_qval < 0.05 | 
  df$precip_qval < 0.05 | 
  df$srad_qval < 0.05

# 4. Filter to see your final candidate list
final_candidates <- subset(df, is_significant == TRUE)

# 5. Export the results
write.csv(df, "RDA_outliers_with_FDR.csv", row.names = FALSE)

#### Gowinda Fisher's exact test ####

# output from gowinda
data <- matrix(c(13011, 943716, 3481, 510908), nrow = 2)
fisher.test(data)

#### NCBI SRA info ####
library(dplyr)
library(tidyr)
library(stringr)

# 1. Load the list from your saved file
files <- readLines("colias-fastq-all.txt")

# 2. Extract pairs and pivot
sra_metadata <- data.frame(full_name = files) %>%
  mutate(
    # Extract ID (e.g., AL2_CKDL230002367-1A_HJJGFCCX2_L6)
    sample_id = str_extract(full_name, "^[^_]+"),
    # Identify which column it goes into (1 -> filename, 2 -> filename2)
    read_num = str_extract(full_name, "(?<=_)[12](?=\\.fq\\.gz$)")
  ) %>%
  pivot_wider(
    names_from = read_num, 
    values_from = full_name, 
    names_prefix = "file_"
  ) %>%
  # 3. Create the exact column structure you provided
  transmute(
    sample_name        = sample_id,
    library_ID         = sample_id,
    title              = "WGS of Colias p. eriphyle: thorax tissue",
    library_strategy   = "WGS",       # Update if this was RNA-Seq or other
    library_source     = "GENOMIC",
    library_selection  = "RANDOM",
    library_layout     = "paired",
    platform           = "ILLUMINA",
    instrument_model   = "Illumina HiSeq 4000",
    design_description = "Standard Novogene library preparation",
    filetype           = "fastq",
    filename           = file_1,
    filename2          = file_2,
    filename3          = "",
    filename4          = "",
    assembly           = "",
    fasta_file         = ""
  )

# 4. Export as a Tab-Delimited file (TSV) for NCBI
# Using na="" ensures the empty columns remain blank for the portal
write.table(sra_metadata, "SRA_Metadata_Final.tsv", 
            sep="\t", row.names=FALSE, quote=FALSE, na="")

write.csv(sra_metadata, "SRA_metadata-final.csv", row.names = F)

#print("Metadata file 'SRA_Metadata_Final.tsv' generated successfully.")

#### NCBI SRA BioSample prep ####
biosample_list <- data.frame(full_name = files) %>%
  mutate(
    # Extract ID: Everything before the first underscore (e.g., AL2)
    sample_name = str_extract(full_name, "^[^_]+"),
    
    # Extract Site: First two letters (e.g., AL)
    site = str_sub(sample_name, 1, 2)
  ) %>%
  # 3. Collapse the two lines (R1 and R2) into one unique line per ID
  distinct(sample_name, site, .keep_all = FALSE)

env <- read.csv("env_by_site.csv")
samples <- read.csv("Colias-samples-2022-sex.all.csv")

biosample_list <- merge(biosample_list, samples, by="sample_name", all.x=T)

# biosample_list <- merge(biosample_list, env, by="site", all=T) %>%
#   dplyr::select(sample_name = sample_name, 
#                 sample_title = paste("Colias p. eriphyle - ", sample_name),
#                 organism = ifelse(site=="AL", "Colias alexandra", ifelse(
#                   site=="Eu", "Colias eurytheme", "Colias philodice eriphyle"
#                 )),
#                 geo_loc_name = paste("USA: ", divide, ": ", site_name),
#                 sex = sex, 
#                 altitude = elevation,
#                 dev_stage = "adult", 
#                 collected_by = "Colorado State University",
#                 lat_long = paste0(round(long, 2), " N ", round(abs(lat), 2), " W")
#   )

biosample_list <- merge(biosample_list, env, by="site", all=T)
 
biosample_list.final <- biosample_list %>% mutate(
    # 1. Create the new columns first
    sample_name = sample_name, # Using your ID column as the unique name
    sample_title = paste("Colias p. eriphyle -", sample_name),
    organism = case_when(
      site == "AL" ~ "Colias alexandra",
      site == "Eu" ~ "Colias eurytheme",
      TRUE         ~ "Colias philodice eriphyle"
    ),
    geo_loc_name = paste0("USA: ", state, ": ", divide, ": ", site.name),
    # Note: latitude is usually 'y' and longitude is 'x'
    # Swap 'lat' and 'long' below if your 'lat' is negative
    lat_lon = paste0(round(lat, 2), " N ", round(abs(long), 2), " W"),
    dev_stage = "adult",
    collected_by = "Colorado State University"
  ) %>%
  # 2. Now select/rename for the final NCBI format
  dplyr::select(
    sample_name, 
    sample_title,
    site_name,
    site,
    organism,
    geo_loc_name,
    sex, 
    altitude = elevation, # Renaming 'elevation' to 'altitude'
    dev_stage, 
    collected_by,
    lat_lon
  )

#write.csv(biosample_list.final, "Biosample_metadata-final.csv", row.names = F)

## re-upload after changes made
biosample_list.final <- read.csv("Biosample_metadata-final.csv") 

## add date
collection_info <- read.csv("Colias-samples-2022-all.csv") %>%
  mutate(site = substr(field_ID, 1, 2)) %>%
  mutate(sample_name = ID, collection_date = as.Date(date, format = "%m/%d/%y")) %>%
  dplyr::select(sample_name, site, collection_date)

biosample_list.final <- merge(biosample_list.final, collection_info, by="sample_name", all.x=T)

#write.csv(biosample_list.final, "Biosample_metadata-final.csv", row.names = F)
