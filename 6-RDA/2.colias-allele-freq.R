## Colias allele freq analysis
setwd("/Users/lilydurkee/OneDrive - Colostate/Grad School/R-Projects-Grad/Colias")

indiv <- read.csv("Data/colias_bioclim-elev.csv")

library(dplyr)
env_by_family <- indiv %>%
  group_by(site, site.1) %>%
  summarize(long = mean(x),
            lat = mean(y),
            tmax = mean(tmax),
            precip = mean(precip),
            srad = mean(srad),
            elevation = mean(elevation))

site_codes <- sites %>% select(site, elevation, pair)
env_by_family <- merge(env_by_family, site_codes, by="site")

#write.csv(file="Data/env_by_family.csv", x=env_by_family)

# Load PLINK allele frequency file
freq <- read.csv("colias_freq_by_pop.csv", header = TRUE, stringsAsFactors = FALSE)

# Load environmental data
env <- read.csv("env_by_site.csv", header = TRUE, stringsAsFactors = FALSE)

freq_env <- merge(freq, env, by="CLST", all=T)

library(dplyr)
library(broom)  # for tidy() output from lm()

# Get list of unique SNPs
snp_list <- unique(freq_env$SNP)

# Initialize empty list to hold results
results <- list()

# Loop through SNPs
for (snp in snp_list) {
  subset_data <- freq_env %>% filter(SNP == snp)
  
  # Skip if missing data or too few rows
 ## if (nrow(subset_data) < 3) next
  
  # Run regressions
  lm_elev <- lm(MAF ~ elevation, data = subset_data)
  lm_precip <- lm(MAF ~ precip, data = subset_data)
  lm_srad <- lm(MAF ~ srad, data = subset_data)
  
  # Tidy results and store
  results[[snp]] <- data.frame(
    SNP = snp,
    elev_p.val = summary(lm_elev)$coefficients["elevation", "Pr(>|t|)"],
    elev_r2 = summary(lm_elev)$r.squared,
    precip_p.val = summary(lm_precip)$coefficients["precip", "Pr(>|t|)"],
    precip_r2 = summary(lm_precip)$r.squared,
    srad_p.val = summary(lm_srad)$coefficients["srad", "Pr(>|t|)"],
    srad_r2 = summary(lm_srad)$r.squared
  )
}

# Combine all results into one data frame
freq_env.results <- do.call(rbind, results)

# write.csv(x=freq_env.results, file="Data/freq_env.results.csv")

# convert to bed file
library(tidyr)
colnames(freq_env.results)
head(freq_env.results)

freq_env.results.bed <- freq_env.results %>%
  mutate(
    scaffold = sub("\\.[^\\.]+$", "", SNP),              # remove the last ".12345"
    pos = sub(".*\\.", "", SNP),                         # keep only the last ".12345"
    pos = as.integer(pos),
    start = pos - 1,
    end = pos
  ) %>%
  select(scaffold, start, end, SNP)


# save to BED
write.table(x=freq_env.results.bed, "freq_env.results.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

library(readr)
# read in bedtools closest analysis
freq_env.gene <- read_tsv("Data/freq_env.results-genes-close.mRNA.bed", col_names = FALSE)
colnames(freq_env.gene)[1:4] <- c("CHR", "START", "END", "SNP")  # adjust if needed
colnames(freq_env.gene)[13] <- "gene_info"
colnames(freq_env.gene)[14] <- "dist"

library(stringr)  
freq_env.gene2 <- freq_env.gene %>%
  separate(gene_info, into = paste0("gene.info", 1:10), sep = ";", fill = "right") %>%
  mutate(across(starts_with("info"), ~ str_trim(.)))

# final file
freq_env.gene.final <- merge(freq_env.gene2, freq_env.results, by="SNP")

#write.csv(x=freq_env.gene.final, "Data/freq_env.gene.csv")

RDA.genes <- read.csv("sig.genes-all.csv") 

RDA.genes.sum <- RDA.genes %>%
  arrange(desc(correlation)) %>%
  distinct(GeneName, .keep_all = TRUE) %>%
  mutate(order = row_number()) %>%
  select(order,product)



RDA.genes$SNP <- paste(RDA.genes$CHROM, RDA.genes$POS, sep = ".")

# merging all genes, this doesn't really work
RDA.freq_env <- merge(freq_env.results, RDA.genes, by="SNP", all=F)

length(unique(freq_env.results$SNP))          # total SNPs in freq-env
length(unique(RDA.genes$SNP))          # total SNPs in RDA
sum(freq_env.results$SNP %in% RDA.genes$SNP)   # total overlapping SNPs

# there is some issue because the merged file is only showing 300 SNPs :(


#### ANALYSIS WITH RDA-IDENTIFIED LOCI ONLY ####
# Load PLINK allele frequency file
RDA.freq <- read.table("Data/colias_freq_by_pop.RDA.frq.strat", header = TRUE)

RDA.freq_env <- merge(RDA.freq, env, by="CLST")

# Get list of unique SNPs
snp_list <- unique(RDA.freq_env$SNP)

# Initialize empty list to hold results
results <- list()

# Loop through SNPs
for (snp in snp_list) {
  subset_data <- RDA.freq_env %>% filter(SNP == snp)
  
  # Skip if missing data or too few rows
  ## if (nrow(subset_data) < 3) next
  
  # Run regressions
  lm_elev <- lm(MAF ~ elevation, data = subset_data)
  lm_precip <- lm(MAF ~ precip, data = subset_data)
  lm_srad <- lm(MAF ~ srad, data = subset_data)
  
  # Tidy results and store
  results[[snp]] <- data.frame(
    SNP = snp,
    elev_p.val = summary(lm_elev)$coefficients["elevation", "Pr(>|t|)"],
    elev_r2 = summary(lm_elev)$r.squared,
    precip_p.val = summary(lm_precip)$coefficients["precip", "Pr(>|t|)"],
    precip_r2 = summary(lm_precip)$r.squared,
    srad_p.val = summary(lm_srad)$coefficients["srad", "Pr(>|t|)"],
    srad_r2 = summary(lm_srad)$r.squared
  )
}

# Combine all results into one data frame
RDA.freq_env.results <- do.call(rbind, results)
write.csv(x=RDA.freq_env.results, file="RDA.freq_env.results.csv")

RDA.genes <- read.csv("sig.genes-all.csv")
RDA.genes$SNP <- RDA.genes$loc

RDA.freq_env.results <- read.csv("RDA.freq_env.results.csv")

RDA.freq_env.all <- merge(RDA.freq_env.results, RDA.genes, by="SNP", all=F)

# count unique gene names
length(unique(RDA.freq_env.all$GeneName))

# graphs
library(patchwork)

#### candidate SNPs/genes to plot ####
plot_list <- list()

snps_to_plot <- c(
#"NC_059548.1.6641660", # drought, ATP binding cassette 
#"NC_059564.1.1288957", # Aromatic L-amino acid decarboxylase involved in melanin pathway
#"NC_059568.1.4445158" (OLD) # TH
"NC_059568.1.4450341", # TH
#"NC_059566.1.2650464", # yellow-f2 gene
"NC_059542.1.4495778", # yellow
#"NC_059558.1.5770039", # PTB 1b (aka PTBN1), hypoxia
"NC_059541.1.7361965", #hsp 70
"NC_059555.1.678453"# cytochrome P450 (CYP6B6), insecticide resistance/adaptation
# NC_059555.1.678527 (OLD)
#"NC_059543.1.5418330", # cytochrome P450 (CYP6K1)
#"NC_059543.1.5420138" (OLD), 
#"NC_059543.1.9407508" #LIP3, starvation/aging
#"NC_059568.1.3041721" #pickpocket gene (PKP), painful stimulii

#"NC_059562.1.2931363" #hsp83
)

plot_titles <- c("TH", "yellow", "HSP70","CYP6B6")
                # "CYP6K1", 
                 #"HSP90")

# subset of data frame
RDA.freq_env.sub <- subset(RDA.freq_env.all, SNP %in% snps_to_plot) %>%
  select(-gene_info)

library(ggplot2)
count <- 1

for(snp in snps_to_plot) {
  df_snp <- RDA.freq_env %>% filter(SNP == snp) 
  
  p <- ggplot(df_snp, aes(x = precip, y = MAF)) +
    geom_point(aes(fill = pair, shape=cat), size = 3) +
    geom_smooth(method = "lm", 
                color="black", size=.75, se=F, linetype="dashed") +
    theme_minimal()+
    theme(plot.title = element_text(face = "italic"))+
    labs(title=plot_titles[[count]], x="Precipitation")+
    scale_shape_manual(values = c(24, 21))+
    scale_fill_manual(values=c("#031A6B", "dodgerblue", "lightblue2",  
                               "#750D37", "firebrick1", "#ff9896",
                               "#AA8EF4"), guide="none")+
    labs(shape="elevation") +
    theme(axis.title = element_blank(),)
  
  
  plot_list[[snp]] <- p
  count <- count+1
  
 # ggsave(paste0("plot_", snp, ".png"), plot = p, width = 6, height = 4)
}

# Combine all plots in a grid
#wrap_plots(plot_list, ncol = 3) +
#  plot_layout(guides = "collect")

library(ggpubr)
plot_grid <- ggarrange(plotlist = plot_list,
                  ncol = 2, nrow = 2,
                  common.legend=TRUE,
                  legend="none")

# Add common axis labels
annotate_figure(plot_grid, 
  bottom = text_grob("Precipitation", size = 16),
  left = text_grob("MAF", size = 16, rot=90))

#### longitude vs. precipitation ####
long <- lm(precip ~ long, data=env)
summary(long)

ggplot(data=env, aes(x=long, y=precip))+
  geom_line(aes(group=pair, color=pair), alpha=.5)+
  geom_point(aes(fill=pair, shape=cat), size = 3) +
  geom_smooth(method = "lm", 
              color="black", size=.75, se=F, linetype="dashed") +

  theme_minimal()+
  scale_shape_manual(values = c(24, 21))+
  scale_fill_manual(values=c("#031A6B", "dodgerblue", "lightblue2",  
                             "#750D37", "firebrick1", "#ff9896",
                             "#AA8EF4"), guide="none")+
  scale_color_manual(values=c("#031A6B", "dodgerblue", "lightblue2",  
                              "#750D37", "firebrick1", "#ff9896",
                              "#AA8EF4"), guide="none")
  

#### gene summary ####
RDA.freq_env.sum <- RDA.freq_env.all %>%
  group_by(GeneName) %>%
  summarize(d = mean(d),
            type=type,
            modifier=modifier,
            predictor = predictor,
            corr.RDA = mean(correlation),
            elev_p.val = mean(elev_p.val),
            precip_p.val = mean(precip_p.val),
            srad_p.val = mean(srad_p.val)) %>%
  arrange(desc(corr.RDA)) %>%
  distinct(GeneName, .keep_all = TRUE)

RDA.corr_precip <- RDA.freq_env.sum %>% filter(precip_p.val < .05)
RDA.corr_elev <- RDA.freq_env.sum %>% filter(elev_p.val < .05)
RDA.corr_srad <- RDA.freq_env.sum %>% filter(srad_p.val < .05)

