library(adegenet)
library(tidyverse)
library(vegan)

## load VCF file
## colias.4x.merged_gatk.rm.relate.SNP.filtered_gatkVQSR2.PASS.8miss.recode.vcf.gl_impute4.1.raw

raw_vcf <- "/scratch/alpine/ldurkee@colostate.edu/colias2023/RDA/colias.4x.merged_gatk.rm.relate.SNP.filtered_gatkVQSR2.PASS.8miss.recode.vcf.gl_impute4.1.rm.WY.raw"

# read the .raw file into R
colias_vcf <- read.PLINK(raw_vcf)

## load bioclim variables and elevation
colias_clim.final <- read.csv("/scratch/alpine/ldurkee@colostate.edu/colias2023/RDA/colias_bioclim-elev.csv")

paste("TRUE means VCF and env orders match")
identical(colias_vcf$ind.names, as.character(colias_clim.final[,1])) 

## subset of variables
colias_env <- colias_clim.final %>% dplyr::select(c(elevation, precip, srad))

## run the RDA
colias_rda <- rda(colias_vcf ~ ., data=colias_env, scale=T)
colias_rda

# save as RDS after running
# saveRDS(colias_rda,file="colias_RDA-env.rds")
# colias_rda <-readRDS("Data/colias_RDA-env.rds")


#### RDA model analyses #### 

## regression coefficients
coef(colias_rda)

paste("proportion of variance explained by each axis")
rda_sum <- summary(colias_rda)
rda_sum$concont$importance[2,]

paste("R-squared")
RsquareAdj(colias_rda)


# eigenvector analysis
rda_sum.axis <- summary(eigenvals(colias_rda, model = "constrained"))
# write.csv(rda.sum, "July24-eigen-sum-env.csv")
rda_sum.axis

# visualize the axes
screeplot(colias_rda)

# variance inflation factor - tells if variables are correlated
# confirm values are < 10
paste("VIF")
vif <- vif.cca(colias_rda)
# write.csv(vif, "July24-vif-MEM.csv")
vif 

# significance
signif.full <- anova.cca(colias_rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full

#### RDA plots ####
# simple plot, axes 1 & 2
# pdf("colias-RDA-simple-full.pdf")
# plot(colias_rda, scaling=3)
# dev.off()

# more informative plots - color coded by site
colias_env$site <- colias_clim.final$site
unique(colias_env$site)

levels(colias_env$site) <- as.factor(c("Babbit Gulch","Carpenter Ranch",
                                       "Dumont Lake","Fountain Valley",
                                       "Harbison Meadows","High Trails Ranch",
                                       "Kebler Pass","D Loukonen Farm",
                                       "Little Laramie","Schreiner Farm",
                                       "North Fork Trail","Orchard Creek Ranch",
                                       "Spring Canyon Park","Soapstone Prairie"))

eco <- colias_env$site

# bg <- c("#750D37","#750D37", "#A31621", "#A31621", "#ff9896","#ff9896",
# "#031A6B","#031A6B", "#306BAC","#306BAC",  "#6F9CEB", "#6F9CEB",
# "#aec7e8","#aec7e8","#918EF4", "#918EF4")

# sh <- rep(c(21, 24), 7)

# site_colors <- bg[match(eco, levels(eco))]
# site_shapes <- sh[match(eco, levels(eco))]

bg <- c("#FF3030","#750D37",
        "#750D37","#B2DFEE", 
        "#AA8EF4", "#B2DFEE", 
        "#ff9896", "#AA8EF4",
        "#031A6B", "#ff9896",
        "#1E90FF", "#FF3030",
        "#1E90FF", "#031A6B")

sh <- c(24, 21, 
        24, 21,
        24, 24,
        24, 21,
        24, 21,
        24, 21,
        21, 21)

site_colors <- bg[match(eco, levels(eco))]
site_shapes <- sh[match(eco, levels(eco))]

pdf("July24-RDA-env.only_1-2-v2.pdf", width=12, height=10)

## plot for axes 1 & 2
plot(colias_rda, type="n", scaling=3, cex.axis=1.2, cex.lab=1.2)
points(colias_rda, display="species", pch=4, cex=.5, col="gray32", scaling=3)  # the SNPs
points(colias_rda, display="sites", pch=site_shapes, cex=2, scaling=3, bg=site_colors)  # the butterflies
text(colias_rda, scaling=3, display="bp", col="darkorchid", cex=1.5)  # the predictors
legend("topright", legend=levels(eco), bty="n", col="gray32", pch=sh, cex=1, pt.bg=bg)
dev.off()

## axis 1-3
pdf("July24-RDA-env.only_1-3-v2.pdf", width=16, height=12)

plot(colias_rda, type="n", scaling=3, cex.axis=1.2, cex.lab=1.2, choices=c(1,3))
points(colias_rda, display="species", pch=4, cex=.5, col="gray32", scaling=3, choices=c(1,3))  # the SNPs
points(colias_rda, display="sites", pch=site_shapes, cex=2, scaling=3, bg=site_colors, choices=c(1,3))  # the butterflies
text(colias_rda, scaling=3, display="bp", col="darkorchid", cex=1.5, choices=c(1,3))  # the predictors
legend("topright", legend=levels(eco), bty="n", col="gray32", pch=sh, cex=1, pt.bg=bg)

dev.off()

#### Looking for candidate loci ####
load.rda <- scores(colias_rda, choices=c(1:3), display="species")

# we are interested in SNPs at the tails of the distributions
# this function will be able to detect outlier SNPs

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# candidate loci
cand1 <- outliers(load.rda[,1],3) 
paste("no. candidate loci on axis 1 = ", length(cand1))

cand2 <- outliers(load.rda[,2],3) 
paste("no. candidate loci on axis 2 = ", length(cand2))

cand3 <- outliers(load.rda[,3],3) 
paste("no. candidate loci on axis 3 = ", length(cand3))

# total number of candidate loci identified
ncand <- length(cand1) + length(cand2) + length(cand3)
ncand

# now, organize the results in one data frame with axis, SNP name, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

## see which environmental vars they are most correlated with

for (i in 1:length(cand_no.dup$snp)) {
  bar <- cand_no.dup[i,]
  cand_no.dup[i,7] <- names(which.max(abs(bar[3:6]))) # gives the variable
  cand_no.dup[i,8] <- max(abs(bar[3:6]))              # gives the correlation
}

colnames(cand_no.dup)[7] <- "predictor"
colnames(cand_no.dup)[8] <- "correlation"

#write.csv(cand_no.dup, "cand-loci-env.only-thresh3.csv")

#### exploring significant genes ####

snps <- read.csv("Data/SNPs-env.only-pred.cor-thresh3-full.csv")
#snps$loc <- sub("_[ATGC]$", "", snps$snp)

#snps.sum <- snps %>% group_by(snp) %>%
#summarize(count=n())

#### Manhattan plot ####
library(ggplot2)
library(dplyr)
library(stringr)

# Step 1: Make chromosome a factor in proper order
# Replace the actual chromosome names in your data with 1:30 and Z
# Assuming snps$chrom has 31 unique names
chroms_sorted <- c(as.character(1:30), "Z")  # desired order
snps$chrom <- factor(snps$chrom, levels = chroms_sorted)

# Step 2: Compute cumulative positions
chrom_info <- snps %>%
  group_by(chrom) %>%
  summarise(chr_len = max(pos)) %>%
  mutate(tot = lag(cumsum(chr_len), default = 0))  # total length of previous chromosomes

snps <- snps %>%
  left_join(chrom_info %>% select(chrom, tot), by = "chrom") %>%
  arrange(chrom, pos) %>%
  mutate(bp_cum = pos + tot)

# Step 3: Get chromosome midpoints for x-axis
axis_snps <- snps %>%
  group_by(chrom) %>%
  summarise(center = (min(bp_cum) + max(bp_cum)) / 2)

#### Genes of interest ####
snps_to_plot <- c(
  "NC_059568.1.4450341", # TH
  "NC_059542.1.4495778", # yellow-like
  "NC_059541.1.7361965", # hsp70
  "NC_059555.1.678453"   # cytochrome P450 (CYP6B6)
)

snps_to_plot <- data.frame(
  loc = c(
    "NC_059568.1.4450341", 
    "NC_059542.1.4495778", 
    "NC_059541.1.7361965", 
    "NC_059555.1.678453"
  ),
  gene = c("TH", "yellow-like", "HSP70", "CYP6B6")
)

# Merge gene names into your snps table
snps_highlight <- merge(snps_to_plot, snps, by="loc")

# Step 4: Plot
ggplot(snps, aes(x = bp_cum, y = abs(loading), color = chrom)) +
  geom_point(alpha = 0.4, size = 1.1) +
  geom_point(
    data = snps_highlight,
    aes(x = bp_cum, y = abs(loading)),
    color = "coral",  # highlight color
    size = 3
  ) +
  # SNP labels
  geom_text(
    data = snps_highlight,
    aes(x = bp_cum, y = abs(loading), label = gene),
    size = 5,
    fontface = "bold",
    color = "coral",
    vjust = c(1.6, -.5, -1, -1),
    hjust=c(1,-.1,1,1)
  ) +
  scale_color_manual(values = rep(c("grey", "black"), length.out = 31)) +
  scale_x_continuous(label = axis_snps$chrom, breaks = axis_snps$center) +
  labs(x = "Chromosome", y = expression(abs(" RDA loading "))) +
  theme_minimal() +
  theme_font1()+
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(hjust = 1, size = 7) # smaller font & angled
  )
