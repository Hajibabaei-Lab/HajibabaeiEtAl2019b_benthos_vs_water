# Teresita M. Porter, Sept. 9, 2019

library(vegan)
library(reshape2)
library(gridExtra)
library(grid)
library(ggplot2)
library(plyr)
library(data.table)
library(stringr)

# Read infile
A<-read.table(file="matrix.csv", head=TRUE, sep=",")

# Select phylum Arthropoda only
B<-A[A$Phylum=="Arthropoda",]

# Narrow down to EPT + O (odonata), do subsetting after rarefaction to retain more data
C<-B[B$Order == "Ephemeroptera" |
      B$Order == "Plecoptera_Insecta" |
      B$Order == "Trichoptera" |
      B$Order == "Odonata",]

C2<-data.frame(C, do.call(rbind, str_split(C$SampleName,"_")))
names(C2)[32:38]<-c("month","year","version","substrate","site","marker","siterep")

# Split siterep into their own columns
C2$site<-substr(C2$siterep, 1, 2)
C2$rep<-substring(C2$siterep, 3)
C2$siterep<-NULL

# Combine substrate+site into their own column
C2$sample<-paste(C2$substrate, C2$site, sep="_")

# Pivot to make matrix for vegan
D<-dcast(C2, Marker_GlobalESV ~ sample, value.var="ESVsize", fun.aggregate = sum)

# Move marker_OTU to row names
rownames(D)<-D$Marker_GlobalESV
D<-D[,-1]

# Transpose to get sites in rows, ESVs in columns
Dt<-t(D)

# Remove columns with only zeros
notnull<-Dt[,colSums(Dt) !=0]

# Remove rows with only zeros
notnull2<-notnull[rowSums(notnull) !=0,]

# Calculate 15th percentile for rrarefy function
percentile<-quantile(rowSums(notnull2), prob=0.15)

# Set random seed for rarefaction
set.seed(12345)

# Rarefy the dataset down to the 15th percentile
df<-rrarefy(notnull2,sample=percentile)

# Convert to presence-absence matrix
df[df>0] <-1

# Do 2-dimensional NMDS
nmds3<-metaMDS(df, k=3, trymax=1000)
# stress = 0.02193911 
# Linear fit, R2 = 0.996

# Create grouping matrix for samples by grabbing row names from above matrix
sample_df<-data.frame(row.names(df))

# Rename the column
names(sample_df)<-"sample"

# Copy column to row names
row.names(sample_df)<-sample_df$sample

# Split first column into their own fields
sample_df[,2:3]<-do.call('rbind', strsplit(as.character(sample_df$sample),'_',fixed=TRUE))

# Remove first column
sample_df<-sample_df[,-1]

# Rename columns
names(sample_df)<-c("substrate","site")

# Grab sites/species scores from NMDS output
site.sc <- data.frame(scores(nmds3, display = "sites"))
#species.sc <- data.frame(scores(nmds2, display = "species"))

# Put it all in one df for ggplot
merged <- merge(site.sc,sample_df,by="row.names")
colnames(merged)[colnames(merged)=="Row.names"] <- "sample"

# Add river
merged$river<-c(rep("Athabasca", 4),rep("Peace",4),rep("Athabasca", 4),rep("Peace",3))

# Create factors
merged$substrate<-factor(merged$substrate, levels=c("B","W"), labels=c("Benthos","Water"))

# Create metadata from rownames 'sample'
# Only keep the samples that are in df
ENV<-merged[,c(1,5:7)]

# Create factors
ENV$substrate<-factor(ENV$substrate, levels=c("Benthos","Water"))
ENV$river<-factor(ENV$river, levels=c("Athabasca", "Peace"))

# Assess dispersion (variance) using ANOVA
# Create distance matrix based on P-A data using jaccard dissimilarity
sor<-vegdist(df, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
#site
bd_substrate<-betadisper(sor, factor(ENV[,2]))
bd_site<-betadisper(sor, factor(ENV[,3]))
bd_river<-betadisper(sor, factor(ENV[,4]))

# check for heterogeneity of beta dispersions within groups
anova(bd_substrate) # n/s
anova(bd_river) # n/s

pdf("BetaDispersion_EPTO.pdf")
par(mfrow=c(2,2))
boxplot(bd_substrate)
boxplot(bd_river)
dev.off()

# Shephards curve and goodness of fit calcs
# Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
pdf("stressplot_EPTO.pdf")
stressplot(nmds3)
gof <-goodness(nmds3)
gof
plot(nmds3, display = "sites", type="n")
points(nmds3, display="sites",cex=2*gof/mean(gof))
dev.off()
#Linear fit, R2 = 1

# Use ADONIS to test significance of groupings (esp. sites, layers, don't expect any diff for expt)

# 1. Test the interaction among all groups first to see how to proceed (use strata so randomizations occur within each site)
adonis(sor~substrate*river, data=ENV, permutations=999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# substrate        1    0.7984 0.79836  2.0873 0.13060  0.002 **
#   river            1    0.6635 0.66354  1.7348 0.10854  0.021 * 
#   substrate:river  1    0.4438 0.44379  1.1603 0.07260  0.282   
# Residuals       11    4.2074 0.38249         0.68826          
# Total           14    6.1131                 1.00000     
# no sig interaction between substrate and river

adonis(sor~substrate, data=ENV, permutations=999)
# Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)  
# substrate  1    0.7984 0.79836  1.9528 0.1306  0.011 *
#   Residuals 13    5.3148 0.40883         0.8694         
# Total     14    6.1131                 1.0000  

adonis(sor~river, data=ENV, permutations=999, strata=ENV$substrate)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# river      1    0.6549 0.65495  1.5599 0.10714  0.031 *
#   Residuals 13    5.4582 0.41986         0.89286         
# Total     14    6.1131                 1.00000   
