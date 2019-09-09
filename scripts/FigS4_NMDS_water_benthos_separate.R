# Teresita M. Porter, Sept. 9, 2019

library(vegan)
library(reshape2)
library(ggplot2)
library(data.table)
library(goeveg) # scree
library(stringr) #str_split

# Read infile
A<-read.table(file="matrix.csv", head=TRUE, sep=",")

# Select phylum Arthropoda only
B<-A[A$Phylum=="Arthropoda",]

# Split SampleName into their own columns
B2<-data.frame(B, do.call(rbind, str_split(B$SampleName,"_")))
names(B2)[32:38]<-c("month","year","version","substrate","site","marker","siterep")

# Combine substrate+siterep into their own column
B2$sample<-paste(B2$substrate,B2$siterep,sep="_")

# Pivot to make matrix for vegan
C<-dcast(B2, Marker_GlobalESV ~ sample, value.var="ESVsize", fun.aggregate = sum)

# Move marker_OTU to row names
rownames(C)<-C$sample
C<-C[,-1]

# Transpose to get sites in rows, ESVs in columns
Ct<-t(C)

# Remove columns with only zeros
notnull<-Ct[,colSums(Ct) !=0]

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

# Create separate plots for benthos and water for supplement FigS4
df_b<-df[grepl("^B_", rownames(df)),]
df_w<-df[grepl("^W_", rownames(df)),]

# # Scree plots to determine number of dimensions to use for NMDS, use k=3
# pdf("Scree_b.pdf")
# # check dims
# dimcheckMDS(df_b)
# dev.off()
# 
# pdf("Scree_w.pdf")
# # check dims
# dimcheckMDS(df_w)
# dev.off()

# Do 2 dimensional NMDS, no environment file for now 
nmds3_b<-metaMDS(df_b, k=3, trymax=100)
# stress = 0.08290066
# Linear fit, R2 = 0.952
nmds3_w<-metaMDS(df_w, k=3, trymax=100)
# stress = 0.09278044
# Linear fit, R2 = 0.949

# Create grouping matrix for samples by grabbing row names from above matrix
sample_df_b<-data.frame(row.names(df_b))
sample_df_w<-data.frame(row.names(df_w))

# Rename the column
names(sample_df_b)<-"sample"
names(sample_df_w)<-"sample"

# Copy column to row names
row.names(sample_df_b)<-sample_df_b$sample
row.names(sample_df_w)<-sample_df_w$sample

# Split first column into their own fields
sample_df_b[,2:3]<-do.call('rbind', strsplit(as.character(sample_df_b$sample),'_',fixed=TRUE))
sample_df_w[,2:3]<-do.call('rbind', strsplit(as.character(sample_df_w$sample),'_',fixed=TRUE))

# Remove first column
sample_df_b<-sample_df_b[,-1]
sample_df_w<-sample_df_w[,-1]

# Rename columns
names(sample_df_b)<-c("substrate","siterep")
names(sample_df_w)<-c("substrate","siterep")

# Grab sites/species scores from NMDS output
site.sc.b <- data.frame(scores(nmds3_b, display = "sites"))
site.sc.w <- data.frame(scores(nmds3_w, display = "sites"))

# Put it all in one df for ggplot
merged_b <- merge(site.sc.b,sample_df_b,by="row.names")
colnames(merged_b)[colnames(merged_b)=="Row.names"] <- "sample"
merged_w <- merge(site.sc.w,sample_df_w,by="row.names")
colnames(merged_w)[colnames(merged_w)=="Row.names"] <- "sample"

# Add river
merged_b$river<-c(rep("Athabasca", 12),rep("Peace",12))
merged_w$river<-c(rep("Athabasca", 12),rep("Peace",12))

# Compile coord for convex hulls
chulls12_b <- ddply(merged_b, .(river), function(merged_b) merged_b[chull(merged_b$NMDS1, merged_b$NMDS2), ])
chulls12_w <- ddply(merged_w, .(river), function(merged_w) merged_w[chull(merged_w$NMDS1, merged_w$NMDS2), ])

# Create scatter plot for BENTHOS
p.b<-ggplot(data=merged_b, aes(x=NMDS1, y=NMDS2)) + 
  ggtitle("Benthos") +
  geom_polygon(data=chulls12_b, aes(x=NMDS1, y=NMDS2, fill=river), alpha=0.5) +
  geom_text(data=merged_b, aes(label=siterep, color=river)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title = element_text(size=12),
    axis.text = element_text(size=12),
    legend.title = element_blank(),
    legend.text = element_text(size=12)) +
  guides(color=FALSE)

# Create scatter plot for WATER
p.w<-ggplot(data=merged_w, aes(x=NMDS1, y=NMDS2)) + 
  ggtitle("Water") +
  geom_polygon(data=chulls12_w, aes(x=NMDS1, y=NMDS2, fill=river), alpha=0.5) +
  geom_text(data=merged_w, aes(label=siterep, color=river)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title = element_text(size=12),
    axis.text = element_text(size=12),
    legend.title = element_blank(),
    legend.text = element_text(size=12)) +
  guides(color=FALSE)

g<-grid.arrange(p.b,p.w,nrow=2)
ggsave("FS4_nmds3_water_benthos_sep.pdf",g, height=8.5, width=8.5, units="in")

# Create metadata from rownames 'sample'
# Only keep the samples that are in df
ENV_b<-merged_b[,c(1,5:7)]
ENV_w<-merged_w[,c(1,5:7)]

# Create factors
ENV_b$river<-factor(ENV_b$river, levels=c("Athabasca","Peace"))
ENV_w$river<-factor(ENV_w$river, levels=c("Athabasca","Peace"))

# Assess dispersion (variance) using ANOVA
# Create distance matrix based on P-A data using jaccard dissimilarity
sor_b<-vegdist(df_b, "bray", binary=TRUE)
sor_w<-vegdist(df_w, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
#site
bd_river_b<-betadisper(sor_b, factor(ENV_b[,4]))
bd_river_w<-betadisper(sor_w, factor(ENV_w[,4]))

# check for heterogeneity of beta dispersions within groups
anova(bd_river_b) # n/s
anova(bd_river_w) # 0.0007742 ***

pdf("BetaDispersion_b_w.pdf")
par(mfrow=c(2,2))
boxplot(bd_river_b, main="Benthos")
boxplot(bd_river_w, main="Water")
dev.off()

# Shephards curve and goodness of fit calcs
# Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
pdf("stressplot_b_w.pdf")
stressplot(nmds3_b, main="Benthos")
gof <-goodness(nmds3_b)
gof
plot(nmds3_b, display = "sites", type="n", main="Benthos")
points(nmds3_b, display="sites",cex=2*gof/mean(gof))

stressplot(nmds3_w, main="Water")
gof <-goodness(nmds3_w)
gof
plot(nmds3_w, display = "sites", type="n", main="Water")
points(nmds3_w, display="sites",cex=2*gof/mean(gof))
dev.off()

# Use ADONIS to test significance of groupings (esp. sites, layers, don't expect any diff for expt)
adonis(sor_b~river, data=ENV_b, permutations=999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# river      1    1.4262 1.42618  4.9992 0.18516  0.001 ***
#   Residuals 22    6.2762 0.28528         0.81484           
# Total     23    7.7024                 1.00000  
# River accounts for 18.5% of variation pval=0.001 from benthos samples

adonis(sor_w~river, data=ENV_w, permutations=999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# river      1    1.2717 1.27168  3.4681 0.13617  0.001 ***
#   Residuals 22    8.0669 0.36668         0.86383           
# Total     23    9.3386                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# River accounts for 15% of variation pval=0.001 from water samples