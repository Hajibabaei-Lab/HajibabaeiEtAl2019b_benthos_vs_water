# Teresita M. Porter, Sept. 9, 2019

library(vegan)
library(reshape2)
library(ggplot2)
library(data.table)
library(goeveg) # scree
library(RColorBrewer)
library(stringr) # str_split

# Read infile
A<-read.table(file="matrix.csv", head=TRUE, sep=",")

# Select phylum Arthropoda only
B<-A[A$Phylum=="Arthropoda",]

# Split SampleName into their own columns from pkg 'stringr'
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

# Scree plots to determine number of dimensions to use for NMDS
#pdf("Scree.pdf")
# check dims
#dimcheckMDS(df)
#dev.off()

# Do 2 dimensional NMDS, no environment file for now 
nmds3<-metaMDS(df, k=3, trymax=100)
# stress = 0.1221286
# Linear fit, R2 = 0.908

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
names(sample_df)<-c("substrate","siterep")

# Grab sites/species scores from NMDS output
site.sc <- data.frame(scores(nmds3, display = "sites"))
#species.sc <- data.frame(scores(nmds2, display = "species"))

# Put it all in one df for ggplot
merged <- merge(site.sc,sample_df,by="row.names")
colnames(merged)[colnames(merged)=="Row.names"] <- "sample"

# Add river

merged$river<-c(rep("Athabasca", 12),rep("Peace",12),rep("Athabasca", 12),rep("Peace",12))

# Create factors
merged$substrate<-factor(merged$substrate, levels=c("B","W"), labels=c("Benthos","Water"))

# Compile coord for convex hulls
chulls12 <- ddply(merged, .(substrate), function(merged) merged[chull(merged$NMDS1, merged$NMDS2), ])

# Create scatter plot highlighting collection method
p1<-ggplot(data=merged, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls12, aes(x=NMDS1, y=NMDS2, fill=substrate), alpha=0.5) +
  geom_text(data=merged, aes(label=siterep, color=substrate)) +
  scale_fill_manual(values=c("#4DAF4A","#377EB8")) +
  scale_color_manual(values=c("#4DAF4A","#377EB8"), guide=FALSE) +
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
    legend.text = element_text(size=12))

# Compile coord for convex hulls
chulls12b <- ddply(merged, .(river), function(merged) merged[chull(merged$NMDS1, merged$NMDS2), ])

brewer.pal(n = 8, name = "Set1")

# Create scatter plot highlighting watershed
p2<-ggplot(data=merged, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls12b, aes(x=NMDS1, y=NMDS2, fill=river), alpha=0.5) +
  geom_text(data=merged, aes(label=siterep, color=river)) +
  scale_fill_manual(values=c("#E41A1C","#984EA3")) +
  scale_color_manual(values=c("#E41A1C","#984EA3"), guide=FALSE) +
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
    legend.text = element_text(size=12))

g<-grid.arrange(p1,p2,nrow=2)
ggsave("F4_nmds.pdf",g, height=8.5, width=8.5, units="in")

# Create metadata from rownames 'sample'
# Only keep the samples that are in df
ENV<-merged[,c(1,5:7)]

# Create factors
ENV$substrate<-factor(ENV$substrate, levels=c("Benthos","Water"))

# Assess dispersion (variance) using ANOVA
# Create distance matrix based on P-A data using jaccard dissimilarity
sor<-vegdist(df, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
#site
bd_substrate<-betadisper(sor, factor(ENV[,2]))
bd_siterep<-betadisper(sor, factor(ENV[,3]))
bd_river<-betadisper(sor, factor(ENV[,4]))

# check for heterogeneity of beta dispersions within groups
anova(bd_substrate) # 0.001588 ** heterogeneous dispersion of beta div but we have balanced design
anova(bd_siterep) # < 2.2e-16 *** as above
anova(bd_river) # 0.0001255 *** as above
pdf("BetaDispersion.pdf")
par(mfrow=c(2,2))
boxplot(bd_substrate)
boxplot(bd_siterep, las=2, cex.axis=0.75) # rotate x-axis labels 90 degrees
boxplot(bd_river)
dev.off()

# Shephards curve and goodness of fit calcs
# Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
pdf("stressplot.pdf")
stressplot(nmds3)
gof <-goodness(nmds3)
gof
plot(nmds3, display = "sites", type="n")
points(nmds3, display="sites",cex=2*gof/mean(gof))
dev.off()
#Linear fit, R2 = 0.908

# Use ADONIS to test significance of groupings (esp. sites, layers, don't expect any diff for expt)

# 1. Test the interaction among all groups first to see how to proceed (use strata so randomizations occur within each site)
adonis(sor~substrate*river*siterep, data=ENV, permutations=999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# substrate          1    1.7200       2       0 0.09168      1
# river              1    1.8108       2       0 0.09652      1
# siterep           22    8.0755       0       0 0.43044      1
# substrate:river    1    0.8870       1       0 0.04728      1
# substrate:siterep 22    6.2676       0       0 0.33408      1
# Residuals          0    0.0000    -Inf         0.00000       
# Total             47   18.7610                 1.00000       
# sig interaction between substrate & river
# no sig interactin between substrate & siterep

adonis(sor~substrate, data=ENV, permutations=999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# substrate  1     1.720 1.72004  4.6431 0.09168  0.001 ***
#   Residuals 46    17.041 0.37045         0.90832           
# Total     47    18.761                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(sor~river, data=ENV, permutations=999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# river      1    1.8108 1.81083  4.9143 0.09652  0.001 ***
#   Residuals 46   16.9501 0.36848         0.90348           
# Total     47   18.7610                 1.00000  

adonis(sor~siterep, data=ENV, permutations=999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# siterep   23    9.8863 0.42984  1.1624 0.52696   0.02 *
#   Residuals 24    8.8746 0.36978         0.47304         
# Total     47   18.7610                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1