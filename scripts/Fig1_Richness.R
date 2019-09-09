# Teresita M. Porter, Sept. 9, 2019

library(vegan)
library(reshape2)
library(gridExtra)
library(grid)
library(ggplot2)
library(plyr)
library(data.table)
library("car")
library(stringr)
library("ggpubr")

# Read infile
A<-read.table(file="matrix.csv", head=TRUE, sep=",")

# Select phylum Arthropoda only
B<-A[A$Phylum=="Arthropoda",]

# Split SampleName into their own columns using pkg "stringr"
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

# Convert to presence-absence
df[df>0] <- 1

# Convert to df
df<-data.frame(df)  

# Get total ESVs per sample
df$sums<-rowSums(df)

# Move rownames to first column
df2<-data.frame(df)
setDT(df2, keep.rownames = TRUE)[]

# Get separate substrate and siterep cols
setDT(df2)[, paste0("S", 1:2) := tstrsplit(rn, "_")]
colnames(df2)[colnames(df2)=="S1"] <- "substrate"
colnames(df2)[colnames(df2)=="S2"] <- "siterep"

# Get separate site and rep cols
df2$site <- str_sub(df2$siterep, 1,2)
df2$rep <- str_sub(df2$siterep, -1)

# Create factors
df2$substrate<-factor(df2$substrate, levels=c("B","W"), labels=c("Benthos","Water"))

#mggplot boxplot
p1<-ggplot(df2) +
#  geom_boxplot(aes(x=df2$site, y=df2$sums, fill=substrate)) +
  geom_point(aes(x=df2$site, y=df2$sums, color=substrate)) +
  labs(x="Sites", y="ESV Richness") +
  facet_wrap(~substrate) +
  scale_color_manual(values=c("#4DAF4A", "#377EB8")) +
  theme(legend.title=element_blank()) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank())

ggsave("F1_ESVrichness.pdf",p1)

# test for normality pkg "ggpubr"
ggdensity(df2$sums, 
          main = "Density plot",
          xlab = "Sample ESV richness")
# not normal

ggqqplot(df2$sums)
# mostly normal

qqPlot(df2$sums)
# mostly normal

#Shapiro-Wilk test of normality
shapiro.test(df2$sums)
# data:  df2$sums
# W = 0.94474, p-value = 0.02477
# sig diff than normal

#paired samples Wilcoxon test (Wilcoxon signed-rank test)
wilcox.test(df2$sums[df2$substrate=="Benthos"], df2$sums[df2$substrate=="Water"], 
            paired = TRUE, alternative = "greater")
# p-value = 3.217e-05, Benthos richness greater than Water richness
