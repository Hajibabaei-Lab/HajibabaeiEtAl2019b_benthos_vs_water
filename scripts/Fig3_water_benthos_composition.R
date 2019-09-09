# Teresita M. Porter, Sept. 9, 2019

library(vegan)
library(reshape2)
library(gridExtra)
library(ggplot2)
library(data.table)
library(stringr)
library(ggrepel)

# Read infile
A<-read.table(file="matrix.csv", head=TRUE, sep=",")

# Select phylum Arthropoda only
B<-A[A$Phylum=="Arthropoda",]

# Split SampleName into their own columns from pkg 'stringr'
B2<-data.frame(B, do.call(rbind, str_split(B$SampleName,"_")))
names(B2)[32:38]<-c("month","year","version","substrate","site","marker","siterep")

# Combine substrate+siterep into their own column
B2$ESV<-paste(B2$substrate,B2$Marker_GlobalESV,sep="_")

# Pivot to make matrix for vegan
C<-dcast(B2, Order+Genus+gBP ~ ESV, value.var="ESVsize", fun.aggregate = sum)

# Merge Order Genus and gBP into one col
C$sample<-paste(C$Order, C$Genus, C$gBP, sep=";")

# Move marker_OTU to row names
row.names(C)<-C$sample
C<-C[,-c(1:3,4908)]

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
df<-as.data.frame(df)  

# Calc colSums
benthos<-colSums(df[grepl("^B_",rownames(df)),])
water<-colSums(df[grepl("^W_",rownames(df)),])

# create df for ggplot scatter plot
merge<-data.frame(cbind(benthos,water))

# move rownames to first column
setDT(merge, keep.rownames = TRUE)[]

# create total col
merge$total<-benthos+water

# Sort by total, descending
merge2 <- merge[order(-total),] 

# Remove rows with total==0
merge2<-merge2[!(merge2$total==0),]

# rename rn
names(merge2)[1]<-"name"

# Split sample into their own columns
merge2[, c("Order", "Genus", "gBP") := tstrsplit(name, ";", fixed=TRUE)]

# Plot species only
merge3<-merge2[merge2$gBP>=0.30,]

p1<-ggplot(merge3, aes(x=water, y=benthos)) +
  geom_point(aes(colour=factor(Order)), position=position_jitter()) +
  labs(x="Water ESVs (log10)", y="Benthos ESVs (log10)", colour = "Order") +
  geom_abline(intercept = 0, slope = 1, linetype = 3) +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  geom_text_repel(aes(x = water, 
                      y = benthos, 
                      label = Genus),
                  data=merge3[merge3$benthos>=2 & merge3$water>=2,],
                  size=3.5) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=11),
        legend.position="bottom")

ggsave("F3_water_benthos.pdf",p1, width=8.5, height=5, units="in")
            