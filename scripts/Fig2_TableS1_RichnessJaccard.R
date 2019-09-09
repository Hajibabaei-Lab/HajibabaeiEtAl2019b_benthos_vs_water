# Teresita M. Porter, Sept. 9, 2019

library(vegan)
library(reshape2)
library(gridExtra)
library(grid)
library(ggplot2)
library(plyr)
library(data.table)
library(stringr) # str_split
library(ggpubr)
library(philentropy) # distance method
library(gdata) # upper triangle
library(Ternary) # ternary plot

# Read infile
A<-read.table(file="matrix.csv", head=TRUE, sep=",")

# Select phylum Arthropoda only
B<-A[A$Phylum=="Arthropoda",]

# Split SampleName into their own columns from pkg 'stringr'
B2<-data.frame(B, do.call(rbind, str_split(B$SampleName,"_")))
names(B2)[32:38]<-c("month","year","version","substrate","site","marker","siterep")

# Get separate substrate and siterep cols
B2$site<-substr(B2$siterep, 1,2)

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

# Convert to df
df3<-as.data.frame(df2)

# Move rn to row names
row.names(df3)<-df3$rn
df3$rn<-NULL

# remove last 5 cols of metadata
range=(length(colnames(df3))-4):length(colnames(df3))
df3<-df3[,-c(range)]

# Subset by benthos/water collection method
df4_b<-df3[grepl("^B_",rownames(df3)),]
df4_w<-df3[grepl("^W_",rownames(df3)),]

# Get colsums
sum_b<-colSums(df4_b)
sum_w<-colSums(df4_w)

# Convert to presence-absence
sum_b[sum_b>0]<-1
sum_w[sum_w>0]<-1

# Append rowsums to df
df5_b<-rbind(df4_b, sum_b)
df5_w<-rbind(df4_w, sum_w)

# rename last rown
rownames(df5_b)[length(rownames(df5_b))]<-"Benthos_PA"
rownames(df5_w)[length(rownames(df5_w))]<-"Water_PA"

# Get total richness stats
vector_b<-as.numeric(df5_b[length(rownames(df5_b)),])
richness_b<-sum(vector_b)
# 1588

vector_w<-as.numeric(df5_w[length(rownames(df5_w)),])
richness_w<-sum(vector_w)
# 658

ratio_richness<-richness_b/richness_w
# 2.413374

######################################################################################
# Custom Jaccard function
# https://stats.stackexchange.com/questions/176613/jaccard-similarity-in-r/303078
# m=1 for long format 
# m=2 for wide format

jaccard <- function(df, margin) {
  if (margin == 1 | margin == 2) {
    M_00 <- apply(df, margin, sum) == 0
    M_11 <- apply(df, margin, sum) == 2
    if (margin == 1) {
      df <- df[!M_00, ]
      JSim <- sum(M_11) / nrow(df)
    } else {
      df <- df[, !M_00]
      JSim <- sum(M_11) / length(df)
    }
    JDist <- 1 - JSim
    return(c(JSim = JSim, JDist = JDist))
  } else break
}
######################################################################################


# Combine water and benthos presence-absence rows
combo<-rbind(df5_b["Benthos_PA",],df5_w["Water_PA",])

# Get jaccard index/dist
jaccard(combo,2)
# JSim     JDist 
# 0.1430025 0.8569975 

# Create ternary plot

cairo_pdf("F2_Ternary.pdf",family="Arial Unicode MS")
TernaryPlot(alab="% ESVs Unique to Benthos \u2192", 
            blab="% ESVs Unique to Water \u2192", 
            clab="% ESVs Shared \u2192",
            point='up', lab.cex=0.8, grid.minor.lines = 0,
            grid.lty='solid', col=rgb(0.9, 0.9, 0.9), grid.col='white', 
            axis.col=rgb(0.6, 0.6, 0.6), ticks.col=rgb(0.6, 0.6, 0.6),
            padding=0.08)

C<-dcast(B2, substrate+siterep~Marker_GlobalESV, value.var="ESVsize", fun.aggregate = sum)

# initialize list (length equal to number of sites)
l <- list()

for (i in unique(C$site)) { # loop through sites

    b<-colnames(C)[as.numeric(C[C$substrate=="B" & C$site==i,])>0]
    w<-colnames(C)[as.numeric(C[C$substrate=="W" & C$site==i,])>0]
    
    # remove substrate and site, first two elements
    b<-b[-c(1:2)]
    w<-w[-c(1:2)]
    
    benthos=setdiff(b,w)
    water=setdiff(w,b)
    shared=intersect(b,w)
    
    l[[i]] <- c(length(benthos), length(water), length(shared))
    
}

AddToTernary(text, l, names(l), cex=0.8, font=2)

dev.off()

