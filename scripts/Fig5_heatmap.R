# Teresita M. Porter, Sept. 9, 2019

library(stringr)
library(reshape2)
library(vegan)
library(ggplot2)
library(data.table)

#####################################################################

# Read in sample x taxonomy table
A<-read.csv(file="matrix.csv", head=TRUE)

# Filter table for Arthropoda only, use for ESVs, order, class
B<-A[A$Phylum=="Arthropoda",]

# For species use sBP>=0.70 (95% correct)
C<-B[B$fBP>=0.20,]

#######################################################
# Create dataframe for ESV rarefaction
## DO NOT pool data across reps
######################################################

# Split SampleName into their own columns from pkg 'stringr'
C2<-data.frame(C, do.call(rbind, str_split(C$SampleName,"_")))
names(C2)[32:38]<-c("month","year","version","substrate","site","marker","siterep")
C2$site<-substr(C2$siterep, 1, 2)
C2$rep<-substring(C2$siterep, 3)
C2$siterep<-NULL

# Combine substrate+siterep into their own column
C2$sample<-paste(C2$substrate, C2$site,sep="_")

# Pivot to make matrix for vegan
D<-dcast(C2, substrate+site+Marker_GlobalESV+Order+Family ~ sample, value.var="ESVsize", fun.aggregate = sum)

# merge ESV+Order+Family into single column
D$sample<-paste(D$substrate, D$site, D$Marker_GlobalESV, D$Order, D$Family, sep="-")

# delete duplicate columns
D$substrate<-NULL
D$site<-NULL
D$Marker_GlobalESV<-NULL
D$Order<-NULL
D$Family<-NULL

#move SampleName to rownames then delete
rownames(D)<-D$sample
D$sample<-NULL

#remove columns with only zeros
notnull<-D[,colSums(D) !=0]

#remove rows with only zeros & edit rownames
notnull2<-notnull[rowSums(notnull) !=0,]

#calculate 15th percentile for rrarefy function
percentile<-quantile(rowSums(notnull2), prob=0.15)

# Rarefy dataset
set.seed(1234)
df<-rrarefy(notnull2, sample=percentile)

# Convert to presence-absence matrix
df[df>0]<-1

# Move rownames to columns
df2<-data.frame(df)
setDT(df2, keep.rownames = TRUE)[]

# Split rn into own columns
df3<-data.frame(df2, do.call(rbind, str_split(df2$rn,"-")))
names(df3)[18:22]<-c("substrate","site","marker","order","family")
df3$rn<-NULL

# calculate rowSums
df3$sum <- rowSums(df3[,1:16])

# pivot to summarize
df4<-dcast(df3, order+family ~ substrate+site, value.var="sum", fun.aggregate = sum)

# create long form for ggplot
long<-melt(df4, id=c("order", "family"))

# split variable into their own fields
long<-data.frame(long, do.call(rbind, str_split(long$variable, "_")))

# rename columns
names(long)[4:6]<-c("ESVs","substrate","site")

# merge order and family into single field
long$taxon<-paste(long$order, long$family, sep="_")

# create factor
long$substrate<-factor(long$substrate, levels=c("B","W"), labels=c("Benthos","Water"))
long$taxon<-as.factor(long$taxon)

# breaks for log scale legend
my_breaks<-c(1, 100, 400)

# Sum total number of ESVs per substrate_site
TotESVs<-cbind(aggregate(ESVs~substrate+site, long, sum), taxon="Total ESVs")

# Count total number of EPTO families per substrate_site
# Create new ESV column with just presence absence data
long$ESVsPA<-long$ESVs
long$ESVsPA[long$ESVsPA>0] <-1
TotFam<-cbind(aggregate(ESVsPA~substrate+site, long, sum), taxon="Total Families")

h<-ggplot(data=long, aes(x=site,y=taxon)) +
  geom_tile(aes(fill=ESVs)) +
  labs(x="Sites") +
  facet_wrap(~substrate) +
  geom_text(aes(label=ESVs), 
            data=TotESVs,
            size=3,
            colour="red") + 
  geom_text(aes(label=ESVsPA), 
            data=TotFam,
            size=3,
            colour="red") + 
  scale_y_discrete(limits = c("","Total Families","","Total ESVs","", rev(levels(long$taxon))),
                  labels=c("","Total Families","","Total ESVs","", rev(levels(long$taxon))),
                  breaks=c("","Total Families","","Total ESVs","", rev(levels(long$taxon)))) +
  scale_fill_gradient(high="red",
                      low="yellow",
                      na.value="lightgrey",
                      name = "ESVs", 
                      trans = "log",
                      breaks = my_breaks, 
                      labels = my_breaks) +
  theme_bw(base_size=8) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(size=9),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8))

ggsave("F5_heatmap.pdf", height=10, width=7, units="in")
