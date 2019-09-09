# Teresita M. Porter, Sept. 9, 2019

library(vegan)
library(ggplot2)
library(reshape2)
library(stringr)
library(data.table)

###################################################################

#read infile prepared by python script
A<-read.table(file="matrix.csv", head=TRUE, sep=",")

# Select phylum Arthropoda only
B<-A[A$Phylum=="Arthropoda",]

#Split up SampleName into their own fields
C<-data.frame(B, do.call(rbind, str_split(B$SampleName,"_")))

# Add new column names
names(C)[32:38]<-c("month","year","version","substrate","site","marker","siterep")

# Pivot to make matrices for vegan, ESVs in rows, sites in columns, reads in cells
C2<-dcast(C, substrate+siterep ~ Marker_GlobalESV, value.var="ESVsize", fun.aggregate = sum)

# Merge substrate + siterep, move to rownames
C2$sample<-paste(C2$substrate,C2$siterep,sep="_")
rownames(C2)<-C2$sample
D<-C2

#Get total number of reads per sample
D$sums<-rowSums(D[,3:4461])

# Create factors
D$substrate<-factor(D$substrate, levels=c("B","W"),
                        labels=c("Benthos","Water"))

#make new df for ggplot boxplot
p1<-ggplot(D) +
  geom_boxplot(aes(D$substrate,D$sums)) +
  labs(x="Collection Method", y="Reads per Sample") +
  scale_y_continuous(label=comma) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

# Prep pres-abs
E<-C2

# Remove non-numeric columns
E<-E[,-c(1:2,4462)]

# Convert to presence-absence
E[E>0] <-1

#Get total number of ESVs per sample
E$sums<-rowSums(E)

#Move rownames into first col and name cols
setDT(E, keep.rownames = TRUE)[]
E2<-data.frame(E, do.call(rbind, str_split(E$rn,"_")))
colnames(E2)[colnames(E2)=="X1"] <- "substrate"
colnames(E2)[colnames(E2)=="X2"] <- "siterep"

# Create factors
E2$substrate<- factor(E2$substrate, levels=c("B","W"), labels=c("Benthos","Water"))

#make new df for ggplot boxplot
p2<-ggplot(E2) +
  geom_boxplot(aes(E2$substrate,E2$sums)) +
  labs(x="Collection Method", y="ESVs per Sample") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

g<-grid.arrange(p1, p2, nrow=1)

ggsave("FS2_boxplots_reads_ESVs.pdf",g)
