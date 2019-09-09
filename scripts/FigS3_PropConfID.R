# Teresita M. Porter, Sept. 9, 2019

library(vegan)
library(reshape2)
library(ggplot2)

###################################################################

# Read infile prepared by python script
A<-read.table(file="matrix.csv", head=TRUE, sep=",")

# Select phylum Arthropoda only
B<-A[A$Phylum=="Arthropoda",]

# Calc total taxa and total taxa confidently id'd
species<-length(unique(B$Species))
species_good<-length(unique(B$Species[B$sBP>=0.70]))
genus<-length(unique(B$Genus))
genus_good<-length(unique(B$Genus[B$gBP>=0.30]))
family<-length(unique(B$Family))
family_good<-length(unique(B$Family[B$fBP>=0.20]))

# create df for ggplot
df<-data.frame("rank"=c("species","species","genus","genus","family","family"),
               "status"=rep(c("all","good"),3),
               "value"=c(species, species_good, genus, genus_good, family, family_good))

# create factors
df$rank = factor(df$rank, levels=c("species","genus","family"), 
                 labels=c("Species","Genus","Family"))
df$status = factor(df$status, levels=c("all","good"),
                   labels=c("All taxa","Confidently identified taxa"))

# create bar plot with two series
p<-ggplot(df, aes(fill=status, y=value, x=rank)) +
  geom_bar(position="dodge",stat="identity") +
  scale_x_discrete(limits = rev(levels(rank))) +
  labs(x="Rank", y="Unique taxa") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        legend.position = "bottom",
        legend.spacing.x = unit(0.5, 'cm'),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        axis.title.x=element_blank())
ggsave("FS3_confidentids.pdf")
