# Teresita M. Porter, Sept. 9, 2019

library(reshape2)
library(plyr)
library(scales)

###################################################################

# Read infile prepared by python script
A<-read.table(file="matrix.csv", head=TRUE, sep=",")

# Get all ESV counts
length(unique(A$Marker_GlobalESV))
#16,841
length(unique(A$Marker_GlobalESV[grepl("^AD",A$Marker_GlobalESV) &
                                 grepl("_B_",A$SampleName)]))
#2,581
length(unique(A$Marker_GlobalESV[grepl("^AD",A$Marker_GlobalESV) &
                                   grepl("_W_",A$SampleName)]))
#1,099
length(unique(A$Marker_GlobalESV[grepl("^BE",A$Marker_GlobalESV) &
                                   grepl("_B_",A$SampleName)]))
#4,831
length(unique(A$Marker_GlobalESV[grepl("^BE",A$Marker_GlobalESV) &
                                   grepl("_W_",A$SampleName)]))
#9,571

# Get all ESV read counts
sum(A$ESVsize)
#5,407,720
sum(A$ESVsize[grepl("^AD",A$Marker_GlobalESV) &
              grepl("_B_",A$SampleName)])
#2,614,558
sum(A$ESVsize[grepl("^AD",A$Marker_GlobalESV) &
                grepl("_W_",A$SampleName)])
#238,697
sum(A$ESVsize[grepl("^BE",A$Marker_GlobalESV) &
                grepl("_B_",A$SampleName)])
#1,774,386
sum(A$ESVsize[grepl("^BE",A$Marker_GlobalESV) &
                grepl("_W_",A$SampleName)])
#780,079

# Select all Arthropoda
B<-A[A$Phylum=="Arthropoda",]

length(unique(B$Marker_GlobalESV))
#4,459
length(unique(B$Marker_GlobalESV[grepl("^AD",B$Marker_GlobalESV) &
                                   grepl("_B_",B$SampleName)]))
#1,735
length(unique(B$Marker_GlobalESV[grepl("^AD",B$Marker_GlobalESV) &
                                   grepl("_W_",B$SampleName)]))
#280
length(unique(B$Marker_GlobalESV[grepl("^BE",B$Marker_GlobalESV) &
                                   grepl("_B_",B$SampleName)]))
#2,398
length(unique(B$Marker_GlobalESV[grepl("^BE",B$Marker_GlobalESV) &
                                   grepl("_W_",B$SampleName)]))
#491

# Get all ESV read counts
sum(B$ESVsize)
#4,399,949
sum(B$ESVsize[grepl("^AD",B$Marker_GlobalESV) &
                grepl("_B_",B$SampleName)])
#2,541,062
sum(B$ESVsize[grepl("^AD",B$Marker_GlobalESV) &
                grepl("_W_",B$SampleName)])
#174,605
sum(B$ESVsize[grepl("^BE",B$Marker_GlobalESV) &
                grepl("_B_",B$SampleName)])
#1,554,853
sum(B$ESVsize[grepl("^BE",B$Marker_GlobalESV) &
                grepl("_W_",B$SampleName)])
#129,429
