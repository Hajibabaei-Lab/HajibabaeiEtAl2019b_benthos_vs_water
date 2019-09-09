# Teresita M. Porter, Sept. 9, 2019

library(stringr)
library(scales)
library(reshape2)
library(vegan)
library(RColorBrewer)

###################################################################
# Edit rarecurve function to remove the horizontal lines
###################################################################

rarecurve2 <- function (x, step = 1, sample, xlab = "Sample Size", ylab = "Species", 
                        label = TRUE, col, lty, ...) 
{
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function accepts only integers (counts)")
  if (missing(col)) 
    col <- par("col")
  if (missing(lty)) 
    lty <- par("lty")
  tot <- rowSums(x)
  S <- specnumber(x)
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  }
  nr <- nrow(x)
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) 
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab, 
       type = "n", ...)
  if (!missing(sample)) {
    abline(v = sample) 
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"), 
                                           y = z, xout = sample, rule = 1)$y)
    #    abline(h = rare, lwd = 0.5) #turn off horizontal lines
  }
  for (ln in seq_along(out)) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
  }
  if (label) { 
    ordilabel(cbind(tot, S), labels = rownames(x), ...)
  }
  invisible(out)
}

#####################################################################

# Read in sample x taxonomy table
A<-read.csv(file="matrix.csv", head=TRUE)

# Filter table for Arthropoda only
B<-A[A$Phylum=="Arthropoda",]

#######################################################
# Create dataframe for rarefaction curves
######################################################

# Split up SampleName for Arthropoda matrix
B2<-data.frame(B, do.call(rbind, str_split(B$SampleName,"_")))

# Filter for substrate only
B2_B<-B2[grepl("_B_",B2$SampleName),]
B2_W<-B2[grepl("_W_",B2$SampleName),]

# Add new column names
names(B2_B)[32:38]<-c("month","year","version","substrate","site","marker","siterep")
names(B2_W)[32:38]<-c("month","year","version","substrate","site","marker","siterep")

# create df for each primer+substrate to a new primer column
B2_B_AD<-B2_B[grepl("AD",B2_B$marker),]
B2_B_BE<-B2_B[grepl("BE",B2_B$marker),]
B2_W_AD<-B2_W[grepl("AD",B2_W$marker),]
B2_W_BE<-B2_W[grepl("BE",B2_W$marker),]

#combine them all into a single dataframe
B3<-rbind(B2_B_AD, B2_B_BE, B2_W_AD, B2_W_BE)

# pivot to make esv matrix (pool across verions, keep only substrate + sites separate)
esv<-dcast(B3, substrate+siterep ~ Marker_GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# merge first 3 colmns into one
esv$sample<-paste(esv$substrate,esv$siterep, sep="_")
esv<-esv[,-c(1:2)]

# move sample to rownames then delete
rownames(esv)<-esv$sample
esv$sample<-NULL

#remove columns with only zeros
esv_notnull<-esv[,colSums(esv) !=0]

#remove rows with only zeros & edit rownames
esv_notnull2<-esv_notnull[rowSums(esv_notnull) !=0,]

#calculate 15th percentile for rrarefy function
esv_15percentile<-quantile(rowSums(esv_notnull2), prob=0.15)
# 15% 
# 2099.35 

###################################################################
##### Plot rarefaction curves
###################################################################

set.seed(1234)

brewer.pal(n = 8, name = "Set1")

#print rarefaction curves for each sample to assess read coverage per sample, indicate 15th percentile
pdf("FS1_Rarecurves.pdf")

esv_rarecurveout<-rarecurve2(esv_notnull2, sample=esv_15percentile, 
                             step=100, xlab="Reads", ylab="ESVs", 
                             col=c(rep("#4DAF4A",24),rep("#377EB8",24)), 
                             cex=0.6, cex.main=0.8, cex.lab=0.6, cex.axis=0.6, 
                             label=FALSE, cex.main=1.25, cex.lab=1.25, 
                             cex.axis=1.25, lty=c(rep(1,312)), lwd=1.5)


esv_rarecurveout<-legend("bottom", ncol=2, col=c("#4DAF4A","#377EB8"), lty=c(1,1), 
               legend=c("Benthos","Water"), bty="n", pt.cex=0.5, 
               y.intersp=.75, lwd=1.5)

dev.off()

