# README

This repository contains the datasets and scripts used to create the figures in:

**Hajibabaei et al., 2019. Watered-down biodiversity?  A comparison of metabarcoding results from DNA extracted from matched water and bulk tissue biomonitoring samples.**

## Data analysis outline

1. Raw reads are available from the NCBI SRA # XXX and were processed using the SCVUC v2.3 pipeline available from https://github.com/Hajibabaei-Lab/SCVUC_COI_metabarcode_pipeline . 

2. Datasets: 
  * ESVS.denoised is a FASTA file of denoised ESVs.
  * matrix.csv contains the taxonomic assignments for the denoised ESVs, as well ad read numbers per sample.

3. Figures were pepared in R as follows:
  * Fig 1 was generated with Fig1_Richness.R from matrix.csv .
  * Fig 2 was generated with Fig2_TableS1_RichnessJaccard.R from matrix.csv . 
  * Fig 3 was generated with Fig3_water_benthos_composition.R from matrix.csv . 
  * Fig 4 was generated with Fig4_NMDS.R from matrix.csvs .
  * Fig 5 was generated with Fig5_heatmap.R from matrix.csv .

4. Supplementary figures were prepared in R as follows:
  * Table S1 and S2 were generated with TableS1_S2.R with matrix.csv .
  * Table S3 was generated with TableS3_EPTO.R with matrix.csv .
  
5. Supplementary figures were prepared in R as follows:
  * Fig S1 was generated with FigS1_rarefaction.R from matrix.csv . 
  * Fig S2 was generated with FigS2_boxplot_outliers.R from matrix.csv . 
  * Fig S3 was generated with FigS3_PropConfID.R from matrixs.csv . 
  * Fig S4 was generated with FigS4_NMDS_water_benthos_separate.R from matrix.csv . 

## Acknowledgements

I would like to acknowledge funding from the Canadian government from the Genomics Research and Development Initiative (GRDI) Ecobiomics project.

Last updated: Sept. 9, 2019.
