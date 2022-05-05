# R Gene Expression Analysis Pipeline

---
This repo contains an R pipeline to plot clustered heatmaps from gene expression data. Four files are needed as input:

### Inputs

---

1. Raw gene expression data as a CSV file:

|Gene|A  |B  |C                         |D          |E          |F          |G          |H          |I          |J          |K          |L          |
|----|---|---|--------------------------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
|1   |6.24234931|6.302402541|6.584037232               |6.645032662|6.591125959|6.536004289|4.465150861|4.175000316|4.587409813|2.813028388|2.597595679|2.805375832|
|2   |5.404286196|5.667918769|5.73193275                |8.686235318|8.841933639|8.728939674|8.99805177 |8.436325549|8.691272211|3.834652157|3.632691162|3.741054209|
|3   |2.648436021|2.68515117|2.727255144               |2.518935413|2.41201574 |2.422645317|2.345908774|2.284605432|2.514958181|5.299665666|5.340035796|5.652512044|


2. A gene annotation CSV file with the following structure:

|Gene|Type|LongName                  |
|----|----|--------------------------|
|1   |XA  |1436799_at                |
|2   |XA  |1436227_at                |
|3   |XA  |1420504_at                |

3. A list of indexes of genes of interest as a single column TXT file with the following structure:

|Gene|
|--|
|63|
|104|
|13|

4. A sample treatment annotation CSV file with the following structure:

|Gene|SampleName|TreatmentGroup|
|----|----------|--------------|
|1   |A         |1             |
|2   |B         |1             |
|3   |C         |1             |

### Outputs

---

The pipeline will take the expression data for the supplied genes of interest and produce two heatmaps detailing expresssion level, gene names, treatments and gene types. One is clustered by each gene and sample:

![]('sample_output/Log2_Row-wise_Expression_Level_Clustered_By_Gene_and_Sample.png')


The other is clustered only by gene:

![]('Log2_Row-wise_Expression_Level_Clustered_By_Gene.png')

### Markdown Report

The repo also contains an Rmd script, 'Generate_Report.Rmd', which generates the heatmaps as part of a PDF report. The PDF also presents information on input files, unique gene types, gene counts, number of samples and number of treatments.
