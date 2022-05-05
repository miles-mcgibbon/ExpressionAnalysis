
# check all required packages are installed before running the script
dependencies <- c("pheatmap","rcartocolor","assertthat","dplyr")
not_installed <- setdiff(dependencies, rownames(installed.packages()))

# function to get user input for package installations
# should work in both RStudio and Rscript
user_input <- function(msg) {
  if (interactive()) {
    response <- readline(msg)
  } else {
    cat(msg);
    response <- readLines("stdin",n=1);
  }
  return(response)
}


# if needed packages are not installed then ask the user if they want to
# install them or quit
if (length(not_installed) != 0){
  for (package in not_installed) {
    print(paste("WARNING: Essential package", package, "is not installed"))
    answer <- user_input("Install now? (y/N)")
    if (answer %in% c('y','Y')){
      install.packages(package, repos='http://cran.uk.r-project.org')
    }
    else {
      print("Exiting...")
      quit(status=1)
    }
  }
}

# load required packages
library(pheatmap)
library(rcartocolor)
library(assertthat)
library(dplyr)

# define where to find experiment files
expression_data_file <- "sample_data/data_all.csv"
gene_annotation_file <- "sample_data/gene_annotation.csv"
sample_annotation_file <- "sample_data/sample_annotation.csv"
gene_list_file <- "sample_data/genelist.txt"

########
# make sure all these files above exist in a loop
########

# for each file above
for (data_file in c(expression_data_file, gene_annotation_file, sample_annotation_file, gene_list_file)) {

  # check that this file exists
  check <- assert_that(file.exists(data_file),

    # give an error if not
    msg = paste(
      "Path '",
      data_file,
      "' does not exist\n Please check all input filepaths are correct!"
    )
  )
}


# read in the expression data csv
data <- read.csv(expression_data_file)

########
# check that all the columns apart from the first column with gene names are numeric
# before applying log2 scaling
########

# get all the columns that are numeric excluding first column containing gene names
numeric_columns <- colnames(select_if(data[, 2:length(colnames(data))], is.numeric))

# get all the columns excluding the first column containing gene names
data_columns <- colnames(data[, 2:length(colnames(data))])

# check that numeric_columns and data_columns are equal and give an error if not
# showing which column(s) contain non-numeric values
check_transformable <- assert_that(length(setdiff(data_columns, numeric_columns)) == 0,
  msg = paste("Expression data in column(s) '",
    paste(unlist(setdiff(data_columns, numeric_columns)), collapse = ", "),
    "' from '",
    expression_data_file,
    "' contain non-numeric values",
    sep = ""
  )
)

########
# apply log2 scaling to the expression data
########

# convert the whole dataframe to a matrix
data_matrix <- as.matrix(data)

# log2 scale all columns except the first (the gene name column)
data_matrix[, 2:length(colnames(data_matrix))] <- log2(data_matrix[, 2:length(colnames(data_matrix))] + 1)

########
# read in the other informative files to select and annotate expression data
########

# read in the gene and sample annotation files
gene_annot <- read.csv(gene_annotation_file)
sample_annot <- read.csv(sample_annotation_file)

# read custom gene list
genelist <- read.table(gene_list_file, header = TRUE)

# ensure genelist is dataframe with only one column
genelist <- genelist[, 1]

# perform soft validation check for duplicates in gene list and give warning
# if any are found
dup_warning <- validate_that(sum(duplicated(genelist)) == 0,
                             msg = paste("WARNING: ", sum(duplicated(genelist)),
                                         " duplicate(s) found in supplied gene list file '",
                                         gene_list_file,
                                         sep = ""
                             )
)

if (dup_warning != TRUE) {
  print(dup_warning)
  print('WARNING: All duplicates will be removed before analysis')
}

# remove any duplicates from the genelist
genelist <- genelist[!duplicated(genelist)]


########
# subset the expression data for only genes supplied in genelist
########

# get matrix of only gene expression data we want using genelist
keep_data <- data_matrix[genelist, ]
keep_data

# turn into dataframe for merging with gene annotation data
keep_data <- as.data.frame(keep_data)

# check that the length of keep_data expression data dataframe is the same
# as the length of the genelist and if not warn the user
gene_check <- validate_that(length(rownames(keep_data)) == length(genelist),
  msg = paste("WARNING: ",
    length(genelist) - length(rownames(keep_data)),
    " gene(s) in the gene list file NOT FOUND in expression data",
    sep = ""
  )
)

if (gene_check != TRUE) {
  print(gene_check)
}

########
# prepare sample annotation file to annotate the heatmaps
########

# check consistent column naming in the supplied sample annotation file
column_check <- assert_that("TreatmentGroup" %in% names(sample_annot), 
                            msg = paste("No TreatmentGroup column found in supplied sample annotation file '",
                                        sample_annotation_file,
                                        "'",
                                        sep = "")
)

# transform sample annotation data for plotting as annotation
# by assigning sample names as rownames to match plotted data
rownames(sample_annot) <- sample_annot$SampleName

# then keep only TreatmentGroup column we need for annotating the plots
sample_annot <- sample_annot[, "TreatmentGroup", drop = FALSE]

# rename TreatmentGroup column to Treatment to take less space when plotting
names(sample_annot)[names(sample_annot) == "TreatmentGroup"] <- "Treatment"
sample_annot$Treatment <- as.character(sample_annot$Treatment)

########
# add gene annotation information to gene expression data
########

# check that the X column with gene number exists in both expression data
# and gene annotation data before merging
gene_column_check <- assert_that("Gene" %in% names(gene_annot),
                                 msg = paste("No gene ID column 'Gene' found in supplied gene annotation file '",
                                             gene_annotation_file,
                                             "'",
                                             sep = "")
)

keep_column_check <- assert_that("Gene" %in% names(keep_data),
                                 msg = paste("No gene ID column 'Gene' found in supplied gene expression file '",
                                             expression_data_file,
                                             "'",
                                             sep = "")
)

# merge the gene annotation dataframe with the gene expression dataframe
keep_data_info <- merge(keep_data, gene_annot, by = "Gene")

# check that no columns have been lost from the expression data due to the merge
all_data_kept_check <- assert_that(length(rownames(keep_data_info)) == length(rownames(keep_data)),
                                   msg = paste(length(rownames(keep_data)) - length(rownames(keep_data_info)),
                                               " gene row(s) lost when merging gene annotation data")
)

########
# create cleaned matrix from merged dataframe for plotting
########

# assign long gene names as merged dataframe rownames for plotting
rownames(keep_data_info) <- keep_data_info$LongName

# rename 'Type' column to 'GeneType' for clearer plot
names(keep_data_info)[names(keep_data_info) == "Type"] <- "GeneType"

# define columns to remove from the dataframe before plotting
drop_cols <- c("X", "GeneType", "Gene", "LongName")

# remove the predefined columns above
plotting_df <- keep_data_info[, !(names(keep_data_info) %in% drop_cols)]

# convert the dataframe to a matrix
plotting_matrix <- as.matrix(plotting_df)

# check that all the data (excluding colnames and rownames) is numeric and can
# be plotted with no errors
plottable_check <- assert_that(is.numeric(plotting_matrix),
                               msg = "Non-numeric values found in gene expression data - this data must be removed before plotting")



########
# create colour palette for different treatments and gene types
# depending on number of treatments and gene types in supplied data
########

# count number of treatments and gene types
treatment_levels <- length(unique(sample_annot$Treatment))
gene_types <- length(unique(keep_data_info$GeneType))

# get total number of distinct colours needed
total_variables <- treatment_levels + gene_types

# get a palette of distinct colours for the total number of variables
# plus one as the carto_pal function always returns grey as the final
# colour
distinct_colours <- carto_pal(total_variables + 1, "Pastel")

# drop grey from the returned colour palette
distinct_colours <- distinct_colours[1:total_variables]

# define empty nested named list to populate with colours
plot_colours <- list(
  "Treatment" = NULL,
  "GeneType" = NULL
)

# populate the plot_colours list with the colour palette
plot_colours[[1]] <- distinct_colours[1:treatment_levels]
plot_colours[[2]] <- distinct_colours[treatment_levels + 1:total_variables]
plot_colours[[2]] <- plot_colours[[2]][!is.na(plot_colours[[2]])]

# add names of treatments and gene types to the nested named list
names(plot_colours[[1]]) <- unique(sample_annot$Treatment)
names(plot_colours[[2]]) <- unique(keep_data_info$GeneType)


########
# create heatmaps of data
########

# create output directory and ignore warnings if it exists already
dir.create('output', showWarnings = FALSE)

# heatmap one with only genes clustered

pheatmap(plotting_matrix,
  scale = "row",
  annotation_row = keep_data_info[, "GeneType", drop = FALSE],
  annotation_col = sample_annot,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  fontsize_row = 6,
  cellheight = 6,
  annotation_colors = plot_colours,
  legend = FALSE,
  main = "Log2 Row-wise Expression Level Clustered By Gene",
  fontsize_col = 12,
  angle_col = 45,
  filename = 'output/Log2_Row-wise_Expression_Level_Clustered_By_Gene.png',
)

# heatmap two with genes and samples clustered

pheatmap(plotting_matrix,
  scale = "row",
  annotation_row = keep_data_info[, "GeneType", drop = FALSE],
  annotation_col = sample_annot,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  fontsize_row = 6,
  treeheight_col = 5,
  annotation_colors = plot_colours,
  legend = FALSE,
  main = "Log2 Row-wise Expression Level Clustered\n By Gene and Sample",
  fontsize_col = 12,
  angle_col = 45,
  cutree_cols = 4,
  filename = 'output/Log2_Row-wise_Expression_Level_Clustered_By_Gene_and_Sample.png'
)

print("Analysis Finished! Heatmaps saved to output/ dir")

