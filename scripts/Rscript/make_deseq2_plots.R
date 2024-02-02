################################
#### Load required packages ####
################################

#Script written by Jacqueline Boccacino

library(stringr)
library(optparse)
library(dplyr)
library(ggplot2)
library(ggrepel)


############################################
#### Define functions that create plots ####
############################################

make_volcano <- function(
    dds_, # DESeq2 output results (i.e. genes with their respective log fold changes, adjusted p-values etc)
    lower_lfc_threshold_, # Lower log fold change threshold (e.g. -2)
    upper_lfc_threshold_, # Upper log fold change threshold (e.g. 2)
    gene_column_id_, # Column containing gene names
    width_, # Width for volcano plot.
    height_, # Height for volcano plot.
    outdir # Output directory where plot should be saved.
) {
  volcano <- ggplot(
    dds_,
    aes(x = LFC,
        y = minus_log10_padj,
        color = status)
  ) +
    geom_point(size = 0.5) +
    facet_wrap(~cell_type_plus_n_pbSamples, scale = "free") +
    theme_classic() +
    xlab(expression("log"[2]~"(Fold change)")) +
    ylab(expression("log"[10]~"(Adjusted p-value)")) +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold")
    ) +
    scale_color_manual(values = c("royalblue1", "snow3", "indianred1")) +
    geom_vline(xintercept = c(lower_lfc_threshold_, upper_lfc_threshold_), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(
        data = dds_[dds_$status == "Downregulated" & dds$show == TRUE,],
        aes(label = gene),
        xlim = c(NA, lower_lfc_threshold_),
        seed = 1,
        size = 2
    ) +
    geom_text_repel(
        data = dds_[dds_$status == "Upregulated" & dds$show == TRUE,],
        aes(label = gene),
        xlim = c(upper_lfc_threshold_, NA),
        seed = 1,
        size = 2
    )
  
    pdf(
        paste0(outdir, "DESeq2_volcano_", dds_comparison, "_design_", dds_design, ".pdf"),
        width = width_,
        height = height_
    )
    print(volcano)
    dev.off()
}

make_barplot <- function(
    df_, # DESeq2 output results (i.e. genes with their respective log fold changes, adjusted p-values etc)
    outdir # Output directory where plot should be saved.
) {
    bar <- ggplot(df, aes(x = cell_type, y = n_genes, fill = status)) +
        geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12)) +
    scale_fill_manual(values = c("royalblue1", "indianred1"),
                    name = "Status") +
    xlab("Cell type") +
    ylab("# Genes")

    pdf(
        paste0(outdir, "DESeq2_barplot_", dds_comparison, "_design_", dds_design, ".pdf"),
        width = 15,
        height = 10
    )
    print(bar)
    dev.off()


}

##############################
#### Get input arguments  ####
##############################

progname <- "make_deseq2_plots.R"

option_list <- list(
#    make_option(c("-i", "--input-type"), 
#        dest = "input_type", 
#        action = "store", 
#        type = "character", 
#        default = "psbulk", 
#        help = "Specify the input type (either 'psbulk' or 'other')."),
    make_option(c("-d", "--deseq2-results"), 
        dest = "dds_file", 
        action = "store", 
        type = "character", 
        default = NA, 
        help = "Path to a csv file containing DESeq2 results table (i.e. DEGs with their respective log fold changes, adjusted p-values etc)."),
    make_option(c("-e", "--deseq2-design"), 
        dest = "dds_design", 
        action = "store", 
        type = "character", 
        default = NA, 
        help = "Design used in DESeq2 analysis that generated the results table specified by -d/--deseq2-results (e.g. dataset+sex)."),
    make_option(c("-a", "--deseq2-comparison"), 
        dest = "dds_comparison", 
        action = "store", 
        type = "character", 
        default = NA, 
        help = "Comparison used in DESeq2 analysis that generated the results table specified by -d/--deseq2-results (e.g. sex_female_vs_male)."),
    make_option(c("-m", "--metadata"), 
        dest = "metadata", 
        action = "store", 
        type = "character", 
        default = NA, 
        help = "Path to a csv file containing metadata."),
    make_option(c("-s", "--metadata-row-names"), 
        dest = "metadata_rownames", 
        action = "store", 
        type = "integer", 
        default = 1, 
        help = "Position of column containing row names in metadata table [default = 1]."),
    make_option(c("-y", "--cell-type-metadata"), 
        dest = "cell_type_metadata", 
        action = "store", 
        type = "character", 
        default = NA, 
        help = "Name of column containing cell types in metadata table specified by -m/--metadata."),
    make_option(c("-n", "--donor-metadata"), 
        dest = "donor_metadata", 
        action = "store", 
        type = "character", 
        default = NA, 
        help = "Name of column containing donors in metadata table specified by -m/--metadata."),
    make_option(c("-l", "--lfc-column-id"), 
        dest = "lfc_column_id", 
        action = "store", 
        type = "character", 
        default = "log2FoldChange", 
        help = "Name of column containing log fold changes in DESeq2 results table specified by -d/--deseq2-results."),
    make_option(c("-1", "--lower-lfc-threshold"), 
        dest = "lower_lfc_threshold", 
        action = "store", 
        type = "double", 
        default = -2.0, 
        help = "Lower log fold change threshold (default = -2.0)."),
    make_option(c("-2", "--upper-lfc-threshold"), 
        dest = "upper_lfc_threshold", 
        action = "store", 
        type = "double", 
        default = 2.0, 
        help = "Upper log fold change threshold (default = 2.0)."),
    make_option(c("-p", "--padj-column-id"), 
        dest = "padj_column_id", 
        action = "store", 
        type = "character", 
        default = "padj", 
        help = "Name of column containing adjusted p-values in DESeq2 results table specified by -d/--deseq2-results."),
    make_option(c("-c", "--cell-type-dds"), 
        dest = "cell_type_dds", 
        action = "store", 
        type = "character", 
        default = NA, 
        help = "Name of column containing cell types in DESeq2 results table specified by -d/--deseq2-results."),
    make_option(c("-g", "--gene-column-id"), 
        dest = "gene_column_id", 
        action = "store", 
        type = "character", 
        default = "gene", 
        help = "Name of column containing gene names in DESeq2 results table specified by -d/--deseq2-results."),
    make_option(c("--volcano-width"), 
        dest = "volcano_width", 
        action = "store", 
        type = "double", 
        default = 30.0, 
        help = "Volcano plot width."),
    make_option(c("--volcano-height"), 
        dest = "volcano_height", 
        action = "store", 
        type = "double", 
        default = 30.0, 
        help = "Volcano plot height."),
    make_option(c("-r", "--genes-to-remove"), 
        dest = "genes_to_remove", 
        action = "store", 
        type = "character", 
        default = NA, 
        help = "Path to a txt file containing genes that should not be displayed in the volcano plot; each line should be a gene."),
    make_option(c("-o", "--output-dir"), 
        dest = "output_dir", 
        action = "store", 
        type = "character", 
        default = NA, 
        help = "Path to output directory.")
)

parser <- OptionParser(
    usage = paste0(progname, " [options] 
#    --input-type psbulk
    --deseq2-results PATH/TO/RESULTS
    --deseq2-design FACTOR1,FACTOR2...
    --deseq2-comparison CONDITION_INTEREST_VS_CONTROL
    --metadata PATH/TO/METADATA 
    --metadata-row-names 1
    --cell-type-metadata METADATA_CELL_TYPE_COLUMN
    --donor-metadata METADATA_DONOR_COLUMN
    --lfc-column-id LFC_COLUMN
    --lower-lfc-threshold X
    --upper-lfc-threshold Y
    --padj-column-id PADJ_COLUMN
    --cell-type-dds DDS_CELL_TYPE_COLUMN
    --gene-column-id gene
    --volcano-width 30
    --volcano-height 30
    --genes-to-remove PATH/TO/FILE
    --output-dir PATH/TO/OUTPUT_DIR"), 
    option_list = option_list)

arguments <- parse_args(parser, positional_arguments = 0)

#input_type <- arguments$options$input_type
dds_file <- arguments$options$dds_file
dds_design <- arguments$options$dds_design
dds_comparison <- arguments$options$dds_comparison
metadata <- arguments$options$metadata
metadata_rownames <- arguments$options$metadata_rownames
cell_type_metadata <- arguments$options$cell_type_metadata
donor_metadata <- arguments$options$donor_metadata
lfc_column_id <- arguments$options$lfc_column_id
lower_lfc_threshold <- arguments$options$lower_lfc_threshold
upper_lfc_threshold <- arguments$options$upper_lfc_threshold
padj_column_id <- arguments$options$padj_column_id
cell_type_dds <- arguments$options$cell_type_dds
gene_column_id <- arguments$options$gene_column_id
volcano_width <- arguments$options$volcano_width
volcano_height <- arguments$options$volcano_height
genes_to_remove <- arguments$options$genes_to_remove
output_dir <- arguments$options$output_dir


#####################################################
#### Check integrity of input files/directories  ####
#####################################################

#if(is.na(input_type)) warning("No input type has been provided. Default 'psbulk' will be used.")
if(is.na(dds_file)) stop("No DESeq2 results table has been provided. Please review your input.")
if(is.na(dds_design)) stop("No DESeq2 design has been provided. Please review your input.")
if(is.na(dds_comparison)) stop("No DESeq2 comparison has been provided. Please review your input.")
if(is.na(metadata)) stop("No metadata table has been provided. Please review your input.")
if(is.na(metadata_rownames)) stop("No metadata row names position has been provided. Please review your input.")
if(is.na(cell_type_metadata)) stop("No cell type column name has been provided for the metadata table. Please review your input.")
if(is.na(donor_metadata)) stop("No donor column name has been provided for the metadata table. Please review your input.")
if(is.na(lfc_column_id)) warning("No log fold change column name has been provided for the DESeq2 results table. The default name (log2FoldChange) will be used.")
if(is.na(lower_lfc_threshold)) warning("No lower fold change threshold has been provided. Default value of -2.0 will be used.")
if(is.na(upper_lfc_threshold)) warning("No upper fold change threshold has been provided. Default value of 2.0 will be used.")
if(is.na(padj_column_id)) warning("No upper fold change threshold has been provided. The default name (padj) will be used.")
if(is.na(cell_type_dds)) stop("No cell type column name has been provided for the DESeq2 results table. Please review your input.")
if(is.na(gene_column_id)) warning("No gene column name has been provided. The default name (gene) will be used.")
if(is.na(volcano_width)) warning("No width for the volcano plot has been provided. No genes will be removed from the plots.")
if(is.na(volcano_height)) warning("No height for the volcano plot has been provided. No genes will be removed from the plots.")
if(is.na(genes_to_remove)) warning("No file containing genes to remove has been provided. No genes will be removed from the plots.")
if(is.na(output_dir)) stop("No output directory has been provided. Please review your input.")

# Check if DESeq2 and metadata files exist.
if(!file.exists(dds_file)) stop(paste0(dds_file, " does not exist. Please review your input."))
if(!file.exists(metadata)) stop(paste0(metadata, " does not exist. Please review your input."))
# Check if output directory exists.
if(!dir.exists(output_dir)) { 
    stop(paste0(output_dir, " does not exist. Please review your input."))
} else {
    # If output directory exists, ensure that last character is a slash symbol.
    if(!(substr(output_dir, nchar(output_dir), nchar(output_dir)) == "/")) {
        output_dir <- paste0(output_dir, "/")
    }
    if(!dir.exists(paste0(output_dir, "plots/"))) {
        cat("Creating", paste0(output_dir, "plots/"), "\n")
        dir.create(paste0(output_dir, "plots/"))
    } else {
        cat(paste0(output_dir, "plots/"), "already exists.\n")
    }
}


#####################################
#### Read & process input files  ####
#####################################

# Read DESeq2 results.
dds <- read.csv(
    dds_file, 
    header = TRUE,
    row.names = 1
)
# Check if the columns specified actually exist in the DESeq2 results table provided.
columns_to_check_in_dds <- c(lfc_column_id, padj_column_id, cell_type_dds, gene_column_id)
for (col in columns_to_check_in_dds) {
    if(!(col %in% colnames(dds))) {
        stop(paste0(col, " is not a column in the DESeq2 results table provided. Please review your input."))
    }
}

# Read metadata table, where rows are pseudo-bulk samples and columns are observations/info on pseudo-bulk samples.
metadata_table <- read.csv(
    metadata, 
    header = TRUE,
    row.names = metadata_rownames,
    check.names = FALSE
)
# Check if the columns specified actually exists in the metadata table provided.
columns_to_check_in_metadata_table <- c(cell_type_metadata, donor_metadata)
for (col in columns_to_check_in_metadata_table) {
    if(!(col %in% colnames(metadata_table))) {
        stop(paste0(col, " is not a column in the metadata table provided. Please review your input."))
    }
}

# Read list of genes to remove.
genelist_remove <- NA
if(!is.na(genes_to_remove)) {
    genelist_remove <- read.csv(
        genes_to_remove, 
        header = FALSE,
        row.names = NULL,
        sep = "\t"
    )
    genelist_remove <- as.vector(genelist_remove[,1])
}

# Check if log fold changes specified as thresholds fall within the range available in the DESeq2 results table.
range_lfc <- range(dds[,lfc_column_id])
for (threshold in c(lower_lfc_threshold, upper_lfc_threshold)) {
    if(!between(threshold, range_lfc[1], range_lfc[2])) {
        stop(paste0(threshold, " is not within the log fold change range in the DESeq2 results table provided. Please choose a different value."))
    }
}

# Create column in the DESeq2 results table that will contain cell type + number of pseudo-bulk samples.
metadata_table$n_pbSamples_perCellType <- sapply(1:nrow(metadata_table), function(x) {
    return(metadata_table[,donor_metadata][metadata_table[,cell_type_metadata] == metadata_table[,cell_type_metadata][x]] %>% unique %>% length)
})
metadata_table$cell_type_plus_n_pbSamples <- paste0(
    metadata_table[,cell_type_metadata], " (n=", metadata_table$n_pbSamples_perCellType, ")"
)
dds$cell_type_plus_n_pbSamples <- metadata_table$cell_type_plus_n_pbSamples[match(dds[,cell_type_dds], metadata_table[,cell_type_metadata])]

# Create new columns based on user-provided columns.
dds$LFC <- as.numeric(dds[,lfc_column_id])
dds$gene <- dds[,gene_column_id]

# Log-transform adjusted p-values.
dds$minus_log10_padj <- -log10(dds[,padj_column_id])

# Create gene classification column (upregulated/downregulated/non-significant).
dds <- dds %>% mutate(status = case_when(
  LFC >= upper_lfc_threshold & minus_log10_padj >= -log10(0.05) ~ "Upregulated",
  LFC <= lower_lfc_threshold & minus_log10_padj >= -log10(0.05) ~ "Downregulated",
  LFC > lower_lfc_threshold & LFC < upper_lfc_threshold ~ "Non-significant",
  minus_log10_padj < -log10(0.05) ~ "Non-significant"
))

# Create column to determine whether gene name should be displayed in volcano plot or not.
if(!is.na(genes_to_remove)) {
    dds <- dds %>% mutate(show = case_when(
    !(gene %in% genelist_remove) ~ TRUE,
    gene %in% genelist_remove & gene != "XIST" & gene != "RPS4Y1" ~ FALSE,
    gene == "XIST" | gene == "RPS4Y1" ~ TRUE
    ))
} else {
    dds$show <- rep(TRUE, times = nrow(dds))
}

# Call function to create volcano & save plot to output file.
make_volcano(
    dds_ = dds,
    lower_lfc_threshold_ = lower_lfc_threshold,
    upper_lfc_threshold_ = upper_lfc_threshold, 
    width_ = volcano_width,
    height_ = volcano_height,
    outdir = paste0(output_dir, "plots/")
)

# Create dataframe that contains the number of upregulated/downregulated genes per cell type.
df <- data.frame(matrix(nrow = 0, ncol = 0))
for(cell_type in dds$cell_type %>% unique) {
    n_up_genes <- dds$gene[dds$status == "Upregulated" & dds$cell_type == cell_type] %>% length
    n_down_genes <- dds$gene[dds$status == "Downregulated" & dds$cell_type == cell_type] %>% length
  
    df <- rbind(df, c(cell_type, "Upregulated", n_up_genes))  
    df <- rbind(df, c(cell_type, "Downregulated", n_down_genes))  
}
colnames(df) <- c("cell_type", "status", "n_genes")
df$n_genes <- as.numeric(df$n_genes) 

# Call function to create bar plot & save plot to output file.
make_barplot(
    df_ = df,
    outdir = paste0(output_dir, "plots/")
)