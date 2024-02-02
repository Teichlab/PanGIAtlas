################################
#### Load required packages ####
################################

#Script written by Jacqueline Boccacino

library(DESeq2)
library(stringr)
library(optparse)
library(AnnotationDbi)
library(org.Hs.eg.db)


##############################################################
#### Define function to perform DGE analysis with DESeq2  ####
##############################################################

perform_DGE <- function(
    cdata, # Count matrix
    mdata, # Metadata table
    columns_for_design_, # Columns for DESeq2 analysis design
    column_to_compare_, # Column that contains the condition of interest and control condition (e.g. sex)
    condition_of_interest_, # Condition we are interested in (e.g. female)
    control_condition_, # Control condition (e.g. male)
    cell_type_id, # Cell type of interest
    outdir # Output directory for saving per-cell type DESeq2-normalised count tables.
) {
    # This chunck of the code is designed like this to prevent the script from bailing out
    # in cases when the design matrix is not full rank; the goal is to just skip to the next
    # cell type instead of just stopping the program...
    tryCatch({

        # Build a design formula for downstream DESeq2 analysis using the columns for design defined above.
        analysis_design <- as.formula(
            paste("~", paste(columns_for_design_, collapse= "+"))
        )
        
        # Build DESeq2 object based on a given design.
        dds <- DESeqDataSetFromMatrix(
            countData = cdata,
            colData = mdata,
            design = analysis_design
        )
        # Make sure the comparison is made against the control condition.
        colData(dds)[,column_to_compare_] <- relevel(
            colData(dds)[,column_to_compare_], 
            ref = control_condition_
        )
        # Estimate size factors.
        dds <- estimateSizeFactors(dds)
        # Run DESeq2.
        dds <- DESeq(dds)
        # Extract results for the comparison of interest.
        res <- results(
            dds, 
            contrast = c(column_to_compare_, condition_of_interest_, control_condition_)
        )
        # Perform log fold chance shrinkage.
        resLFC <- lfcShrink(
            dds, 
            res = res,
            coef = paste(column_to_compare_, condition_of_interest_, "vs", control_condition_, sep = "_"),
            type = "ashr"
        )
        # Filter out non-significant results according to a given adjusted p-value threshold.
        resLFC <- as.data.frame(resLFC[!is.na(resLFC$padj),])
        # Add cell type to dataframe.
        resLFC$cell_type <- rep(cell_type_id, times = nrow(resLFC))

        # Save cell type-specific DESeq2-normalised count table to output file.
        write.csv(
            counts(dds, normalized=TRUE),
            paste0(
                outdir, "DESeq2_", column_to_compare_, "_", condition_of_interest_, "_vs_", control_condition_, "_",
                cell_type_id, "_design_", paste(columns_for_design_, collapse = "_"), "_normalisedCounts.csv"
            )
        )

        # Return DGE results.
        return(resLFC)
    }, error = function(e){
        cat(
            cell_type, "failed. Reason:\n", conditionMessage(e), "\n"
        )
    })
}


##############################
#### Get input arguments  ####
##############################

progname <- "run_dge_deseq2_per_celltype.r"

option_list <- list(
    make_option(c("-x", "--count-matrix"), 
        dest = "counts", 
        action = "store", 
        type = "character", 
        default = NA, 
        help = "Path to input count matrix."),
    make_option(c("-l", "--sample-location"), 
        dest = "sample_loc", 
        action = "store", 
        type = "character", 
        default = "rows", 
        help = "Whether samples in count matrix are rows or columns [default = rows]."),
    make_option(c("-m", "--metadata"), 
        dest = "metadata", 
        action = "store", 
        type = "character", 
        default = NA, 
        help = "Path to input metadata table."),
    make_option(c("-n", "--metadata-row-names"), 
        dest = "metadata_rownames", 
        action = "store", 
        type = "integer", 
        default = 1, 
        help = "Position of column containing row names in metadata table [default = 1]."),
    make_option(c("-o", "--output-dir"), 
        dest = "output_dir", 
        action = "store", 
        type = "character", 
        default = NA, 
        help = "Path to output directory."),
    make_option(c("-p", "--output-prefix"), 
        dest = "output_prefix", 
        action = "store", 
        type = "character", 
        default = NA, 
        help = "Prefix for output files."),
    make_option(c("-c", "--cell-type-column-id"), 
        dest = "cell_type_column_id", 
        action = "store", 
        type = "character", 
        default = NA, 
        help = "Column name in AnnData.obs that corresponds to cell types/states/identities."),
    make_option(c("-d", "--columns-for-design"), 
        dest = "columns_for_design", 
        action = "store", 
        type = "character", 
        default = NA, 
        help = "Comma-separated list of AnnData.obs column names to be used in the design of DESeq2 analysis."),
    make_option(c("-r", "--column-to-compare"), 
        dest = "column_to_compare", 
        action = "store", 
        type = "character", 
        default = NA, 
        help = "Column name in AnnData.obs that contains the conditions of interest that should be compared."),
    make_option(c("-1", "--condition-of-interest"), 
        dest = "condition_of_interest", 
        action = "store", 
        type = "character", 
        default = NA, 
        help = "Name of condition of interest; must be present in the column specified by -r/--column-to-compare."),
    make_option(c("-2", "--control-condition"), 
        dest = "control_condition", 
        action = "store", 
        type = "character", 
        default = NA, 
        help = "Name of control condition; must be present in the column specified by -r/--column-to-compare.")
)

parser <- OptionParser(
    usage = paste0(progname, " [options] 
    --count-matrix PATH/TO/COUNTS 
    --count-location rows
    --metadata PATH/TO/METADATA 
    --metadata-row-names 1
    --output-dir PATH/TO/DIR
    --output-prefix PREFIX
    --cell-type-column-id CELL_TYPE 
    --columns-for-design COLUMN1,COLUMN2,COLUMN3 
    --column-to-compare COLUMN 
    --condition-of-interest INTEREST 
    --control-condition CONTROL"), 
    option_list = option_list)

arguments <- parse_args(parser, positional_arguments = 0)

counts <- arguments$options$counts
sample_loc <- arguments$options$sample_loc
metadata <- arguments$options$metadata
metadata_rownames <- arguments$options$metadata_rownames
output_dir <- arguments$options$output_dir
output_prefix <- arguments$options$output_prefix
cell_type_column_id <- arguments$options$cell_type_column_id
columns_for_design <- arguments$options$columns_for_design
column_to_compare <- arguments$options$column_to_compare
condition_of_interest <- arguments$options$condition_of_interest
control_condition <- arguments$options$control_condition

# Check if sample position (rows/cols) in count matrix has been provided.
if(is.na(sample_loc)) cat("No sample location (rows/cols) has been provided. Proceeding with default (rows).")

# Check if counts and metadata files exist.
if(!file.exists(counts)) stop(paste0(counts, " does not exist. Please review your input."))
if(!file.exists(metadata)) stop(paste0(metadata, " does not exist. Please review your input."))
# Check if output directory exists.
if(!dir.exists(output_dir)) { 
    stop(paste0(output_dir, " does not exist. Please review your input."))
} else {
    # If output directory exists, ensure that last character is a slash symbol.
    if(!(substr(output_dir, nchar(output_dir), nchar(output_dir)) == "/")) {
        output_dir <- paste0(output_dir, "/")
    }
    # Create output directory for normalised counts files.
    norm_counts_oudir <- paste0(output_dir, "normalised_counts/")
    if(!dir.exists(norm_counts_oudir)) {
        cat("Creating", norm_counts_oudir, "\n")
        dir.create(norm_counts_oudir)
    } else {
        cat(norm_counts_oudir, "already exists.\n")
    }
    # Create output directory for DEG tables.
    deg_oudir <- paste0(output_dir, "degs/")
    if(!dir.exists(deg_oudir)) {
        cat("Creating", deg_oudir, "\n")
        dir.create(deg_oudir)
    } else {
        cat(deg_oudir, "already exists.\n")
    }
}


#####################################
#### Read & process input files  ####
#####################################

# Get format of count matrix file (e.g. 'csv', 'txt'...)
elements_in_counts_path <- unlist(str_split(counts, "\\."))
counts_file_format <- elements_in_counts_path[length(elements_in_counts_path)]

# Read count matrix 
count_matrix <- NA
if(counts_file_format == "csv") {
    count_matrix <- read.csv(
        counts, 
        header = TRUE,
        row.names = 1,
        check.names = FALSE
    )
} else {
    count_matrix <- read.table(
        counts, 
        header = TRUE,
        row.names = 1,
        check.names = FALSE
    )
}

# Transpose it so that rows are genes and columns are pseudo-bulk samples.
if(sample_loc == "rows") {
    count_matrix <- as.data.frame(t(count_matrix))
}

# Read metadata table, where rows are pseudo-bulk samples and columns are observations/info on pseudo-bulk samples.
metadata_table <- read.csv(
    metadata, 
    header = TRUE,
    row.names = metadata_rownames,
    check.names = FALSE
)
# Check if the cell type column specified actually exists in the metadata table provided.
if(!(cell_type_column_id %in% colnames(metadata_table))) {
    stop(paste0(cell_type_column_id, " is not a column in the metadata table provided. Please review your input."))
}
# Check if column to compare actually exists in the metadata table provided.
if(!(column_to_compare %in% colnames(metadata_table))) {
    stop(paste0(column_to_compare, " is not a column in the metadatcoua table provided. Please review your input."))
} else { 
    # If the column to compare exists, then check if the condition of interest specified exists under that column.
    if(!(condition_of_interest %in% metadata_table[,column_to_compare])) {
        stop(paste0(condition_of_interest, " is not present under ", column_to_compare, ". Please review your input."))
    }
    # Also check if the control condition specified exists under that column.
    if(!(control_condition %in% metadata_table[,column_to_compare])) {
        stop(paste0(control_condition, " is not present under ", column_to_compare, ". Please review your input."))
    }
}
# Ensure that row names in metadata table have the same order as column names in the count matrix.
metadata_table <- metadata_table[colnames(count_matrix),]
# Ensure that the elements under column_to_compare do not have any spaces or dashes - replace those with underscores!
metadata_table[,column_to_compare] <- gsub(" ", "_", metadata_table[,column_to_compare])
metadata_table[,column_to_compare] <- gsub("-", "_", metadata_table[,column_to_compare])

# Create list from user-defined comma-separated columns.
columns_for_design <- unlist(
    str_split(columns_for_design, ",")
)
for(col in columns_for_design) {
    # Check if columns for design actually exist in the metadata table provided.
    if(!(col %in% colnames(metadata_table))) stop(paste0(col, " is not a column in the metadata table provided. Please review your input."))
    # Ensure that the elements under column do not have any spaces or dashes - replace those with underscores!
    metadata_table[,col] <- gsub(" ", "_", metadata_table[,col])
    metadata_table[,col] <- gsub("-", "_", metadata_table[,col])
}

# Ensure that condition_of_interest and control_condition do not have any spaces - replace those with underscores!
condition_of_interest <- gsub(" ", "_", condition_of_interest)
condition_of_interest <- gsub("-", "_", condition_of_interest)
control_condition <- gsub(" ", "_", control_condition)
control_condition <- gsub("-", "_", control_condition)

# Create directories for specified design.
# Normalised counts:
norm_counts_design_oudir <- paste0(output_dir, "normalised_counts/", paste(columns_for_design, collapse = "_"), "/")
if(!dir.exists(norm_counts_design_oudir)) {
    cat("Creating", norm_counts_design_oudir, "\n")
    dir.create(norm_counts_design_oudir)
} else {
    cat(norm_counts_design_oudir, "already exists.\n")
}
# DEGs:
deg_design_oudir <- paste0(output_dir, "degs/", paste(columns_for_design, collapse = "_"), "/")
if(!dir.exists(deg_design_oudir)) {
    cat("Creating", deg_design_oudir, "\n")
    dir.create(deg_design_oudir)
} else {
    cat(deg_design_oudir, "already exists.\n")
}

# Create dataframe to store per-cell type DGE results.
all_dge <- data.frame(matrix(nrow = 0, ncol = 0))

# Iterate over cell types present in metadata table based on user-defined cell type column ID...
loop <- 1
for(cell_type in unique(metadata_table[,cell_type_column_id])) {
    print(
        paste0("Processing ", cell_type, " (", loop, " out of ", length(unique(metadata_table[,cell_type_column_id])), ")...")
    )

    # If for a given cell type there are not enough pseudo-bulk samples, skip that cell type.
    n_psamples_cell_type <- nrow(metadata_table[metadata_table[,cell_type_column_id] == cell_type,])
    if(n_psamples_cell_type < 2) {
        print(paste0(cell_type, " do not have enough pseudo-bulk samples (n=", n_psamples_cell_type, "). Skipping..."))
        loop <- loop + 1
        next
    }

    # Subset count matrix to keep only pseudo-bulk samples that correspond to a given cell type.
    cdata_subset <- count_matrix[,colnames(count_matrix) %in% rownames(metadata_table)[metadata_table[,cell_type_column_id] == cell_type]]

    # Subset metadata to keep only pseudo-bulk samples that correspond to a given cell type.
    metadata_subset <- metadata_table[metadata_table[,cell_type_column_id] == cell_type,]

    # If for a given cell type both control and condition of interest have less than 2 pseudo-bulk samples,
    # skip that cell type - DESeq2 will not do an 1-vs-1 comparison!
    n_psamples_control <- nrow(metadata_subset[metadata_subset[,column_to_compare] == control_condition,])
    n_psamples_interest <- nrow(metadata_subset[metadata_subset[,column_to_compare] == condition_of_interest,])
    if(n_psamples_control < 2 & n_psamples_interest < 2) {
        print(paste0("For ", cell_type, ", neither of the conditions has enough pseudo-bulk samples to be compared. Skipping..."))
        loop <- loop + 1
        next
    }

    # Perform DGE in a given cell type.
    cell_type_dge <- perform_DGE(
        cdata = cdata_subset,
        mdata = metadata_subset,
        columns_for_design_ = columns_for_design,
        column_to_compare_ = column_to_compare,
        condition_of_interest_ = condition_of_interest,
        control_condition_ = control_condition,
        cell_type_id = cell_type,
        outdir = norm_counts_design_oudir
    )

    # Add gene names as new column to avoid formatting issues when putting all results from different cell types together.
    cell_type_dge$gene <- rownames(cell_type_dge)

    # Add cell type results to dataframe.
    all_dge <- rbind(
        all_dge, 
        cell_type_dge
    )

    # Update counter.
    loop <- loop + 1
}

# Get gene symbols if input gene IDs are in the ENSEMBL format.
if(TRUE %in% grepl("ENS", all_dge$gene)) {
    
    # Create column for ENSEMBL IDs.
    all_dge$ensembl <- all_dge$gene

    # Get correspondence between ENSEMBL IDs and gene symbols.
    correspondence <- AnnotationDbi::select(
        org.Hs.eg.db, 
        keys = all_dge$ensembl, 
        columns = "SYMBOL", 
        keytype = "ENSEMBL"
    )

    # Create column with gene symbols.
    all_dge$gene <- correspondence$SYMBOL[match(all_dge$ensembl, correspondence$ENSEMBL)]

    # Replace NAs with "Unknown" for cases in which the gene symbols are missing.
    all_dge$gene[is.na(all_dge$gene)] <- "Unknown"
}

# Write unfiltered results to output csv file.
write.csv(
    all_dge, 
    paste0(
        deg_design_oudir, output_prefix, "DESeq2_", column_to_compare, "_", condition_of_interest, "_vs_", control_condition, "_",
        cell_type_column_id, "_design_", paste(columns_for_design, collapse = "_"), "_DEGs_noFilter.csv"
    )
)

# Filter results & write them to output csv file.
write.csv(
    all_dge[all_dge$padj <= 0.05,], 
    paste0(
        deg_design_oudir, output_prefix, "DESeq2_", column_to_compare, "_", condition_of_interest, "_vs_", control_condition, "_",
        cell_type_column_id, "_design_", paste(columns_for_design, collapse = "_"), "_DEGs_padj0.05_filter.csv"
    )
)