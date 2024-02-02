################################
#### Load required packages ####
################################

#Script written by Jacqueline Boccacino

import scanpy as sc
import pandas as pd
import decoupler as dc
import argparse


###############################
#### Define script's main  ####
###############################

def main():

    # Parse arguments.
    parser = argparse.ArgumentParser()

    # Input h5ad file.
    parser.add_argument(
        '-i', 
        '--input-file', 
        '--input', 
        required = True,
        dest = 'input_h5ad',
        help = 'Path to an input h5ad file.'
    )
    # Output h5ad file.
    parser.add_argument(
        '-o', 
        '--output-prefix'
        '--output', 
        required = True,
        dest = 'output_prefix', 
        help = 'Path + prefix to output files.'
    )
    # Minimum number of cells per donor for a specified cell type (default = 10).
    parser.add_argument(
        '-m', 
        '--min-number-cells',
        required = False,
        default = 10,
        dest = 'min_number_cells', 
        help = 'Minimum number of cells per donor for a specified cell type (default = 10).'
    )
    # Minimum number of counts per donor for a specified cell type (default = 0).
    parser.add_argument(
        '-n', 
        '--min-number-counts',
        required = False,
        default = 0,
        dest = 'min_number_counts', 
        help = 'Minimum number of counts per donor for a specified cell type (default = 0).'
    )
    # Column name in AnnData.obs that corresponds to donor IDs.
    parser.add_argument(
        '-d', 
        '--donor-column-id',
        required = True,
        dest = 'donor_column_id', 
        help = 'Column name in AnnData.obs that corresponds to donor IDs.'
    )
    # Column name in AnnData.obs that corresponds to cell types.
    parser.add_argument(
        '-t', 
        '--cell-type-column-id',
        required = True,
        dest = 'cell_type_column_id', 
        help = 'Column name in AnnData.obs that corresponds to cell types/states/identities.'
    )
    # Comma-separated list of metadata columns that should be kept after pseudo-bulking.
    parser.add_argument(
        '-k', 
        '--keep-columns',
        required = False,
        dest = 'keep_columns', 
        help = 'Comma-separated list of metadata columns that should be kept after pseudo-bulking.'
    )
    group = parser.add_mutually_exclusive_group()
    # Determine that the 'raw' slot exists (e.g. adata.raw/adata.raw.X)
    group.add_argument(
        '--raw-slot', 
        action = 'store_true',
        dest = 'raw_slot', 
        default = True, 
        required = False,
        help = 'This flag tells the program that adata.raw exists [default].'
    )
    # Determine that downsampling should not be done.
    group.add_argument(
        '--no-raw-slot', 
        action = 'store_false',
        dest = 'raw_slot', 
        default = False, 
        required = False,
        help = 'This flag tells the program that adata.raw does not exist.'
    )

    # Get arguments.
    arguments = parser.parse_args()

    # Create a list with columns that should be used in DESeq2's analysis design.
    keep_columns = arguments.keep_columns.split(',')
    print(keep_columns)

    # Read inout h5ad file.
    adata = sc.read_h5ad(arguments.input_h5ad)

    # Add 'counts' and 'normalised' layers to AnnData object with raw and normalised counts, respectively.
    if arguments.raw_slot == True:
        adata.layers['counts'] = adata.raw.X
        adata.layers['normalised'] = adata.X
    else:
        adata.layers['counts'] = adata.X

    # Add donor-cell type column to adata.obs.
    adata.obs['donor_celltype'] = adata.obs[arguments.donor_column_id].astype(str) + '_' + adata.obs[arguments.cell_type_column_id].astype(str)

    # Perform pseudo-bulking with the sum mode.
    pdata = dc.get_pseudobulk(
        adata,
        sample_col = arguments.donor_column_id,
        groups_col = arguments.cell_type_column_id,
        layer = 'counts',
        mode = 'sum',
        min_cells = int(arguments.min_number_cells),
        min_counts = int(arguments.min_number_counts)
    )

    # Make sure that all columns in pdata.obs have string type.
    for column in pdata.obs.columns:
        pdata.obs[column] = pdata.obs[column].astype(str)

    # Add donor-cell type column to pdata.obs.
    pdata.obs['donor_celltype'] = pdata.obs[arguments.donor_column_id].astype(str) + '_' + pdata.obs[arguments.cell_type_column_id].astype(str)

    # Filter original adata to keep only donor-cell type combinations that are present in pdata.
    adata = adata[adata.obs['donor_celltype'].isin(pdata.obs['donor_celltype'])]

    # Get dataframe from AnnData's obs slot (before pseudo-bulking).
    metadata_before_psbulk = pd.DataFrame(
        adata.obs,
        index = adata.obs.index,
        columns = adata.obs.columns
    )
    # The above dataframe will have all cells that correspond to the pseudo-bulk samples, and we want to keep
    # only one line per donor-cell type combination.
    metadata_before_psbulk = metadata_before_psbulk.drop_duplicates(subset = 'donor_celltype', keep = "last")

    # Get dataframe from AnnData's obs slot (after pseudo-bulking).
    metadata_after_psbulk = pd.DataFrame(
        pdata.obs,
        index = pdata.obs.index,
        columns = pdata.obs.columns
    )

    # Add columns to pdata's metadata.
    for column in keep_columns:
        # Only add column if it is not already present!
        if column not in metadata_after_psbulk.columns:
            metadata_after_psbulk = pd.merge(
                metadata_after_psbulk,
                metadata_before_psbulk[['donor_celltype', column]],
                left_on = 'donor_celltype', 
                right_on= 'donor_celltype',
                how = 'left'
            )

    # Set pseudo-bulk sample names as row names.
    metadata_after_psbulk = metadata_after_psbulk.set_index('donor_celltype')

    # Save count matrix to output csv file.
    count_matrix_file = arguments.output_prefix + '_count_matrix.csv'
    pd.DataFrame(data = pdata.X, index = pdata.obs_names, columns = pdata.var_names).to_csv(count_matrix_file)

    # Save metadata to output csv file.
    metadata_file = arguments.output_prefix + '_metadata.csv'
    metadata_after_psbulk.to_csv(metadata_file)

if __name__ == "__main__":
    main()

