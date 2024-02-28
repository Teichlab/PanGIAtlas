import re
import os
import pandas as pd
import decoupler as dc
import omnipath as op
import liana as li

def get_cellpairs(W, fct):
    # reset index and separate by & to source and target, keep fct
    cell_pairs = W.reset_index().rename(columns={'ct':'source_target'}).copy()
    cell_pairs[['source', 'target']] = cell_pairs['source_target'].str.split('&', expand=True)
    cell_pairs = cell_pairs[['source', 'target', fct]].copy()
    
    return cell_pairs

def alphanumeric_sort_key(s):
    """
    Convert the input string into a list of integers and strings for sorting.
    """
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]


def enrich_and_format(dea_df, cell_group, groupby, stat, net):
    
    ct_df = dea_df[dea_df[groupby] == cell_group].copy()
    mat = ct_df[[stat]].T.rename(index={stat: cell_group})
        
    # Infer pathway activities with ulm
    tf_acts, tf_pvals = dc.run_ulm(mat=mat, net=net)
    tf_acts = tf_acts.T.rename(columns={cell_group: 'activity'})
    tf_pvals = tf_pvals.T.rename(columns={cell_group: 'pvalue'})

    act = tf_acts.merge(tf_pvals, left_index=True, right_index=True)
    act['fdr'] = dc.p_adjust_fdr(act['pvalue'])
    return act


def check_file_exists(func):
    def wrapper(*args, **kwargs):
        path = kwargs.pop('path', None)
        if path is None:
            raise ValueError("The 'path' argument is required.")

        if not os.path.exists(path):
            data = func(*args, **kwargs)
            data.to_csv(path)
            return data
        else:
            return pd.read_csv(path, index_col=0)
    return wrapper


@check_file_exists
def get_ppi(*args, **kwargs):
    ppis = op.interactions.OmniPath().get(genesymbols = True)

    ppis['mor'] = ppis['is_stimulation'].astype(int) - ppis['is_inhibition'].astype(int)
    ppis = ppis[(ppis['mor'] != 0) & (ppis['curation_effort'] >= 3) & ppis['consensus_direction']]

    ppis = ppis[['source_genesymbol', 'mor', 'target_genesymbol']]
    ppis.columns = ['source', 'mor', 'target']
    return ppis

@check_file_exists
def get_lr_progeny(*args, **kwargs):
    net = dc.get_progeny(organism='human', top=5000)
    lr_pairs = li.resource.select_resource('consensus')
    lr_progeny = li.rs.generate_lr_geneset(lr_pairs, net, lr_sep="&")
    return lr_progeny

@check_file_exists
def get_regulons(*args, **kwargs):
    return dc.get_collectri(organism='human', *args, **kwargs)