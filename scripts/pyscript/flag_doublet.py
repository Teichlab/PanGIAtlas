#!/usr/bin/env python
"""Flag doublet by scrublet

Usage: program [options] <input> <output>

Options:
  --debug           print debug information
  --profile         print profile information
  --filter <str>    only consider cells that have "True" values in this variable of obs
  --groupby <str>   run scrublet for each group specified by this variable of obs
  <input>           input h5ad
  <output>          output csv
"""


import logging
import signal
import sys
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import gc
import pandas as pd
import scanpy as sc
import sctk as sk

def main(args):
    logging.debug(args)

    input_h5ad = args["input"]
    output_csv = args["output"]
    filter_var = args["filter"]
    groupby = args["groupby"]

    ad = sc.read(input_h5ad)
    if groupby and groupby in ad.obs.columns:
        groups = ad.obs[groupby].cat.categories
        dfs = []
        for i, grp in enumerate(groups):
            print(grp)
            k_grp = ad.obs[groupby] == grp
            if filter_var and filter_var in ad.obs.columns:
                filter_value = ad.obs[filter_var]
                if filter_value.dtype.kind == "b":
                    k_filter = filter_value.values
                elif filter_value.dtype.kind == "O":
                    k_filter = filter_value.values == "True"
                else:
                    raise TypeError(f"{filter_var} is not a boolean or categorical")
                k = k_grp & k_filter
            else:
                k = k_grp
            if k.sum() > 10:
                ad1 = ad[k]
                try:
                    grp_df = sk.run_scrublet(ad1, inplace=False)
                    grp_df.to_csv(output_csv.replace(".csv", f".{grp}.csv"))
                    dfs.append(grp_df)
                except:
                    logging.warn(f"{grp}: scrublet failed")
                del ad1
                gc.collect()
        df = pd.concat(dfs)
    else:
        if filter_var and filter_var in ad.obs.columns:
            filter_value = ad.obs[filter_var]
            if filter_value.dtype.kind == "b":
                k = filter_value.values
            elif filter_value.dtype.kind == "O":
                k = filter_value.values == "True"
            else:
                raise TypeError(f"{filter_var} is not a boolean or categorical")
            ad1 = ad[k]
        else:
            ad1 = ad

        df = sk.run_scrublet(ad1, inplace=False)

    df.to_csv(output_csv)

    return 0


if __name__ == '__main__':
    from docopt import docopt
    args = docopt(__doc__)
    args = {k.lstrip('-<').rstrip('>'):args[k] for k in args}
    try:
        if args.get('debug'):
            logLevel = logging.DEBUG
        else:
            logLevel = logging.WARN
        logging.basicConfig(
            level=logLevel,
            format='%(asctime)s; %(levelname)s; %(funcName)s; %(message)s',
            datefmt='%y-%m-%d %H:%M:%S')
        if args.get('profile'):
            import cProfile
            cProfile.run('main(args)')
        else:
            main(args)
    except KeyboardInterrupt:
        logging.warning('Interrupted')
        sys.exit(1)
