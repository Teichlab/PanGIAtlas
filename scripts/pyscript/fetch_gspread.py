#!/usr/bin/env python3
"""Fetch google spreadsheet

Usage: fetch_gspread.py [options] <url> <output>

Options:
  --debug            print debug information
  --profile          print profile information
  -s, --sheet <str>  work sheet name [default: Sheet1]
  <url>              google spreadsheet url
  <output>           output csv
"""


import logging
import signal
import sys
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

import os
import time
from datetime import datetime
import pandas as pd
import gspread as gs


def main(args):
    logging.debug(args)

    url = args["url"]
    sheet = args["sheet"]
    output = args["output"]

    gsc = gs.service_account()
    sh = gsc.open_by_url(url)
    ws = sh.worksheet(sheet)

    df = pd.DataFrame(ws.get_all_records())

    if os.path.exists(output):
        logging.warning(f"{output} exists, overwriting in 3 seconds")
        for i in range(3):
            time.sleep(1)

    outdir = os.path.dirname(output)
    os.makedirs(outdir, exist_ok=True)
    df.to_csv(output)
    logging.debug("fetch done")

    with open(os.path.join(outdir, "fetch_metadata.log"), "a+") as f:
        print(f"{output} fetched on {datetime.now()}", file=f)
    logging.debug("log done")

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
