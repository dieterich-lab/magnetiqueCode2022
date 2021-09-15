#! /usr/bin/env python3

"""Group k-fold cross-validation of MAGNet data
"""

import sys
import os
import yaml
import argparse
import logging

import pandas as pd

import utils

logger = logging.getLogger(__name__)


executable = os.path.join(os.getcwd(), 'scaden_cross_validation.py')
        

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Wrapper for scaden_cross_validation.py""")
    
    parser.add_argument('config', help="""Configuration file.""")
    
    utils._add_logging_options(parser)
    utils._add_sbatch_options(parser)
    args = parser.parse_args()
    utils._update_logging(args)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)
    
    # check that we are using the right environment for scaden
    utils._check_programs_exist(['scaden'])
    
    call = not args.do_not_call
    do_not_call_str = ""
    if not call:
        do_not_call_str = "--do-not-call"
        
    try:
        os.makedirs(config['output_dir'])
    except FileExistsError:
        # directory already exists
        logger.warning(f"{config['output_dir']} already exists, will not continue.")
        return
    
    try:
        folds = pd.read_csv(config['folds_file'], sep='\t')
    except FileNotFoundError:
        logger.warning(f"<folds_file> No such file: {config['folds_file']}.")
        return
    
    k = folds.shape[1]
    all_folds = [str(f) for f in list(range(k))]
    logger.info(f"Performing {k}-fold cross-validation.")
    
    # iterate through folds
    for (fold, test_data) in folds.iteritems():
        cmd = "{} {} {} {}".format(
            executable,
            args.config,
            fold,
            do_not_call_str)
        utils._check_sbatch(cmd, args=args)
        

if __name__ == '__main__':
    main()
