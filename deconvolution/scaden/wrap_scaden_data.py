#! /usr/bin/env python3

"""Prepare data for Scaden input loom -> h5ad
"""

import sys
import os
import glob
import argparse
import logging

import loompy as lp

import utils


logger = logging.getLogger(__name__)

executable = os.path.join(os.getcwd(), 'scaden_data.py')

    

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Wrapper. Write Scaden input files loom -> h5ad, 
                                     allowing to split by sample and library (sc, sn, etc.).""")

    parser.add_argument('input_dir')
    parser.add_argument('output_dir')
    
    parser.add_argument('--fmt', help="Output format", choices=['h5ad', 'txt'], default='txt')
    parser.add_argument('--layer', help="""Which layer to use for the matrix. If X, then use data matrix. Note: although
                        layers are stored chunked and block-compressed, it seems that X is assigned from dense by loompy
                        when slicing, while we need to convert layers. The current behaviour might be unexpected 
                        if things change.""", default='X')
    parser.add_argument('--no-split', help="Files are no split per sample and/or library.", action='store_true')
    parser.add_argument('--split-key-sample', help="Sample key, used unless [--no-split] is passed.", default='Individual')
    parser.add_argument('--split-key-library', help="""Library key, used unless [--no-split] is passed, and if the key
                        is found in attributes. In that case, make sure it is not meant for something else.""", default='Library')
    parser.add_argument('--pattern', help="Input format file pattern. Currently only loom is supported.", default='loom')
    
    
    utils._add_logging_options(parser)
    utils._add_sbatch_options(parser)
    args = parser.parse_args()
    utils._update_logging(args)
    
    
    if args.pattern != 'loom':
        logger.critical(f"{args.pattern} currently not supported. Use loom.")
        return
    
    try:
        output_dir = os.path.join(args.output_dir, 'input')
        os.makedirs(output_dir)
    except FileExistsError:
        # directory already exists
        logger.warning(f"{output_dir} already exists. Files are overwritten by default.")

    samples = {}
    filen = {}
    for fname in glob.glob(os.path.join(args.input_dir, f"*.{args.pattern}")):
        sc_data = os.path.basename(fname).split('_')[0] # default underscore with dataset name in first place
        filen[sc_data] = fname
        with lp.connect(fname) as ds:
            samples[sc_data] = list(set(ds.ca[args.split_key_sample]))

    all_args = (f"--fmt {args.fmt} --layer {args.layer} "
                f"--split-key-sample {args.split_key_sample} --split-key-library {args.split_key_library}")
    if args.no_split:
        all_args = f"--no-split {all_args}"
    
    # pass sample list anyway, [--no-split] overrides if set
    for sc_data, sample_list in samples.items():
        cmd = "{} {} {} {} --samples {} {}".format(
            executable,
            sc_data,
            filen[sc_data],
            output_dir,
            ' '.join(sample_list),
            all_args)
        utils._check_sbatch(cmd, args=args)
        

        
if __name__ == '__main__':
    main()
