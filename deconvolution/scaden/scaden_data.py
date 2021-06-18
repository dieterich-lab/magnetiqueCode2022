#! /usr/bin/env python3

"""Prepare data for Scaden input
"""

import sys
import os
import glob
import argparse
import logging
import subprocess

import anndata

import pandas as pd
import loompy as lp
import numpy as np
import scipy as sp

logger = logging.getLogger(__name__)


xstr = lambda s: '' if s is None else str(s)


def get_filename(sc_data, tup_str, to_append, ext):
    
    ret = f"{sc_data}"
    for t1, t2 in tup_str:
        s = f"{xstr(t1)}{xstr(t2)}"
        if s:
            ret = f"{ret}-{s}"
    t = ''
    if to_append is not None:
        t = f"_{to_append}"
    return f"{ret}{t}.{ext}"
    

def write_txt(ds, layer, sc_data, n_attributes, output_dir, 
              cluster_key, gene_key, 
              sample_str='s', sample=None, 
              library_str=None, library=None, mask=None):
    
    ext = 'txt'
    
    if mask is None:
        mask = np.array([True]*n_attributes) # h5py?
    
    att = [(sample_str, sample), (library_str, library)]
    loc = get_filename(sc_data, att, 'celltypes', ext)
    fname_celltypes = os.path.join(output_dir, loc)
    loc = get_filename(sc_data, att, 'counts', ext)
    c_celltypes = os.path.join(output_dir, loc)

    pd.DataFrame(ds.ca[mask][cluster_key]).to_csv(fname_celltypes, 
                                                  sep="\t",
                                                  header=['Celltype'],
                                                  index=False)
    
    # dense, transpose
    if layer == 'X':
        try:
            X = ds[:,mask]
        except: # TypeError: Indexing arrays must have integer dtypes - bug see e.g. https://github.com/h5py/h5py/issues/1847
            X = ds[:,:][:,mask]
    else:
        try:
            X = ds.layers[layer][:,mask]
        except: 
            X = ds.layers[layer][:,:][:,mask]
            
    df = pd.DataFrame(X)
    df.index = ds.ra[gene_key]
    df.columns = ds.ca[mask]['obs_names']
    df.T.to_csv(c_celltypes, 
                sep='\t', 
                header=True, 
                index=True, 
                float_format='%.5f')
            
    
def write_h5ad(ds, layer, sc_data, n_attributes, output_dir, 
               cluster_key, gene_key, 
               sample_str='s', sample=None, 
               library_str=None, library=None, mask=None):
    
    ext = 'h5ad'
    
    if mask is None:
        mask = np.array([True]*n_attributes) # h5py?
        
    # dense, transpose, save as sparse below
    if layer == 'X':
        try:
            X = ds[:,mask].T
        except: # TypeError: Indexing arrays must have integer dtypes - bug see e.g. https://github.com/h5py/h5py/issues/1847
            X = ds[:,:][:,mask].T
    else: 
        try:
            X = ds.layers[layer][:,mask].T
        except: 
            X = ds.layers[layer][:,:][:,mask].T
        
    obs = pd.DataFrame([ds.ca[mask]['obs_names'], ds.ca[mask][cluster_key]]).T
    obs.rename(columns={1:'Celltype'}, inplace=True)
    obs.set_index(0, drop=True, inplace=True, verify_integrity=True)
    obs.index.name = None 

    # there are in some cases duplicate genes?
    var = pd.DataFrame(ds.ra[gene_key])
    try:
        var.set_index(0, drop=True, inplace=True, verify_integrity=True)
    except:
        m = var.duplicated()
        logger.warning(f"{sum(m)} duplicated genes removed from {sc_data}, sample {sample}")
        var = var[~m]
        X = X[:,~m]
        var.set_index(0, drop=True, inplace=True, verify_integrity=True)
    var.index.name = None 

    # generate AnnData
    adata = anndata.AnnData(sp.sparse.csr_matrix(X),
                            obs=obs,
                            var=var)
    
    att = [(sample_str, sample), (library_str, library)]
    loc = get_filename(sc_data, att, None, ext)

    adata.write(os.path.join(output_dir, loc))
        
        
def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Write Scaden input files.""")

    parser.add_argument('sc_data')
    parser.add_argument('fname')
    parser.add_argument('output_dir')
    
    parser.add_argument('--samples', nargs="+", required=True) # ignored if [--no-split]
    
    parser.add_argument('--fmt', help="Output format", choices=['h5ad', 'txt'], default='txt')
    parser.add_argument('--layer', help="""Which layer to use for the matrix. If X, then use data matrix. Note: although
                        layers are stored chunked and block-compressed, it seems that X is assigned from dense by loompy
                        when slicing, while we need to convert layers. The current behaviour might be unexpected 
                        if things change.""", default='X')
    parser.add_argument('--no-split', help="Files are no split per sample and/or library.", action='store_true')
    parser.add_argument('--split-key-sample', help="Sample key, used unless [--no-split] is passed.", default='Individual')
    parser.add_argument('--split-key-library', help="""Library key, used unless [--no-split] is passed, and if the key
                        is found in attributes. In that case, make sure it is not meant for something else.""", default='Library')
    
    parser.add_argument('--cluster-key', help="""We use standard nomenclature""", default='ClusterID')
    parser.add_argument('--gene-key', help="""We use standard nomenclature""", default='Accession')

    args = parser.parse_args()
    
    # don't use a context manager here...
    ds = lp.connect(args.fname)
    n_cells = ds.shape[1]
    
    if args.no_split:
        if args.fmt == 'txt':
            write_txt(ds, args.layer, args.sc_data, n_cells, args.output_dir, 
                      args.cluster_key, args.gene_key, sample_str=None)
        else:
            write_h5ad(ds, args.layer, args.sc_data, n_cells, args.output_dir, 
                       args.cluster_key, args.gene_key, sample_str=None)
    else:
        for sample in args.samples:
            mask = ds.ca[args.split_key_sample] == sample
            if args.split_key_library in ds.ca.keys():
                libs = np.unique(ds.ca[mask][args.split_key_library]) 
                for lib in libs:
                    mask = mask & (ds.ca[args.split_key_library] == lib)
                    if args.fmt == 'txt':
                        write_txt(ds, args.layer, args.sc_data, n_cells, args.output_dir, 
                                  args.cluster_key, args.gene_key, 
                                  sample_str='s', sample=sample,
                                  library_str=None, library=lib, mask=mask)
                    else:
                        write_h5ad(ds, args.layer, args.sc_data, n_cells, args.output_dir, 
                                   args.cluster_key, args.gene_key, 
                                   sample_str='s', sample=sample,
                                   library_str=None, library=lib, mask=mask)
            else:
                if args.fmt == 'txt':
                    write_txt(ds, args.layer, args.sc_data, n_cells, args.output_dir, 
                              args.cluster_key, args.gene_key, 
                              sample_str='s', sample=sample, mask=mask)
                else:
                    write_h5ad(ds, args.layer, args.sc_data, n_cells, args.output_dir, 
                               args.cluster_key, args.gene_key, 
                               sample_str='s', sample=sample, mask=mask)

    ds.close()

        
if __name__ == '__main__':
    main()
