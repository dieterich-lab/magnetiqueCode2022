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

import scipy as sp

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Wrapper write Scaden input files.""")

    parser.add_argument('scData')
    parser.add_argument('fname')
    parser.add_argument('locDir')
    #parser.add_argument('sample')
    parser.add_argument('--samples', nargs="+", required=True)
    parser.add_argument('--fmt', choices=['h5ad', 'txt'], default='txt')
    
    args = parser.parse_args()

    ds = lp.connect(args.fname)
    for sample in args.samples:
        mask = ds.ca['Individual'] == sample # does not work ds[:,mask] on .scpy bug version
        
        if args.fmt == 'txt':
        
            loc = '{}Sample{}_celltypes.txt'.format(args.scData, sample)
            fname_celltypes = os.path.join(args.locDir, 'input', loc)
            loc = '{}Sample{}_counts.txt'.format(args.scData, sample)
            c_celltypes = os.path.join(args.locDir, 'input', loc)
        
            pd.DataFrame(ds.ca[mask]['Cluster']).to_csv(fname_celltypes, 
                                                        sep="\t",
                                                        header=['Celltype'],
                                                        index=False)
            df = pd.DataFrame(ds[:,mask])
            df.index = ds.ra['gene_ids']
            df.columns = ds.ca[mask]['obs_names']
            df.T.to_csv(c_celltypes, 
                        sep='\t', 
                        header=True, 
                        index=True, 
                        float_format='%.5f')
            
        else:
            
            # transpose, save as sparse below
            X = ds[:,mask].T
            
            obs = pd.DataFrame([ds.ca[mask]['obs_names'], ds.ca[mask]['Cluster']]).T
            obs.rename(columns={1:'Celltype'}, inplace=True)
            obs.set_index(0, drop=True, inplace=True, verify_integrity=True)
            obs.index.name = None 

            # there are in some cases duplicate genes?
            var = pd.DataFrame(ds.ra['gene_ids'])
            try:
                var.set_index(0, drop=True, inplace=True, verify_integrity=True)
            except:
                m = var.duplicated()
                logger.warning('{} duplicated genes removed from {}, sample {}'.format(sum(m), args.scData, sample))
                var = var[~m]
                X = X[:,~m]
                var.set_index(0, drop=True, inplace=True, verify_integrity=True)
            var.index.name = None 
    
            # generate AnnData
            adata = anndata.AnnData(sp.sparse.csr_matrix(X),
                                    obs=obs,
                                    var=var)
            loc = '{}Sample{}.h5ad'.format(args.scData, sample)
            adata.write(os.path.join(args.locDir, 'input', loc))

    ds.close()
            
    #with lp.connect(args.fname) as ds:
        #mask = ds.ca['Individual'] == args.sample
        #loc = '{}Sample{}_celltypes.txt'.format(args.scData, args.sample)
        #fname_celltypes = os.path.join(args.locDir, 'input', loc)
        #loc = '{}Sample{}_counts.txt'.format(args.scData, args.sample)
        #c_celltypes = os.path.join(args.locDir, 'input', loc)
        ## max batch size = num cells to avoid splitting samples
        ## slicing would be equally fine ds[:,mask], and ds.ca[mask]['Cluster']
        ## if not faster?
        #for (ix, selection, view) in ds.scan(items=mask, axis=1, batch_size=ds.shape[1]):
            ## cell types
            #pd.DataFrame(view.ca['Cluster']).to_csv(fname_celltypes, 
                                                    #sep="\t",
                                                    #header=['Celltype'],
                                                    #index=False)
            ## counts
            #df = pd.DataFrame(view[:,:])
            #df.index = view.ra['gene_ids']
            #df.columns = view.ca['obs_names']
            #df.T.to_csv(c_celltypes, 
                        #sep='\t', 
                        #header=True, 
                        #index=True, 
                        #float_format='%.5f')
        
        
if __name__ == '__main__':
    main()
