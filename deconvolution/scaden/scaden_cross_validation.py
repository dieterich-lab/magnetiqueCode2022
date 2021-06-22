#! /usr/bin/env python3

"""Group k-fold cross-validation of MAGNet data
"""

import sys
import os
import yaml
import argparse
import logging

import itertools

import anndata
import pandas as pd

sys.path.append("/prj/MAGE/analysis/deconvolution/scaden")
import utils

logger = logging.getLogger(__name__)


batch_script = os.path.join(os.getcwd(), 'run-scaden-gpus')


def simulate(parent, subdir_name, lfiles, prefix, config, call):
    
    subdir = os.path.join(parent, subdir_name)
    subdir_input = os.path.join(subdir, 'input')
    os.makedirs(subdir_input, exist_ok=True)
    for src in lfiles:
        dst = os.path.join(subdir_input, os.path.basename(src))
        if os.path.lexists(dst):
            os.remove(dst)
        os.symlink(src, dst)
    
    cmd = (f"scaden simulate --cells {config['simulate']['cells_per_sample']} "
            f"--n_samples {config['simulate']['n_samples']} --data {subdir_input} "
            f"--pattern \"{config['pattern']}\" --data-format {config['fmt']} "
            f"--prefix {prefix} --out {subdir}")
    for u in config['simulate']['unknown_cell_type']:
        cmd = f"{cmd} --unknown {u}"
    return utils._check_srun(cmd, call=call)


def process(bulk_file, parent, prefix, config, call):
    
    # currently Scaden only inputs bulk (prediction data) as txt ... 
    # we should modify that as well... 
    adata = anndata.read_h5ad(f"{bulk_file}.{config['fmt']}")
    bulk_txt = f"{bulk_file}.txt"
    obs_names = anndata.AnnData(var=pd.DataFrame(index=adata.obs['ds'].values))
    obs_names.var_names_make_unique()
    X = adata.X.T # see to_df method of anndata...
    pd.DataFrame(X, 
                 index=adata.var_names.values, 
                 columns=obs_names.var_names.values).to_csv(bulk_txt, 
                                                            index=True, 
                                                            header=True, 
                                                            sep='\t',
                                                            float_format='%.10f')
    
    training_data = f"{prefix}.{config['fmt']}"
    processed_data = f"{prefix}_processed.{config['fmt']}"
    cmd = (f"scaden process {training_data} {bulk_txt} "
            f"--var_cutoff {config['process']['var_cutoff']} --processed_path {processed_data}")
    return utils._check_srun(cmd, call=call), processed_data
        

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Cross-validation using predefined folds.""")
    
    parser.add_argument('config', help="""Configuration file.""")
    parser.add_argument('fold', help="""Fold number.""", type=str)

    utils._add_logging_options(parser)
    utils._add_sbatch_options(parser)
    args = parser.parse_args()
    utils._update_logging(args)
    
    config = yaml.load(open(args.config), Loader=yaml.FullLoader)
    
    try:
        folds = pd.read_csv(config['folds_file'], sep='\t')
    except FileNotFoundError:
        logger.warning(f"<folds_file> No such file: {config['folds_file']}.")
        return
    k = folds.shape[1]
    all_folds = [str(f) for f in list(range(k))]
    
    call = not args.do_not_call
    
    # process files
    test_files = folds[args.fold].values
    test_files = [d for d in test_files if d==d] # remove any NaNs due to unequal fold length
    train_data = folds[list(filter(lambda x: x != args.fold, all_folds))].values
    train_files = list(itertools.chain.from_iterable(train_data))
    train_files = [d for d in train_files if d==d] 

    # create local directory for input files, test and trainign data, model, etc.
    fold_dir = os.path.join(config['output_dir'], f"fold{args.fold}")
    os.makedirs(fold_dir, exist_ok=True)
    
    # simulate bulk using test data
    prefix = (f"test_bulk_{config['simulate']['cells_per_sample']}cells_"
              f"{config['simulate']['n_samples']}samples")
    # full path, no extension
    bulk_file = os.path.join(fold_dir, 'bulk', prefix)
    _ = simulate(fold_dir, 'bulk', test_files, bulk_file, config, call)
    
    # write cell composition as txt
    adata = anndata.read_h5ad(f"{bulk_file}.{config['fmt']}")
    true_types = f"{bulk_file}_cell_composition.txt"
    adata.obs[adata.uns['cell_types']].to_csv(true_types, 
                                              header=True, 
                                              index=False, 
                                              sep='\t', 
                                              float_format='%.10f')
    
    # training: simulate, process, train, predict
    prefix = (f"train_bulk_{config['simulate']['cells_per_sample']}cells_"
              f"{config['simulate']['n_samples']}samples")
    bulk_training = os.path.join(fold_dir, 'train', prefix)
    _ = simulate(fold_dir, 'train', train_files, bulk_training, config, call)
    _, processed_data = process(bulk_file, fold_dir, bulk_training, config, call)
    
    # now we need to wrap call to GPUs for training...
    fold_dir_model = os.path.join(fold_dir, 'model')
    export = {'TRAIN': 1,
              'MODEL_DIR': fold_dir_model,
              'BATCH_SIZE': config['train']['batch_size'],
              'LEARNING_RATE': config['train']['learning_rate'], 
              'STEPS': config['train']['steps'],
              'SEED': config['train']['seed'],
              'TRAINING_DATA': processed_data}
    job_id_train = utils._check_sbatch(batch_script,
                                       call=call,
                                       use_slurm=True,
                                       partitions='gpu',
                                       num_gpus=config['gres']['num'],
                                       gpu_name=config['gres']['name'],
                                       export=export)
    
    # ... and prediction
    prefix = (f"predicted_bulk_{config['simulate']['cells_per_sample']}cells_"
              f"{config['simulate']['n_samples']}samples_"
              f"cell_composition.txt")
    prediction_data = os.path.join(fold_dir, prefix)
    export = {'TRAIN': 0,
              'MODEL_DIR': fold_dir_model,
              'SEED': config['train']['seed'],
              'OUTNAME': prediction_data,
              'BULK': f"{bulk_file}.txt"}
    _ = utils._check_sbatch(batch_script, 
                            call=call,
                            use_slurm=True,
                            partitions='gpu',
                            num_gpus=config['gres']['num'],
                            gpu_name=config['gres']['name'],
                            export=export,
                            dependencies=[job_id_train])
       

if __name__ == '__main__':
    main()
