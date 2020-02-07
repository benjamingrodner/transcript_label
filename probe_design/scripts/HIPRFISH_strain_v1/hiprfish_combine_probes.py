import os
import re
import glob
import argparse
import subprocess
import numpy as np
import pandas as pd
import dask
import dask.dataframe as dd
from dask.distributed import Client, LocalCluster
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, generic_dna
from matplotlib import pyplot as plt
from Bio.Blast.Applications import NcbiblastnCommandline

os.environ['OMP_NUM_THREADS'] = "16"
os.environ['MKL_NUM_THREADS'] = "16"
os.environ['DASK_NUM_THREADS'] = "16"

###############################################################################################################
# HiPR-FISH-strain : combine probes
###############################################################################################################

def cm_to_inches(x):
    return(x/2.54)

def initialize_client():
    cluster = LocalCluster(n_workers = 100, threads_per_worker = 1)
    client = Client(cluster)
    return(client)

def get_molecule_id(x):
    mid = re.sub('_probes.csv', '', os.path.basename(x))
    return(mid)

def calculate_source(df):
    source_count = df.iloc[:,1:].nunique()
    source_unique = source_count > 1
    max_design_level_numeric = 7 - source_unique.sum()
    design_level = source_count.index[max_design_level_numeric]
    design_target = df[design_level].values[0]
    source_count.iloc[:max_design_level_numeric + 1] = df.iloc[0,1:max_design_level_numeric + 2].astype(str)
    source_count.iloc[max_design_level_numeric + 1:] = 'NA'
    source_count['max_design_level_numeric'] = max_design_level_numeric
    source_count['max_design_level'] = design_level
    source_count['max_design_target'] = design_target
    return(source_count)

def write_to_hdf(df, filename):
    key = '{}/{}_{}'.format(df.design_level.values[0], df.design_level.values[0], df.design_target.values[0])
    df['probe_id'] = np.arange(df.shape[0])
    df.to_hdf(filename, key, format = 'table')
    return(pd.Series('complete'))

def plot_probe_specificity_scatter(probes_summary, theme_color, filenamem):
    taxonomic_levels = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    fig = plt.figure()
    fig.set_size_inches(cm_to_inches(5), cm_to_inches(5))
    for i in range(len(taxonomic_levels)):
        plt.plot(i + 0.1*np.random.random(probes_summary.shape[0]), probes_summary.loc[:, taxonomic_levels[i]].values, 'o', markersize = 1, alpha = 0.5, color = (0,0.5,1))
    plt.xticks(np.arange(6), taxonomic_levels, rotation = 90, fontsize = 8, color = theme_color)
    plt.ylabel('Specificy', fontsize = 8, color = theme_color)
    plt.tick_params(direction = 'in', length = 2, colors = theme_color)
    plt.subplots_adjust(left = 0.25, bottom = 0.25, right = 0.98, top = 0.98)
    plt.yscale('log')
    plt.savefig(filename, dpi = 300, transparent = True)
    return

def combine_probes(primer3_dir, util_dir, probes_summary_dir):
    cluster = LocalCluster(n_workers = 100, threads_per_worker = 1)
    client = Client(cluster)
    # primer3_dir = '/workdir/hs673/Runs/V1/Samples/HIPRFISH_4/11_27_2018/primer3'
    # util_dir = '/workdir/hs673/Runs/V1/Samples/HIPRFISH_4/11_27_2018/utilities'
    probes_filenames = '{}/*_probes.csv'.format(primer3_dir)
    blast_lineage = pd.read_csv('{}/blast_lineage.tab'.format(util_dir), sep = '\t')
    taxonomic_levels = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    blast_lineage_slim = blast_lineage.loc[:, ['molecule_id'] + taxonomic_levels]
    probes = dd.read_csv(probes_filenames)
    probes['molecule_id'] = probes.source.apply(get_molecule_id, meta = ('str'))
    probes = probes.merge(blast_lineage_slim, on = 'molecule_id', how = 'left')
    probes['superkingdom'] = 2
    extended_taxonomic_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'molecule_id']
    probes_taxa = probes.loc[:, ['seq'] + extended_taxonomic_levels]
    probes_summary = probes_taxa.groupby('seq').apply(calculate_source, meta = [('superkingdom', 'int'), ('phylum','int'),('class','int'),('order','int'),('family','int'), ('genus','int'), ('species','int'), ('molecule_id', 'str'), ('max_design_level_numeric', 'int'), ('max_design_level', 'str'), ('max_design_target', 'int')])
    probes_properties = probes.loc[:,['seq', 'length', 'Tm', 'GC', 'N', 'self_any_th', 'self_end_th', 'hair-pin', 'quality']]
    probes_summary = probes_summary.compute()
    probes_summary = probes_summary.reset_index()
    probes_properties = probes_properties.drop_duplicates().compute()
    probes_summary = probes_summary.merge(probes_properties, on = 'seq', how = 'left')
    client.close()
    probe_summary_filename = '{}/probes_summary.h5'.format(probes_summary_dir)
    probes_summary.loc[:,'max_design_target'] = probes_summary.max_design_target.astype(str)
    taxonomic_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'molecule_id']
    for i in range(8):
        probes_summary_working_design_level = probes_summary.loc[probes_summary.max_design_level_numeric >= i, :]
        probes_summary_working_design_level.loc[:,'design_level'] = taxonomic_levels[i]
        probes_summary_working_design_level.loc[:,'design_target'] = probes_summary_working_design_level.loc[:, taxonomic_levels[i]]
        probes_summary_working_design_level.groupby(['design_level', 'design_target']).apply(write_to_hdf, probe_summary_filename)
    return

###############################################################################################################
# main function
###############################################################################################################


def main():
    parser = argparse.ArgumentParser('Design FISH probes for a complex microbial community')
    parser.add_argument('primer3_dir', type = str, help = 'Input FASTA file containing 16S sequences')
    args = parser.parse_args()
    data_dir = os.path.split(args.primer3_dir)[0]
    util_dir = '{}/utilities'.format(data_dir)
    probes_summary_dir = '{}/probes_summary'.format(data_dir)
    settings_filename = '{}/primer3_settings.txt'.format(args.primer3_dir)
    try:
        os.mkdir(probes_summary_dir)
    except:
        pass
    combine_probes(args.primer3_dir, util_dir, probes_summary_dir)
    return

if __name__ == '__main__':
    main()
