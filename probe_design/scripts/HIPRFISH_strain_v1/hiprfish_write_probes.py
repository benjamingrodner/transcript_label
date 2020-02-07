import os
import re
import glob
import argparse
import subprocess
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

def write_design_level_target_probes(probes, probes_dir):
    design_level = probes.design_level.values[0]
    design_target = probes.design_target.values[0]
    design_level_target_dir = '{}/{}/{}'.format(probes_dir, design_level, design_target)
    if not os.path.exists(design_level_target_dir):
        os.makedirs(design_level_target_dir)
    for i in range(probes.shape[0]):
        probe_seq = SeqRecord(Seq(probes.seq.values[i]).reverse_complement(), id = '{}_{}_{}'.format(design_level, design_target, probes.probe_id.values[i]), description = '')
        probe_fasta_filename = '{}/{}/{}/{}.fasta'.format(probes_dir, design_level, design_target, probes.index[i])
        SeqIO.write(probe_seq, probe_fasta_filename, 'fasta')
    return(probes.shape)

def write_probes(probe_summary_filename, probes_dir):
    cluster = LocalCluster(n_workers = 100, threads_per_worker = 1)
    client = Client(cluster)
    probes = dd.read_hdf(probe_summary_filename, '*/*', mode = 'r').compute()
    client.close()
    probes_write = probes.groupby(['design_level', 'design_target']).apply(write_design_level_target_probes, probes_dir)
    # for target_rank in taxonomic_levels:
    #     design_level_dir = '{}/{}'.format(probes_dir, target_rank)
    #     if not os.path.exists(design_level_dir):
    #         os.makedirs(design_level_dir)
    #     probes = dd.read_hdf(probe_summary_filename, '{}/*'.format(target_rank)).compute()
    #     probes_write = probes.groupby(['design_level', 'design_target']).apply(write_design_level_target_probes, probes_dir)
    #     # probes = dd.read_hdf(probe_summary_filename, '{}/*'.format(target_rank)).compute()
    #     # for i in range(probes.shape[0]):
    #     #     design_level = probes.design_level.values[i]
    #     #     design_target = probes.design_target.values[i]
    #     #     probe_id = probes.index[i]
    #     #     probe_seq = SeqRecord(Seq(probes.seq.values[i]).reverse_complement(), id = '{}_{}_{}'.format(design_level, design_target, probe_id), description = '')
    #     #     probe_fasta_filename = '{}/{}/{}/{}.fasta'.format(probes_dir, design_level, design_target, probe_id)
    #     #     SeqIO.write(probe_seq, probe_fasta_filename, 'fasta')
    return

###############################################################################################################
# main function
###############################################################################################################

def main():
    parser = argparse.ArgumentParser('Design FISH probes for a complex microbial community')
    parser.add_argument('probe_summary_filename', type = str, help = 'Input FASTA file containing 16S sequences')
    args = parser.parse_args()
    primer3_dir = os.path.split(args.probe_summary_filename)[0]
    data_dir = os.path.split(primer3_dir)[0]
    probes_dir = '{}/probes'.format(data_dir)
    evaluate_dir = '{}/evaluate'.format(data_dir)
    if not os.path.exists(probes_dir):
        os.makedirs(probes_dir)
    if not os.path.exists(evaluate_dir):
        os.makedirs(evaluate_dir)
    write_probes(args.probe_summary_filename, probes_dir)
    return

if __name__ == '__main__':
    main()
