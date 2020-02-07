import os
import re
import glob
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Alphabet import IUPAC, generic_dna
os.environ['OMP_NUM_THREADS'] = '16'
os.environ['MKL_NUM_THREADS'] = '16'
os.environ['DASK_NUM_THREADS'] = '16'

###############################################################################################################
# HiPR-FISH Probe Design Pipeline
###############################################################################################################

###############################################################################################################
# Workflow functions
###############################################################################################################

def write_to_hdf(df, evaluation_dir):
    design_level = df.design_level.values[0]
    design_target = df.design_target.values[0]
    probe_nid = df.probe_nid.values[0]
    filename = '{}/{}/{}_probe_evaluation.h5'.format(evaluation_dir, design_level, design_target)
    key = '{}_{}/{}_{}_{}'.format(design_level, design_target, design_level, design_target, probe_nid)
    df.to_hdf(filename, key, format = 'table')
    return

def calculate_mch(df):
    qseq = df.qseq
    sseq = df.sseq
    if qseq != sseq:
        snp_indices = np.where(np.array(list(qseq)) != np.array(list(sseq)))[0]
        diffs = np.diff(snp_indices)
        mch = np.max(np.append(diffs,[snp_indices[0], len(qseq) - 1 - snp_indices[-1]]))
    else:
        mch = len(qseq)
    return(mch)

def sub_slash(str):
    return(re.sub('/', '_', str))

def get_design_info(probe_name):
    design_level, design_target, probe_nid = probe_name.split('_')
    return(pd.Series({'probe_id': probe_name, 'design_level': design_level, 'design_target': design_target, 'probe_nid': probe_nid}))

def get_blast_lineage_slim(blast_lineage_filename):
    blast_lineage = pd.read_table(blast_lineage_filename, dtype = {'staxids':str})
    lineage_columns = ['molecule_id', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    blast_lineage_slim = blast_lineage[lineage_columns]
    return(blast_lineage_slim)

def evaluate_probes(design_target_dir, evaluation_dir, blast_lineage_filename):
    probes_blast_filenames = glob.glob('{}/*.blast.out'.format(design_target_dir))
    blast_lineage_slim = get_blast_lineage_slim(blast_lineage_filename)
    for f in probes_blast_filenames:
        probes_blast = pd.read_table(f, header = None)
        probes_blast.columns = ['probe_id', 'molecule_id', 'pid', 'qcovhsp', 'length', 'mismatch', 'gapopen', 'probe_start', 'probe_end', 'molecule_start', 'molecule_end', 'evalue', 'bitscore', 'staxids', 'qseq', 'sseq']
        probes_design_info = probes_blast.probe_id.drop_duplicates().apply(get_design_info)
        probes_blast['mch'] = probes_blast.loc[:,['qseq', 'sseq']].apply(calculate_mch, axis = 1)
        probes_blast['molecule_id'] = probes_blast.molecule_id.apply(sub_slash)
        probes_blast = probes_blast.merge(blast_lineage_slim, on = 'molecule_id', how = 'left')
        probes_blast = probes_blast.merge(probes_design_info, on = 'probe_id', how = 'left')
        write_to_hdf(probes_blast, evaluation_dir)
    return

###############################################################################################################
# main function
###############################################################################################################


def main():
    parser = argparse.ArgumentParser('Design FISH probes for a complex microbial community')
    # input blast filename
    parser.add_argument('probes_blast_complete_filename', type = str, help = 'Input file containing blast results')
    parser.add_argument('-d', '--design_level', dest = 'design_level', type = str)
    args = parser.parse_args()

    design_level_dir, probe_blast_complete_basename = os.path.split(args.probes_blast_complete_filename)
    design_target = re.sub('_blast_complete.txt', '', probe_blast_complete_basename)
    design_target_dir = '{}/{}'.format(design_level_dir, design_target)
    probes_dir = os.path.split(design_level_dir)[0]
    data_dir = os.path.split(probes_dir)[0]
    probe_evaluation_filename = re.sub('_blast_complete.txt', '_probe_evaluation.h5', args.probes_blast_complete_filename)
    probe_evaluation_complete_filename = '{}/evaluate/{}/{}_probe_evaluation_complete.txt'.format(data_dir, args.design_level, design_target)
    evaluation_dir = '{}/evaluate'.format(data_dir)
    if not os.path.exists(evaluation_dir):
        os.makedirs(evaluation_dir)
    blast_lineage_filename = '{}/utilities/blast_lineage.tab'.format(data_dir)
    if not os.path.exists(probe_evaluation_complete_filename):
        evaluate_probes(design_target_dir, evaluation_dir, blast_lineage_filename)
    file = open(probe_evaluation_complete_filename, 'w')
    file.write('Probe evaluation level is complete.')
    file.close()
    return

if __name__ == '__main__':
    main()
