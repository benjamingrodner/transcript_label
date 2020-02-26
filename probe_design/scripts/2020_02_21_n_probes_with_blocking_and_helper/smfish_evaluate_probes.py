import argparse
import pandas as pd
import subprocess
import os
import multiprocessing
import glob
import re
import itertools
from SetCoverPy import setcover
import numpy as np
import random
from ete3 import NCBITaxa
from SetCoverPy import setcover
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, generic_dna
from Bio.Blast.Applications import NcbiblastnCommandline
from ete3 import NCBITaxa
import re
import time
import smtplib
from joblib import Parallel, delayed
os.environ['OMP_NUM_THREADS'] = '1'
###############################################################################################################
# HiPR-FISH Probe Design Pipeline
###############################################################################################################

###############################################################################################################
# Workflow functions
###############################################################################################################

def send_confirmation():
    server = smtplib.SMTP('smtp.gmail.com', 587)
    server.starttls()
    server.login('proudquartzprogramcheck@gmail.com', 'Shprogramcheck')
    server.sendmail('proudquartzprogramcheck@gmail.com', 'hs673@cornell.edu', 'Subject: {}\n\n{}'.format('Program finished', ''))

def mch_test_v4(df, blast_output_filename):
    qseq_array = df.qseq.values
    sseq_array = df.sseq.values
    mch_array = np.zeros(len(qseq_array))
    for i in range(len(qseq_array)):
        qseq = qseq_array[i]
        sseq = sseq_array[i]
        try:
            if qseq != sseq:
                snp_indices = np.where(np.array(list(qseq)) != np.array(list(sseq)))[0]
                diffs = np.diff(snp_indices)
                mch_array[i] = np.max(np.append(diffs,[snp_indices[0], len(qseq) - 1 - snp_indices[-1]]))
            else:
                mch_array[i] = len(qseq)
        except TypeError:
            num_line = subprocess.call(['wc', '-l', blast_output_filename])
            print(num_line)
            print(blast_output_filename)
            print(df.ix[i,:])
    return(mch_array)

def evaluate_feature_probes(feature_probe_directory, store):
    probe_filenames = glob.glob(feature_probe_directory + '/*.fasta')
    for filename in probe_filenames:
        blast_output_filename = filename + '.blast.out'
        try:
            chunk = pd.read_table(blast_output_filename, header = None)
            # print(chunk.shape)
            chunk.columns = ['probe_id', 'feature', 'pid', 'qcovhsp', 'length', 'mismatch', 'gapopen', 'probe_start', 'probe_end', 'molecule_start', 'molecule_end', 'evalue', 'bitscore', 'staxids', 'qseq', 'sseq']
            chunk['mch'] = mch_test_v4(chunk, blast_output_filename)
            chunk_name = re.search('probe_[0-9]+', filename).group(0)
            # print(chunk_name)
            # print(chunk_name)
            store.append(chunk_name, chunk, index=False, data_columns = True, min_itemsize = {'sseq': 30, 'qseq': 30}, complib = 'lzo', complevel = 9)
        except pd.errors.EmptyDataError:
            pass

def evaluate_probes(sim_dir, feature):
    feature_probe_directory = '{}/primer3/{}'.format(sim_dir, feature)
    probe_blast_hdf5_filename = '{}/blast/{}_probe_evaluation.h5'.format(sim_dir, feature)
    if not os.path.exists(probe_blast_hdf5_filename):
        store = pd.HDFStore(probe_blast_hdf5_filename, complib = 'lzo', complevel = 9)
        evaluate_feature_probes(feature_probe_directory, store)
        store.close()
        # for blastfile in glob.glob('{}/*.blast.out'.format(feature_probe_directory)):
        #     os.remove(blastfile)
    return

###############################################################################################################
# main function
###############################################################################################################


def main():
    parser = argparse.ArgumentParser('Design FISH probes for a complex microbial community')
    # input blast filename
    parser.add_argument('design_probe_filename', type = str, help = 'Input file containing blast results')

    args = parser.parse_args()

    sim_primer3_dir, feature_probe_filename = os.path.split(args.design_probe_filename)
    sim_dir = os.path.split(sim_primer3_dir)[0]
    feature = re.sub('.int', '', feature_probe_filename)
    print('Evaluating probes for {}'.format(feature))
    evaluation_complete_filename = '{}/blast/{}_probe_evaluation_complete.txt'.format(sim_dir, feature)
    if not os.path.exists(evaluation_complete_filename):
        evaluate_probes(sim_dir, feature)
    file = open(evaluation_complete_filename, 'w')
    file.write('Feature {} probe evaluation is done.'.format(feature))
    file.close()
    return

if __name__ == '__main__':
    main()
