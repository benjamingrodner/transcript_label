# evaluate_probes
# Purpose:
# get info on blast results and evaluate each probe for specificity
# Ben Grodner and Hao Shi

import pandas as pd
import re
import numpy as np
import subprocess

################################################################################
# Functions
################################################################################


def id_correction(str):
    match = re.search("NR_[0-9]+.[0-9]", str)
    match = match.group(0)
    return(match)


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
                mch_array[i] = np.max(
                    np.append(diffs, [snp_indices[0], len(qseq) - 1 - snp_indices[-1]]))
            else:
                mch_array[i] = len(qseq)
        except TypeError:
            num_line = subprocess.call(['wc', '-l', blast_output_filename])
            print(num_line)
            print(blast_output_filename)
            print(df.ix[i, :])
    return(mch_array)


################################################################################
# Variables
################################################################################
# Snakemake variables
probe_blast_filename = snakemake.input[0]
blast_lineage_strain_filename = snakemake.input[1]
probe_evaluation_filename = snakemake.output[0]

# # Test variables
# probe_blast_filename = "blast_output/synthetic_2019-07-18/probes/1280_consensus.blast.out"
# blast_lineage_strain_filename = "blast_output/synthetic_2019-07-18/blast_lineage_strain.csv"
# probe_evaluation_filename = "probe_evaluation/synthetic_2019-07-18/1280_consensus.csv"

################################################################################
# Script
################################################################################
# Import lineage information
blast_lineage_df = pd.read_csv(blast_lineage_strain_filename, sep="\t")
lineage_columns = ['molecule_id', 'superkingdom', 'phylum',
                   'class', 'order', 'family', 'genus', 'species', 'strain']
blast_lineage_slim = blast_lineage_df[lineage_columns]

# Import probe blast results
# store = pd.HDFStore(probe_evaluation_filename, complib='lzo', complevel=9)
chunk = pd.read_csv(probe_blast_filename, sep="\t", header=None)
chunk.columns = ['probe_id', 'molecule_id', 'pid', 'qcovhsp', 'length', 'mismatch',
                 'gapopen', 'probe_start', 'probe_end', 'molecule_start', 'molecule_end',
                 'evalue', 'bitscore', 'staxids', 'qseq', 'sseq']

# Determine mch for each blast result
chunk['mch'] = mch_test_v4(chunk, probe_blast_filename)

# Merge probe blast results with molecule id lineage information
chunk['molecule_id'] = chunk.molecule_id.apply(id_correction)
chunk = chunk.merge(blast_lineage_slim, on='molecule_id', how='left')

# Write merged results to csv
chunk.to_csv(probe_evaluation_filename, sep='\t')

# chunk_name = re.search('probe_[0-9]*', filename).group(0)
# store.append(chunk_name, chunk, index=False, data_columns=True, min_itemsize={
#              'sseq': 30, 'qseq': 30}, complib='lzo', complevel=9)
# store.close()
