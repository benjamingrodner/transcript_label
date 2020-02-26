import os
import re
import glob
import argparse
import numpy as np
import pandas as pd
from ete3 import NCBITaxa
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, generic_dna
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
import tables

probe_evaluation_complete_filename = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/runs/2020_02_13_gfp_plasmid/sample_008/blast/GFP_probe_evaluation_complete.txt'
probe_summary_info_filename = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/runs/2020_02_13_gfp_plasmid/simulation/DSGN_008/GFP_probe_selection.csv'
probe_evaluation_filename = re.sub('_complete.txt', '.h5', probe_evaluation_complete_filename)
def get_probes(probe_evaluation_filename, input_probe_directory):
    # get feature name
    target_feature = re.sub('_probe_evaluation.h5', '', os.path.basename(probe_evaluation_filename))
    # read probe design/evaluation file to pandas DataFrame
    probes = pd.read_table('{}/{}.int'.format(input_probe_directory, target_feature), skiprows = 3, header = None, delim_whitespace = True)
    # column names
    probes.columns = ['probe_id', 'seq', 'p_start', 'ln', 'N', 'GC', 'Tm', 'self_any_th', 'self_end_th', 'hairpin', 'quality']
    # write target feature to probes DataFrame
    probes.loc[:,'target_feature'] = target_feature
    return(target_feature, probes)

sim_blast_dir, feature_evaluaton_filename = os.path.split(probe_evaluation_filename)
sim_dir = os.path.split(sim_blast_dir)[0]
sim_primer3_dir = '{}/primer3'.format(sim_dir)
# Get probes and feature name from evaluation file
target_feature, probes = get_probes(probe_evaluation_filename, sim_primer3_dir)

probe_selection = pd.read_csv(probe_summary_info_filename)

# print(target_feature, probes)
# build a fasta for the feature
feature_fasta_filename = '{}/{}.fasta'.format(sim_primer3_dir, target_feature)
# initiate a dataframe to summarize all probes
probe_summary = pd.DataFrame()
# iterate through the probe ids
# for probe_idx in probes.probe_id:
    # probe_name = '{}_{}/{}_{}_{}'.format(target_rank, target_taxon, target_rank, target_taxon, probe_idx)
    # generate a new probe name
probe_idx =

probe_name = 'probe_' + str(probe_idx)
    # probe_name = 'probe_' + str(probe_id_idx)
    # read the hdf file with blast and evaluation data
    probe_blast = pd.read_hdf(probe_evaluation_filename, probe_name)
# # Terminal
# probes[probes['probe_id'] == 0]
# probe_blast['qseq']
    # write the melting temp and gc content of each off target match
    probe_blast['melting_temp'] = calculate_tm(probe_blast, Na, dnac1oligo)
    probe_blast['GC_count'] = calculate_gc_count(probe_blast)
    # probe_blast.loc[:,'target_taxon_hit'] = (probe_blast[target_rank].values.astype(str) == str(target_taxon))
    # probe_blast.loc[:,'target_taxon_hit_full_match'] = (probe_blast[target_rank].values.astype(str) == str(target_taxon))*(probe_blast.pid.values >= 99.9)*(probe_blast.qcovhsp.values >= 99.9)
    # probe_summary = probe_summary.append(probe_blast_summarize(probe_blast, max_continuous_homology = max_continuous_homology, taxon_abundance = target_taxon_abundance, target_rank = target_rank), ignore_index = True, sort = False)
    # Write boolean defining if the blast result is on or off target
    probe_blast.loc[:,'target_feature_hit'] = (probe_blast.feature.values.astype(str) == str(target_feature))
    probe_blast.loc[:,'target_feature_hit_full_match'] = (probe_blast.feature.values.astype(str) == str(target_feature))*(probe_blast.pid.values >= 99.9)*(probe_blast.qcovhsp.values >= 99.9)
    # Create a list of the off-target bindings
    probe_blast_off_target = probe_blast[probe_blast.feature.values.astype(str) != str(target_feature)]
    # Filter the probes to a temporary DataFrame using blast on target rate and off-target {bitscore, max melting temperature, and GC content}
    probe_blast_off_target_worst = probe_blast_off_target.loc[(probe_blast_off_target['mch'] > max_continuous_homology) |
                                                            (probe_blast_off_target['bitscore'] > bitscore_thresh) |
                                                            (probe_blast_off_target['melting_temp'] > mt_cutoff) |
                                                            (probe_blast_off_target['GC_count'] > ot_gc_cutoff)]
    # Append the probe blast analysis to the summary of all probes
    probe_summary = probe_summary.append(probe_blast_summarize(probe_blast, max_continuous_homology, probe_blast_off_target_worst), ignore_index = True, sort = True)
