# split_taxon_probe_file
# Purpose:
# Make individual fasta files for each probe
# Ben Grodner and Hao Shi

from Bio import SeqIO
import re
import pandas as pd
import os

################################################################################
# Functions
################################################################################

################################################################################
# Variables
################################################################################
# # Snakemake variables
# probe_design_input_filename = snakemake.input[0]
#
# probe_design_output_filename = snakemake.output[0]
# probe_fasta_filename = snakemake.output[1]

# # Test variables
probe_design_input_filename = "1280_consensus.int"

probe_design_output_filename = "test/1280_consensus/probe_design.int"
probe_fasta_filename = "test/1280_consensus/sequences.int"


################################################################################
# Script
################################################################################
# Iterate over all input files
probes = pd.read_table(probe_design_input_filename, skiprows=3, header=None, delim_whitespace=True)
probes.columns = ['probe_num', 'seq', 'start', 'length', 'N', 'GC',
                  'Tm', 'self_any_th', 'self_end_th', 'hair-pin', 'quality']

# Create a list of
probes_list = [SeqRecord(Seq(probes['seq'][i]).reverse_complement(), id=str(
    probes['probe_num'][i]), description='') for i in range(0, probes.shape[0])]
probes_fasta_filenames = [probe_dir + '/' + taxon + '/' + taxon_consensus_filename +
                          '.probe_' + str(probe_num) + '.fasta' for probe_num in range(probes.shape[0])]
for i in range(probes.shape[0]):
    SeqIO.write(probes_list[i], probes_fasta_filenames[i], 'fasta')
