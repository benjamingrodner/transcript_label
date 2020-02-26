# combine_fasta_files
# Purpose:
# Put all consensus sequences in one file and prepare a primer3 input
# Ben Grodner and Hao Shi

from Bio import SeqIO
import re
import pandas as pd

################################################################################
# Functions
################################################################################


def combine_fasta_files(consensus_fasta_filename, num_rec, taxon,
                        consensus_combined_fasta_filename, primer3_input_filename,
                        probe_design_buffer_start, probe_design_buffer_end):
    for record in SeqIO.parse(consensus_fasta_filename, 'fasta'):
        if num_rec > 1:
            # If multiple consensus seqs,
            # rename consensus sequences with their taxon and a unique id
            record.id = taxon + '_' + record.id + '_consensus'
        else:
            # Rename consensus sequences with their taxon
            record.id = taxon + '_consensus'
        record.description = ''
        # Write sequences to a single output fasta file
        with open(consensus_combined_fasta_filename, 'a') as output_handle:
            SeqIO.write(record, output_handle, 'fasta')
        # Create a row for the primer3 input
        primer3_record = ['SEQUENCE_ID=' + str(record.id),
                          'SEQUENCE_TEMPLATE=' + str(record.seq).upper(),
                          'SEQUENCE_INCLUDED_REGION=' + str(probe_design_buffer_start) + ',' +
                          str(len(record.seq) - probe_design_buffer_end - probe_design_buffer_start),
                          'P3_FILE_FLAG=1', 'PRIMER_EXPLAIN_FLAG=1', '=']
        # Append new row to primer3 input file
        pd.DataFrame(primer3_record).to_csv(primer3_input_filename,
                                            header=None, index=False, mode='a', sep=' ')


################################################################################
# Variables
################################################################################
# Snakemake variables
consensus_fasta_list = snakemake.input[:]
consensus_combined_fasta_filename = snakemake.output[0]
primer3_input_filename = snakemake.output[1]
probe_design_buffer_start = snakemake.params[0]
probe_design_buffer_end = snakemake.params[1]

# # Test variables
# consensus_fasta_list = ["taxon_consensus/synthetic_2019-07-18/1280.consensus.fasta", "taxon_consensus/synthetic_2019-07-18/1747.consensus.fasta",
#                         "taxon_consensus/synthetic_2019-07-18/562.consensus.fasta", "taxon_consensus/synthetic_2019-07-18/1282.consensus.fasta"]
# # consensus_fasta_list = ["taxon_consensus/synthetic_2019-07-18/1280.grouped.fasta"]
# primer3_input_filename = "test.fasta"

################################################################################
# Script
################################################################################
# Iterate over all input files
for consensus_fasta_filename in consensus_fasta_list:
    # Determine the taxon
    taxon = re.search(r'\w+(?=.consensus.fasta)', consensus_fasta_filename)
    taxon = taxon.group(0)

    # Read the fasta file and determine how many sequences there are
    conseq = SeqIO.parse(consensus_fasta_filename, 'fasta')
    num_rec = sum(1 for record in conseq)

    # Write the sequences to a single file, labeled with their taxon and/or consensus id
    # Write the input file for primer3
    combine_fasta_files(consensus_fasta_filename, num_rec, taxon,
                        consensus_combined_fasta_filename, primer3_input_filename,
                        probe_design_buffer_start, probe_design_buffer_end)
