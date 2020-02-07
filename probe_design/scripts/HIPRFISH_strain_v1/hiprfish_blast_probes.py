import argparse
import subprocess
import os
from Bio.Blast.Applications import NcbiblastnCommandline
import re
import glob
import pandas as pd

###############################################################################################################
# HiPR-FISH : blast probes
###############################################################################################################

###############################################################################################################
# Workflow functions
###############################################################################################################

def blast_taxon_individual_probe(infile, blast_database, blast_output_filename):
    out_format = '6 qseqid sseqid pident qcovhsp length mismatch gapopen qstart qend sstart send evalue bitscore staxids qseq sseq'
    return_code = subprocess.check_call(['blastn', '-db', blast_database, '-query', infile, '-out', blast_output_filename, '-outfmt', out_format, '-task', 'blastn-short', '-max_hsps', '1', '-max_target_seqs', '100000', '-strand', 'minus', '-evalue', '100', '-num_threads', '1'])
    return(return_code)

def blast_taxon_probes(taxon_probe_directory, blast_database):
    probe_filenames = glob.glob('{}/*.fasta'.format(taxon_probe_directory))
    probe_blast_results = pd.DataFrame(index = probe_filenames, columns = ['blast_return_code'])
    for filename in probe_filenames:
        blast_output_filename = '{}.blast.out'.format(filename)
        return_code = blast_taxon_individual_probe(filename, blast_database, blast_output_filename)
        probe_blast_results.loc[filename, 'blast_return_code'] = return_code
    return(probe_blast_results)

###############################################################################################################
# main function
###############################################################################################################

def main():
    parser = argparse.ArgumentParser('Blast FISH probes designed for a complex microbial community')
    parser.add_argument('input_probe_fasta_dir', type = str, help = 'Input file containing all probes designed by primer3')
    parser.add_argument('oriented_pacbio_filename', type = str, help = 'Input FASTA file containing full length 16S sequences of the complex microbial community')
    args = parser.parse_args()

    design_level_dir = os.path.split(args.input_probe_fasta_dir)[0]
    probes_dir, design_target = os.path.split(design_level_dir)
    blast_complete_filename = '{}_blast_complete.txt'.format(args.input_probe_fasta_dir)
    blast_return_code_filename = '{}_blast_return_code.csv'.format(args.input_probe_fasta_dir)
    probe_blast_results = blast_taxon_probes(args.input_probe_fasta_dir, args.oriented_pacbio_filename)
    probe_blast_results.to_csv(blast_return_code_filename, index = None)
    file = open(blast_complete_filename, 'w')
    file.write('Probe blast for taxon {} is done.'.format(design_target))
    file.close()
    return

if __name__ == '__main__':
    main()
