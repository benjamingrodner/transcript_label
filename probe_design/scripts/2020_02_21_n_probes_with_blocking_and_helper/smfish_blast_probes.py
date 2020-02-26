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

def blast_feature_individual_probe(infile, blast_database, blast_output_filename, molecule):
    # read in probe information
    # print('    Blasting ' + os.path.basename(infile))
    out_format = '6 qseqid sseqid pident qcovhsp length mismatch gapopen qstart qend sstart send evalue bitscore staxids qseq sseq'
    # try:
    #     os.path.exists(blast_output_filename)
    #     return_code = 2
    # except:
        # blastn_cline = NcbiblastnCommandline(cmd = 'blastn', query = infile, db = blast_database, outfmt = '"' + out_format + '"', out = blast_output_filename, task = 'blastn-short', max_hsps = 1, max_target_seqs = 100000, strand = 'minus', evalue = 1000, num_threads = 4)
        # result = blastn_cline()
    if molecule == 'RNA':
        return_code = subprocess.check_call(['blastn', '-db', blast_database, '-query', infile, '-out', blast_output_filename, '-outfmt', out_format, '-task', 'blastn-short', '-max_hsps', '1', '-max_target_seqs', '100000', '-strand', 'minus', '-evalue', '100', '-num_threads', '1'])
    elif molecule == 'DNA':
        return_code = subprocess.check_call(['blastn', '-db', blast_database, '-query', infile, '-out', blast_output_filename, '-outfmt', out_format, '-task', 'blastn-short', '-max_hsps', '1', '-max_target_seqs', '100000', '-strand', 'plus', '-evalue', '100', '-num_threads', '1'])
    # return_code = subprocess.check_call(['blastn', '-db', blast_database, '-query', infile, '-out', blast_output_filename, '-outfmt', out_format, '-task', 'blastn-short', '-max_hsps', '1', '-max_target_seqs', '100000', '-evalue', '100', '-num_threads', '1'])
    return(return_code)

def blast_feature_probes(feature_probe_directory, blast_database, molecule):
    probe_filenames = glob.glob(feature_probe_directory + '/*.fasta')
    probe_blast_results = pd.DataFrame(index = probe_filenames, columns = ['blast_return_code'])
    for filename in probe_filenames:
        blast_output_filename = filename + '.blast.out'
        return_code = blast_feature_individual_probe(filename, blast_database, blast_output_filename, molecule)
        probe_blast_results.loc[filename, 'blast_return_code'] = return_code
    return(probe_blast_results)


###############################################################################################################
# main function
###############################################################################################################

def main():
    parser = argparse.ArgumentParser('Blast FISH probes designed for a complex microbial community')
    # input 16S full length sequence
    # input probes file

    parser.add_argument('custom_blast_db', type = str, help = 'Input FASTA file containing full length 16S sequences of the complex microbial community')

    parser.add_argument('input_probes_filename', type = str, help = 'Input file containing all probes designed by primer3')
    parser.add_argument('-mol', '--molecule', dest = 'molecule', type = str, default = 'RNA', help = 'Is the target DNA or RNA (i.e. are we hybridizing to non-coding or coding)?')

    args = parser.parse_args()

    sim_primer3_dir, feature_probes_filename = os.path.split(args.input_probes_filename)
    feature = re.sub('.int', '', feature_probes_filename)
    feature_probe_directory = '{}/{}'.format(sim_primer3_dir, feature)
    blast_complete_filename = '{}/{}_probe_blast_complete.txt'.format(sim_primer3_dir, feature)
    if os.path.exists(blast_complete_filename):
        print('I am skipping blasting')
    else:
        print('I am here doing blasting')
        probe_blast_results = blast_feature_probes(feature_probe_directory, args.custom_blast_db, args.molecule)
        probe_blast_results.to_csv('{}/{}_blast_return_code.csv'.format(sim_primer3_dir, feature))
        file = open(blast_complete_filename, 'w')
        file.write('Feature {} probe blast is done.'.format(feature))
        file.close()
    return

if __name__ == '__main__':
    main()
