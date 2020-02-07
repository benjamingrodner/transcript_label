import argparse
import subprocess
import os
from Bio.Blast.Applications import NcbiblastnCommandline
import re
import glob
import pandas as pd

###############################################################################################################
# HiPR-FISH-strain : design probes
###############################################################################################################

def design_probes(primer3_input_filename, primer3_settings_filename, primer3_output_filename):
    return_code = subprocess.check_call(['/programs/primer3-2.3.5/src/primer3_core', '-p3_settings_file', primer3_settings_filename, '-output', primer3_output_filename, '-format_output', primer3_input_filename])
    probe_int_filename = re.sub('_primer3_input.txt', '.int', primer3_input_filename)
    probe_csv_filename = re.sub('_primer3_input.txt', '_probes.csv', primer3_input_filename)
    probes = pd.read_table(probe_int_filename, skiprows = 3, header = None, delim_whitespace = True)
    probes.columns = ['probe_id', 'seq', 'start', 'length', 'N', 'GC', 'Tm', 'self_any_th', 'self_end_th', 'hair-pin', 'quality']
    probes['source'] = probe_csv_filename
    probes.to_csv(probe_csv_filename, index = None)
    return

###############################################################################################################
# main function
###############################################################################################################

def main():
    parser = argparse.ArgumentParser('Blast FISH probes designed for a complex microbial community')

    parser.add_argument('primer3_input_filename', type = str, help = 'Input FASTA file containing full length 16S sequences of the complex microbial community')

    args = parser.parse_args()
    primer3_dir = os.path.split(args.primer3_input_filename)[0]
    primer3_settings_filename = '{}/primer3_settings.txt'.format(primer3_dir)
    output_filename = re.sub('_primer3_input.txt', '_primer3_output.txt', args.primer3_input_filename)
    os.chdir(primer3_dir)
    design_probes(args.primer3_input_filename, primer3_settings_filename, output_filename)
    return

if __name__ == '__main__':
    main()
