import argparse
import pandas as pd
import os
import glob
import re
import numpy as np
import random
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, generic_dna
from Bio.Blast.Applications import NcbiblastnCommandline
from matplotlib import pyplot as plt
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC



design_dir = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/runs/2020_02_13_gfp_plasmid/simulation/DSGN_008'
design_dir_folder = os.path.split(design_dir)[1]
probe_selection_order_format_filename = '{}/{}_full_length_probe_selection_order_format.xlsx'.format(design_dir, design_dir_folder)
probe_selection_order_format = pd.read_excel(probe_selection_order_format_filename)
# pip install xlrd
# for seq in probe_selection_order_format.Sequence:
sink_df = pd.DataFrame()
for index, full_probe in probe_selection_order_format.iterrows():
    seq = full_probe.Sequence
    probe_name = full_probe.Name
    if 'Probe' in probe_name:
    # seq = probe_selection_order_format.Sequence[1]
        # print(seq)
        probe_split = re.findall(r'(?<=\s)\w+', seq)
        if len(probe_split) == 1:
            probe = probe_split[0]
            # print(probe_name, probe)
        elif len(probe_split) == 2:
            probe = max(probe_split, key = len)
            # print(probe_name, probe)
        elif (len(probe_split) == 3) or (len(probe_split) == 4):
            probe = max(probe_split[0:2], key = len)
            # print(probe_name, probe)
        # elif len(probe_split) == 4:
        #     print(name, probe)
        else:
            print('No spaces between flanking regions/spacers and probe.')
        seq_spaces_removed = re.sub(' ', '', seq)
        # print(seq_spaces_removed)
        sink_probe_section = probe[:-7]
        # # print(sink_probe_section)
        sink_flank_section = re.search(r'\w{5}(?=' + sink_probe_section + ')', seq_spaces_removed).group(0)
        # # print(sink_flank_section.group(0), '\n-')
        sink = sink_flank_section + sink_probe_section
        # print(probe_name)
        probe_name = re.findall(r'(?<=\.)\w+', probe_name)[1]
        # print(probe_name)
        sink_name = 'Sink.' + probe_name
        # print(sink_name)
        sink_df = sink_df.append({'Name': sink_name, 'Sequence':sink}, ignore_index = True)

sink_df['Scale'] = '25nm'
sink_df['Purification'] = 'STD'
probe_selection_order_format = probe_selection_order_format.append(sink_df)
print(probe_selection_order_format)
probe_selection_order_format.to_excel(probe_selection_order_format_filename, index = False)
