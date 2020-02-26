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

###############################################################################################################
# HiPR-FISH Probe Design Pipeline
###############################################################################################################

###############################################################################################################
# Workflow functions
###############################################################################################################

def cm_to_inches(x):
    return(x/2.54)

def add_spacer(feature_best_probes, sim_primer3_dir, output_filename):
    probes = pd.read_csv(feature_best_probes)
    probes['seqrcsa'] = ''
    for i in range(0, probes.shape[0]):
        probe_seq = Seq(probes['seq'][i], generic_dna)
        rrna_file = '{}/{}.fasta'.format(sim_primer3_dir, str(probes.loc[i, 'target_feature']))
        rrna_seq = SeqIO.read(rrna_file, 'fasta').seq
        sstart = probes['p_start'][i]
        probe_length = probes['ln'][i]
        send = sstart + probe_length - 1
        right_spacer = rrna_seq[sstart - 2] + rrna_seq[sstart - 3] + rrna_seq[sstart - 4]
        left_spacer = rrna_seq[send + 2] + rrna_seq[send + 1] + rrna_seq[send]
        probe_seq_rcsa = str(left_spacer).upper() + str(probe_seq.reverse_complement()) + str(right_spacer).upper()
        probes.ix[i, 'seqrcsa'] = probe_seq_rcsa
    probes.loc[:,'quadg'] = (probes['seqrcsa'].str.upper().str.count('GGGG')) + (probes['seqrcsa'].str.upper().str.count('GGGGG'))
    probes = probes[probes['quadg'] == 0]
    probes.to_csv(output_filename, index = False)
    return

def convert_numeric_barcode(num, nbit):
    code = re.sub('0b', '', format(num, '#0' + str(nbit+2) + 'b'))
    return(code)

def count_number_of_bits(binary_barcode):
    bin_list = list(binary_barcode)
    bin_int_list = [int(index) for index in bin_list]
    return(np.sum(bin_int_list))

# def generate_full_probes(design_dir, bot, plf = 'T', primer = 'T', primerset = 'B', barcode_selection = 'MostSimple'):
#     design_id = os.path.basename(design_dir)
#     RP = ['GATGATGTAGTAGTAAGGGT',
#           'AGGTTAGGTTGAGAATAGGA',
#           'AGGGTGTGTTTGTAAAGGGT',
#           'TTGGAGGTGTAGGGAGTAAA',
#           'AGAGTGAGTAGTAGTGGAGT',
#           'ATAGGAAATGGTGGTAGTGT',
#           'TGTGGAGGGATTGAAGGATA']
#     nbit = len(RP)
#     if plf == 'T':
#         # probe_length_filler = 'GGAATCGATGGTGCACTGCT'
#         probe_length_filler = 'GTCTATTTTCTTATCCGACG'
#     else:
#         probe_length_filler = ''
#     if primer == 'T':
#         if primerset == 'A':
#             forward_primer = 'CGCGGGCTATATGCGAACCG'
#             reverse_primer = 'GCGTTGTATGCCCTCCACGC'
#             # TAATACGACTCACTATAGGGCGTGGAGGGCATACAACGC
#         elif primerset == 'B':
#             forward_primer = 'CGATGCGCCAATTCCGGTTC'
#             reverse_primer = 'CAACCCGCGAGCGATGATCA'
#             # TAATACGACTCACTATAGGGTGATCATCGCTCGCGGGTTG
#         elif primerset == 'C':
#             forward_primer = 'GTTGGTCGGCACTTGGGTGC'
#             reverse_primer = 'CCACCGGATGAACCGGCTTT'
#             # TAATACGACTCACTATAGGGAAAGCCGGTTCATCCGGTGG
#     else:
#         forward_primer = ''
#         reverse_primer = ''
#     feature_best_probes_sa_filenames = glob.glob('{}/*_probe_selection_sa.csv'.format(design_dir))
#     feature_best_probes_sa_filenames.sort()
#     feature_best_probes_list = [pd.read_csv(filename) for filename in feature_best_probes_sa_filenames]
#     feature_best_probes_filtered_list = [x for x in feature_best_probes_list if x.blast_on_target_rate.max() > bot]
#     oligo_list = []
#     blocking_probe_list = []
#     barcodes = pd.DataFrame(np.arange(1,2**nbit))
#     barcodes.columns = ['NumericBarcode']
#     barcodes['BinaryBarcode'] = barcodes.NumericBarcode.apply(convert_numeric_barcode, args = (7,))
#     barcodes['NumberBits'] = barcodes.BinaryBarcode.apply(count_number_of_bits)
#     if barcode_selection == 'MostSimple':
#         barcodes_sorted = barcodes.sort_values(by = ['NumberBits', 'NumericBarcode'])
#     elif barcode_selection == 'MostComplex':
#         barcodes_sorted = barcodes.sort_values(by = ['NumberBits', 'NumericBarcode'], ascending = [False, False])
#     elif barcode_selection == 'Random':
#         barcodes_sorted = barcodes.reindex(np.random.permutation(barcodes.index))
#     elif barcode_selection == 'Sequential':
#         barcodes_sorted = barcodes.copy()
#     for i in range(min(127,len(feature_best_probes_filtered_list))):
#         probes = feature_best_probes_filtered_list[i]
#         probes = probes[(probes['blast_on_target_rate'] > bot)]
#         assigned = barcodes_sorted.NumericBarcode.values[i]
#         if barcodes_sorted.NumberBits.values[i] > 2:
#             barcode_repeat = np.round(15/barcodes_sorted.NumberBits.values[i]).astype(int)
#         else:
#             barcode_repeat = 15
#         if probes.shape[0] > 0:
#             for k in range(barcode_repeat):
#                 probes = probes.sort_values(by = 'quality', ascending = True).reset_index().drop(columns = ['index'])
#                 for i in range(0, probes.shape[0]):
#                     tarseq = probes['seqrcsa'][i]
#                     if 'N' not in list(tarseq):
#                         code = np.asarray(list(re.sub('0b', '', format(assigned, '#0' + str(7+2) + 'b'))), dtype = np.int64)
#                         indices = np.where(code == 1)[0]
#                         if len(indices) > 2:
#                             indices = np.append(indices, indices[0])
#                             subcode = np.zeros((len(indices)-1, nbit), dtype = np.int64)
#                             for j in range(0, len(indices) - 1):
#                                 subcode[j, indices[j]] = 1
#                                 subcode[j, indices[j+1]] = 1
#                             oligo = [[str(probes['target_feature'][i]), probes['p_start'][i], probes['ln'][i], tarseq, probes['Tm'][i], probes['GC'][i], re.sub('\[|\]| ','',str(code)), assigned, probes['probe_id'][i], str(probes['target_feature'][i]) + '_' + str(assigned) + '_' + re.sub('\[|\]| ','',str(subcode[k])) + '_' + str(probes['probe_id'][i]), forward_primer + RP[np.where(subcode[k] == 1)[0][0]] + tarseq + RP[np.where(subcode[k] == 1)[0][1]] + reverse_primer] for k in range(0, len(subcode))]
#                         elif len(indices) == 2:
#                             oligo = [[str(probes['target_feature'][i]), probes['p_start'][i], probes['ln'][i], tarseq, probes['Tm'][i], probes['GC'][i], re.sub('\[|\]| ','',str(code)), assigned, probes['probe_id'][i], str(probes['target_feature'][i]) + '_' + str(assigned) + '_' + re.sub('\[|\]| ','', str(code)) + '_' + str(probes['probe_id'][i]), forward_primer + RP[np.where(code == 1)[0][0]] + tarseq + RP[np.where(code == 1)[0][1]] + reverse_primer]]
#                         else:
#                             oligo = [[str(probes['target_feature'][i]), probes['p_start'][i], probes['ln'][i], tarseq, probes['Tm'][i], probes['GC'][i], re.sub('\[|\]| ','',str(code)), assigned, probes['probe_id'][i], str(probes['target_feature'][i]) + '_' + str(assigned) + '_' + re.sub('\[|\]| ','',str(code)) + '_' + str(probes['probe_id'][i]), forward_primer + RP[np.where(code == 1)[0][0]] + tarseq + probe_length_filler + reverse_primer]]
#                         oligo_list = oligo_list + oligo
#     oligo_df = pd.DataFrame(oligo_list)
#     oligo_df.columns = ['target_feature', 'p_start', 'ln', 'rna_seq', 'Tm', 'GC', 'code', 'numeric_code', 'probe_numeric_id', 'probe_id', 'probe_full_seq']
#     return(oligo_df)


def generate_full_probes(design_dir, target_table_filename, blast_database, mch, bitscore_thresh, ot_gc_cutoff, mt_cutoff, dnaconc, sod, molecule, probe_selection_level):
# # Terminal
# data_dir = '/workdir/bmg224/hiprfish/transcript_label/probe_design/runs/2020_02_05_hao_zhou_mouse_gut/'
# design_dir = '/workdir/bmg224/hiprfish/transcript_label/probe_design/runs/2020_02_05_hao_zhou_mouse_gut/simulation/DSGN_002'
# # input_dir = '/workdir/bmg224/hiprfish/transcript_label/probe_design/runs/2020_02_05_hao_zhou_mouse_gut/sample_002/inputs/'
# target_table_filename = glob.glob('{}/{}/inputs/*.csv'.format(data_dir, sample))
    # Collect the filenames for the selected probes
    # feature_best_probes_filenames = glob.glob('{}/*_probe_selection.csv'.format(design_dir))
    #     # Read all the tables to a list
    # feature_best_probes_list = [pd.read_csv(filename) for filename in feature_best_probes_filenames]
    #     # Concatenate all the selected probes
    # probe_selection = pd.concat(feature_best_probes_list)
    # print(probe_selection.shape[0])
    # print(probe_selection.target_feature)
    # Read in the table with the naming and flanking region info
    target_table = pd.read_csv(target_table_filename)
    # print(target_table.target_feature)
    # target_table = pd.read_csv('/workdir/bmg224/hiprfish/mobile_elements/probe_design/runs/2020_02_13_gfp_plasmid/sample_004/inputs/target_table_004.csv')
    # Merge the info table with the probe selection
    # probe_selection = probe_selection.merge(target_table, on = 'target_feature', how = 'left')
    # Start a new dataframe to recieve full length seqs and names
    full_length_names_and_seqs = pd.DataFrame()
    # Define flanking region sequences
    flanking_regions = pd.DataFrame({'readout':['F3','F7'], 'readout_sequence':['GGGAGAATGAGGTGTAATGT','TTGGAGGTGTAGGGAGTAAA']})
    nucleotide = ['a', 't', 'g', 'c']
    spacers = [' ', 'a', 't', 'g', 'c']
    spacer_df = pd.DataFrame(columns = ['spacer_name', 'spacer_1', 'spacer_2'])
    for i in nucleotide:
        for j in nucleotide:
            spacers.append(i + j)

    for i in spacers:
        for j in spacers:
            spacer_name = i + '_' + j
            spacer_name = re.sub(' ', '', spacer_name)
            spacer_1 = i
            spacer_2 = j
            temp = pd.DataFrame({'spacer_name': [spacer_name], 'spacer_1': [spacer_1], 'spacer_2': [spacer_2]})
            spacer_df = spacer_df.append(temp)


    full_length_names_and_seqs = pd.DataFrame()
    probe_selection = pd.DataFrame()
    for feature in target_table['target_feature'].unique():
        feature_best_probes_filename = '{}/{}_probe_selection.csv'.format(design_dir, feature)
        feature_best_probes = pd.read_csv(feature_best_probes_filename)
        feature_best_probes = feature_best_probes.merge(target_table, on = 'target_feature', how = 'left')
        num_readouts = target_table[target_table['target_feature'] == feature]['num_readouts'].values[0]
        for i in range(num_readouts):
            for probe_num in range(int(np.floor(i * feature_best_probes.shape[0] / num_readouts)), int(np.floor((i + 1) * feature_best_probes.shape[0] / num_readouts))):
            # for probe_num in range(int(np.floor(i * probe_selection.shape[0] / num_readouts)), int(np.floor((i + 1) * probe_selection.shape[0] / num_readouts))):
                probe_selection_temp = feature_best_probes.iloc[probe_num, :]
                readout = 'readout_' + str(i+1)
                name = 'Probe.' + probe_selection_temp[readout] + '.' + probe_selection_temp['target_name'] + '_' + probe_selection_temp['p_start'].astype(int).astype(str) + '.' + probe_selection_temp[readout]
                probe_id = probe_selection_temp['probe_id']
                print('-\nfeature\n',feature,'\nname\n',name,'\nprobe_id\n',probe_id,'\nsequence\n',probe_selection_temp['seq'])
                readout_sequence = flanking_regions[flanking_regions['readout'] == probe_selection_temp[readout]]['readout_sequence'].values[0]
                # temp = pd.DataFrame({'Name':[name], 'probe_id':[probe_id]})
                spacer_fasta_filename = '{}/{}_spacers/{}.fasta'.format(design_dir, feature, probe_id)
                blast_output_filename = re.sub('.fasta', '.blast.out', spacer_fasta_filename)
                seq_df = pd.DataFrame()
                for index, spacer in spacer_df.iterrows():
                    spacer_sum = spacer['spacer_1'] + spacer['spacer_2']
                    spacer_length = len(spacer_sum) - spacer_sum.count(' ')
                    diff = len(probe_selection_temp['seq']) - (20 - spacer_length)
                    if diff < 0:
                        diff = 0
                    # print(diff, spacer_length, spacer['spacer_1'], spacer['spacer_2'])
                    sequence = readout_sequence[diff:] + ' ' + spacer['spacer_1'] + ' ' + probe_selection_temp['seq'] + ' ' + spacer['spacer_2'] + ' ' + readout_sequence
                    # if probe_selection_temp.probe_id == 681:
                    #     print('spacer_name',spacer['spacer_name'],'seq',len(probe_selection_temp['seq']),'spacer',spacer_length,'diff',diff, 'Sequence', sequence)
                    temp_seq = pd.DataFrame({'spacer_name':[spacer['spacer_name']], 'Sequence':[sequence]})
                    seq_df = seq_df.append(temp_seq)
                spacer_seqs = spacer_df.merge(seq_df, on = 'spacer_name', how = 'left')
                if not os.path.exists(spacer_fasta_filename) or not os.path.exists(blast_output_filename):
                    full_length_fasta = []
                    for index, spacer in spacer_df.iterrows():
                        spacer_sum = spacer['spacer_1'] + spacer['spacer_2']
                        spacer_length = len(spacer_sum) - spacer_sum.count(' ')
                        diff = len(probe_selection_temp['seq']) - (20 - spacer_length)
                        if diff < 0:
                            diff = 0
                        # print(diff, spacer_length, spacer['spacer_1'], spacer['spacer_2'])
                        sequence = readout_sequence[diff:] + ' ' + spacer['spacer_1'] + ' ' + probe_selection_temp['seq'] + ' ' + spacer['spacer_2'] + ' ' + readout_sequence
                        record = SeqRecord(Seq(sequence), id = spacer['spacer_name'])
                        full_length_fasta.append(record)
                    if not os.path.exists('{}/{}_spacers'.format(design_dir, feature)):
                        os.mkdir('{}/{}_spacers'.format(design_dir, feature))
                    # name = re.sub('.fna','',name)
                    # spacer_fasta_filename = '{}/{}_spacers/{}.fasta'.format(design_dir, feature, probe_id)
                    # combined_filename
                    with open(spacer_fasta_filename, 'w') as output_handle:
                        SeqIO.write(full_length_fasta, output_handle, 'fasta')
                    # print(spacer_fasta_filename)
                    # blast_output_filename = re.sub('.fasta', '.blast.out', spacer_fasta_filename)
                    blast_probes(spacer_fasta_filename, blast_output_filename, blast_database, molecule)
                # print(blast_output_filename)
                probes_blast = pd.read_csv(blast_output_filename, header = None, delim_whitespace = True)
                probes_blast.columns = ['probe_id', 'feature', 'pid', 'qcovhsp', 'length', 'mismatch', 'gapopen', 'probe_start', 'probe_end', 'molecule_start', 'molecule_end', 'evalue', 'bitscore', 'staxids', 'qseq', 'sseq']
                # print(probes_blast.qseq.values[0])
                probes_blast_significant = get_significant_blast_homologies(probes_blast, mch = mch, mt_cutoff = mt_cutoff, ot_gc_cutoff = ot_gc_cutoff, bitscore_thresh = bitscore_thresh, dnaconc = dnaconc, sod = sod)
                # print('probes_blast_significant spacer design', probes_blast_significant.shape[])
                # probe_id_to_feature = probe_selection[['probe_id','target_feature']]
                # probes_blast_significant = probes_blast_significant.merge(probe_id_to_feature, on = 'probe_id', how = 'left')
                # probes_blast_significant_ot = pd.DataFrame()
                # for index, probe in feature_best_probes.iterrows():
                #     feature = probe['target_feature']
                #     target_table = pd.read_csv(target_table_filename)
                #     # probe_selection_level = target_table[target_table['target_feature'] == feature]['selection_level'].iloc[0]
                    # print(probe_selection_level, probe_selection_level == 'molecule')
                if probe_selection_level == 'genome_bin':
                    total_transcript_table = pd.read_csv(total_transcript_table_filename)
                    feature_transcript_info = total_transcript_table[total_transcript_table['Gene'] == str(feature)]
                    phylogenetic_marker_genes = ['dnaG', 'frr', 'infC', 'nusA', 'pgk', 'pyrG', 'rplA', 'rplB', 'rplC', 'rplD', 'rplE', 'rplF', 'rplK', 'rplL', 'rplM', 'rplN',
                    'rplP', 'rplS', 'rplT', 'rpmA', 'rpoB', 'rpsB', 'rpsC', 'rpsE', 'rpsI', 'rpsJ', 'rpsK', 'rpsM', 'rpsS', 'smpB', 'tsf']
                    feature_genome_bin = total_transcript_table[total_transcript_table['Genome_bin'] == feature_transcript_info['Genome_bin']]
                    feature_genome_bin_marker = feature_genome_bin[feature_genome_bin['Preferred_name'].isin(phylogenetic_marker_genes)]
                    genome_bin_feature_list = feature_genome_bin_marker['Gene'].values
                    probes_blast_significant_ot_temp =  probe_blast[~probe_blast.feature.isin(genome_bin_feature_list)]
                    probes_blast_significant_ot = probes_blast_significant_ot.append(probes_blast_significant_ot_temp)
                elif probe_selection_level == 'molecule':
                    # # print('True')
                    probes_blast_significant_ot = probes_blast_significant[probes_blast_significant['feature'] != feature]
                    # print('on target\n',probes_blast_significant[probes_blast_significant['feature'] == feature].shape[0])
                    # probes_blast_significant_ot_temp = probes_blast_significant[probes_blast_significant['feature'] != feature]
                    # # print(probes_blast_significant_ot_temp.columns)
                    # probes_blast_significant_ot = probes_blast_significant_ot.append(probes_blast_significant_ot_temp)
                # probes_blast_significant_ot = probes_blast_significant[probes_blast_significant['target_feature'] != probes_blast_significant['feature']]
                blocking_probes_filename = '{}/{}_blocking_probes.csv'.format(design_dir, feature)
                try:
                    blocking_probes = pd.read_csv(blocking_probes_filename)
                except:
                    blocking_probes = pd.DataFrame(columns = ['sseq'])
                # blocking_probes = blocking_probes.merge(probe_id_to_feature, on = 'probe_id')
                # print('probes_blast_significant spacer design', probes_blast_significant.shape,'probes_blast_significant_ot spacer design', probes_blast_significant_ot.shape)
                probes_blast_significant_ot_unblocked = probes_blast_significant_ot[~probes_blast_significant_ot['sseq'].isin(blocking_probes['sseq'])]
                count_df = pd.DataFrame()
                for index, spacer in spacer_seqs.iterrows():
                    count = probes_blast_significant_ot_unblocked[probes_blast_significant_ot_unblocked['probe_id'] == spacer['spacer_name']].shape[0]
                    spacer_sum = spacer['spacer_1'] + spacer['spacer_2']
                    spacer_length = len(spacer_sum) - spacer_sum.count(' ')
                    spacer_name = spacer['spacer_name']
                    temp_count = pd.DataFrame({'spacer_name': [spacer_name], 'count': [count], 'length': [spacer_length]})
                    count_df = pd.concat([count_df, temp_count], ignore_index = True)
                    # count_df = count_df.append(temp_count)
                spacer_seqs = spacer_seqs.merge(count_df, on = 'spacer_name', how = 'left')
                spacer_seqs.sort_values(by = ['count','length'], ascending = [True, True], inplace = True)
                # spacer_choice_sequence = spacer_df[spacer_df['count'] == spacer_df['count'].min()]
                # spacer_choice_sequence = spacer_choice_sequence[spacer_choice_sequence['length'] == spacer_choice_sequence['length'].min()]

                # spacer_choice_sequence.iloc[0,:]
                spacer_choice_sequence = spacer_seqs.iloc[0,:]['Sequence']
                print('probe_blast_off_target\n',probes_blast_significant_ot_unblocked[probes_blast_significant_ot_unblocked['probe_id'] == spacer_seqs.iloc[0,:]['spacer_name']]['sseq'])
                # print('probe_id',probe_id,probes_blast_significant_ot_unblocked[probes_blast_significant_ot_unblocked['probe_id'] == spacer_seqs.iloc[0,:]['spacer_name']]['qseq'])
                print('spacer_choice_sequence\n',spacer_choice_sequence)
                temp_spacer_choice_sequence = pd.DataFrame({'probe_id':[probe_id], 'Sequence':[spacer_choice_sequence], 'Name':[name]})
                full_length_names_and_seqs = full_length_names_and_seqs.append(temp_spacer_choice_sequence)


                # diff = len(probe_selection_temp['seq']) - 20
                # full_length_sequence = readout_sequence[diff:] + ' ' + probe_selection_temp['seq'] + ' ' +  readout_sequence
                # full_length_names_and_seqs = pd.concat([full_length_names_and_seqs, temp])
        feature_best_probes = feature_best_probes.merge(full_length_names_and_seqs, on = 'probe_id', how = 'left')
        probe_selection = probe_selection.append(feature_best_probes)
    print('probe_selection from full probes\n',probe_selection[['probe_id','Name', 'Sequence']])
    # print(probe_selection.Name)
    # probe_selection = probe_selection.merge(full_length_names_and_seqs, on = 'probe_id', how = 'left')
    # probe_selection.columns
    # # Generate a name column that includes flanking region, target, and probe location on the target molecule info
    # probe_selection['Name'] = 'Probe.' + probe_selection['readout_1'] + '.' + probe_selection['target_name'] + '_' + probe_selection['p_start'].astype(int).astype(str) + '.' + probe_selection['readout_1']
    # # Create a DataFrame for the flanking regions
    # flanking_regions = pd.DataFrame({'readout':['F3','F7'], 'readout_sequence':['GGGAGAATGAGGTGTAATGT','TTGGAGGTGTAGGGAGTAAA']})
    # # Merge the flanking regions with the probe_selection
    # probe_selection = probe_selection.merge(flanking_regions, left_on = 'readout_1', right_on = 'readout', how = 'left')
    # # Create a new column with the merged sequences and flanking regions
    # probe_selection['Sequence'] = probe_selection['readout_sequence'] + ' ' + probe_selection['seq'] + ' ' +  probe_selection['readout_sequence']
    # Add IDT columns
    probe_selection.loc[:, 'Scale'] = '25nm'
    probe_selection.loc[:, 'Purification'] = 'STD'
    # Check if the length is greater than 60
    for seq in probe_selection.Sequence.values:
        if len(seq) - seq.count(' ') > 60:
            probe_selection.loc[probe_selection.Sequence.values == seq, 'Scale'] = '100nm'
# probe_selection_full_length_filename = '{}/probe_selection_full_length.csv'.format(design_dir)
# probe_selection.to_csv(probe_selection_full_length_filename)
    # helper_probes_filename = glob.glob('{}/*_helper_probes.csv'.format(design_dir))
    # # Read all the tables to a list
    # helper_probes_list = [pd.read_csv(filename) for filename in helper_probes_filename]
    # # Concatenate all the selected probes
    # helper_probes = pd.concat(helper_probes_list)
    # # Read in the table with the naming and flanking region info
    # # Add IDT columns
    # helper_probes.loc[:, 'Scale'] = '25nm'
    # helper_probes.loc[:, 'Purification'] = 'STD'
    # target_table_filename = '{}/'.format(os.path.split(os.path.split(design_dir)[0])[0])
    # target_table_filename = '{}/*.csv'.format(input_dir)
    # probe_selection[''] = re.search('(?<=NODE)',probe_selection['target_feature'])
    # probe_selection_order_format = probe_selection[['Name','Sequence','Scale','Purification']]
    # probe_selection_order_format.to_excel('{}/full_length_probe_selection_order_format.xlsx'.format(design_dir))
    return(probe_selection)

def generate_blocking_probes(design_dir):
        # Initiate DataFrames to output the off target binding of selected probes
        probe_off_target_summary = pd.DataFrame()
        probe_blast_off_target_all = pd.DataFrame()
        blocking_probes_filenames = glob.glob('{}/*_blocking_probes.csv'.format(design_dir))
        try:
            blocking_probes_list = [pd.read_csv(file) for file in blocking_probes_filenames]
            blocking_probes = pd.concat(blocking_probes_list)
        except:
            blocking_probes = pd.DataFrame(columns = [])
        # Iterate through all selected probes
        for probe_idx in best_probes.probe_id:
            # ensure id is an integer
            probe_idx = int(probe_idx)
            # Get probe names
            probe_info = best_probes.loc[best_probes.probe_id == probe_idx]
            print('-\nProbe {}\n-\n{}'.format(probe_idx, probe_info[['seq','off_target_worst_count']]))
            # Get the hdf probe name
            probe_name = 'probe_' + str(probe_idx)
            # Read in the probe blast results
            probe_blast = pd.read_hdf(probe_evaluation_filename, probe_name)
            # # write the melting temp and gc content of each off target match
            probe_blast['melting_temp'] = calculate_tm(probe_blast, Na, dnac1oligo)
            probe_blast['GC_count'] = calculate_gc_count(probe_blast)
            probe_blast_off_target = probe_blast[probe_blast.feature.values.astype(str) != str(target_feature)]
            probe_blast_off_target.sort_values(['qcovhsp', 'mch', 'bitscore', 'melting_temp', 'GC_count'], ascending = [False, False, False, False, False], inplace = True)
            probe_blast_off_target_all = probe_blast_off_target_all.append(probe_blast_off_target)
            # Check if the probe needs blocking probes
            if probe_info.blocking.values[0]:
                # # Get the hdf probe name
                # probe_name = 'probe_' + str(probe_idx)
                # # Read in the probe blast results
                # probe_blast = pd.read_hdf(probe_evaluation_filename, probe_name)
                # # # write the melting temp and gc content of each off target match
                # probe_blast['melting_temp'] = calculate_tm(probe_blast, Na, dnac1oligo)
                # probe_blast['GC_count'] = calculate_gc_count(probe_blast)

                # # Create a list of the off-target bindings
                # probe_blast_off_target = probe_blast[probe_blast.feature.values.astype(str) != str(target_feature)]
                # # Sort the off target probes by off-target binding parameters
                # probe_blast_off_target.sort_values(['qcovhsp', 'mch', 'bitscore', 'melting_temp', 'GC_count'], ascending = [False, False, False, False, False], inplace = True)
                # print(probe_blast_off_target[['qcovhsp','mch','evalue','bitscore','melting_temp','GC_count']].iloc[0:20,:])
                # Create a list of the off-target bindings
                # probe_blast_off_target = probe_blast[probe_blast.feature.values.astype(str) != str(target_feature)]
                # probe_blast_off_target.sort_values(['qcovhsp', 'mch', 'bitscore', 'melting_temp', 'GC_count'], ascending = [False, False, False, False, False], inplace = True)
                # Filter the probes to a temporary DataFrame using blast on target rate and off-target {bitscore, max melting temperature, and GC content}
                probe_blast_off_target_worst = probe_blast_off_target.loc[(probe_blast_off_target['mch'] > max_continuous_homology) |
                                                                        (probe_blast_off_target['bitscore'] > bitscore_thresh) |
                                                                        (probe_blast_off_target['melting_temp'] > mt_cutoff) |
                                                                        (probe_blast_off_target['GC_count'] > ot_gc_cutoff)]
                print('There are {} blocking probes for probe {}'.format(probe_blast_off_target_worst.shape[0],probe_idx))
                blocking_probes = blocking_probes.append(probe_blast_off_target_worst)

def generate_blocking_probes(design_dir, bplc, plf = 'T', primer = 'T', primerset = 'B', barcode_selection = 'MostSimple'):
    design_dir_folder = os.path.split(design_dir)[1]
    blocking_probes_filenames = glob.glob('{}/*_off_target_summary_info.csv'.format(design_dir))
    probes_filename = '{}/{}_primerset_{}_barcode_selection_{}_full_length_probes.csv'.format(design_dir, design_dir_folder, primerset, barcode_selection)
    probes = pd.read_csv(probes_filename)
    blocking_probes_list = []
    blocking_probes_list = [pd.read_csv(filename) for filename in blocking_probes_filenames]
    blocking_probes = pd.concat(blocking_probes_list, sort = True)
    blocking_probes = blocking_probes.loc[blocking_probes.length.values >= bplc]
    blocking_probes.to_csv('{}/{}_primerset_{}_barcode_selection_{}_full_length_blocking_probes.csv'.format(design_dir, design_dir_folder, primerset, barcode_selection))
    blocking_probes_order_format = blocking_probes[['sseq', 'length']]
    blocking_probes_order_format = blocking_probes_order_format.drop_duplicates()
    blocking_probes_order_format = blocking_probes_order_format.sort_values(by = 'length', ascending = False)
    blocking_probes_order_format.insert(0, 'blocking_probe_name', ['blocking_probe_{}'.format(i) for i in range(blocking_probes_order_format.shape[0])])
    blocking_probes_order_format = blocking_probes_order_format.assign(Amount = '25nm', Purification = 'STD')
    blocking_probes_order_format.loc[blocking_probes_order_format.length.values < 15, 'Amount'] = '100nm'
    blocking_probes_order_format = blocking_probes_order_format[['blocking_probe_name', 'sseq', 'Amount', 'Purification', 'length']]
    blocking_probes_order_format.to_excel('{}/{}_primerset_{}_barcode_selection_{}_full_length_blocking_probes_order_format.xlsx'.format(design_dir, design_dir_folder, primerset, barcode_selection))
    probe_length_filler = 'GTCTATTTTCTTATCCGACG'
    # 'GTCTATTTTCTTATCCGACGTGTTG'
    if primerset == 'A':
        forward_primer = 'CGCGGGCTATATGCGAACCG'
        reverse_primer = 'GCGTTGTATGCCCTCCACGC'
        # TAATACGACTCACTATAGGGCGTGGAGGGCATACAACGC
    elif primerset == 'B':
        forward_primer = 'CGATGCGCCAATTCCGGTTC'
        reverse_primer = 'CAACCCGCGAGCGATGATCA'
        # TAATACGACTCACTATAGGGTGATCATCGCTCGCGGGTTG
    elif primerset == 'C':
        forward_primer = 'GTTGGTCGGCACTTGGGTGC'
        reverse_primer = 'CCACCGGATGAACCGGCTTT'
        # TAATACGACTCACTATAGGGAAAGCCGGTTCATCCGGTGG
    else:
        forward_primer = ''
        reverse_primer = ''
    blocking_probe_seq_list = []
    for i in range(blocking_probes_order_format.shape[0]):
        bp_seq = blocking_probes_order_format.sseq.values[i]
        if len(bp_seq) == 15:
            probe_length_filler = 'GTCTATTTTCTTATCCGACGTGTTG'
        elif len(bp_seq) == 16:
            probe_length_filler = 'GTCTATTTTCTTATCCGACGTGTT'
        elif len(bp_seq) == 17:
            probe_length_filler = 'GTCTATTTTCTTATCCGACGTGT'
        elif len(bp_seq) == 18:
            probe_length_filler = 'GTCTATTTTCTTATCCGACGTG'
        elif len(bp_seq) == 19:
            probe_length_filler = 'GTCTATTTTCTTATCCGACGT'
        elif len(bp_seq) >= 20:
            probe_length_filler = 'GTCTATTTTCTTATCCGACG'
        blocking_probe = [[forward_primer + probe_length_filler + bp_seq + probe_length_filler + reverse_primer]]
        blocking_probe_seq_list = blocking_probe_seq_list + blocking_probe
    blocking_probe_seq_list = pd.DataFrame(blocking_probe_seq_list)
    blocking_probe_seq_list.columns = ['blocking_sequence']
    # print(blocking_probe_seq_list.shape[0])
    blocking_probe_seq_list['blocking_sequence'].str.upper().to_csv('{}/{}_primerset_{}_barcode_selection_{}_full_length_blocking_probe_sequences.txt'.format(design_dir, design_dir_folder, primerset, barcode_selection), header = False, index = False)
    return

# def write_final_probes_fasta(probes, design_dir, primerset, barcode_selection):
#     design_dir_folder = os.path.split(design_dir)[1]
#     probes_fasta_filename = '{}/{}_primerset_{}_barcode_selection_{}_full_length_probes.fasta'.format(design_dir, design_dir_folder, primerset, barcode_selection)
#     probes_list = [SeqRecord(Seq(probes['probe_full_seq'][i]), id = str(probes['probe_id'][i]), description = '') for i in range (0, probes.shape[0])]
#     SeqIO.write(probes_list, probes_fasta_filename, 'fasta')
#     return
# def write_final_probes_fasta(probes, design_dir, primerset, barcode_selection):
def write_final_probes_fasta(probes, probes_fasta_filename):
    # design_dir_folder = os.path.split(design_dir)[1]
    # probes_fasta_filename = '{}/full_length_probe_selection.fasta'.format(design_dir)
    # print(probes[['probe_id','Sequence']].iloc[0],)
    probes_list = [SeqRecord(Seq(probes['Sequence'].values[i]), id = str(probes['probe_id'].values[i]), description = '') for i in range (probes.shape[0])]
    print(probes_list)
    # probes_list = [SeqRecord(Seq(probes['Sequence'][i]), id = str(probes['probe_id'][i]), description = '') for i in range (0, probes.shape[0])]
    # print(probes_list)
    # SeqIO.write(probes_list, probes_fasta_filename, 'fasta')
    with open(probes_fasta_filename, 'w') as output_handle:
        SeqIO.write(probes_list, output_handle, 'fasta')

    # WRite non-full length probes to fasta
# probes_fasta_filename = '{}/probe_selection.fasta'.format(design_dir)
# probes_list = [SeqRecord(Seq(probes['seq'][i]), id = str(probes['Name'][i]), description = '') for i in range (0, probes.shape[0])]
# SeqIO.write(probes_list, probes_fasta_filename, 'fasta')

    return

def write_final_unique_probes_fasta(probes, design_dir, primerset, barcode_selection):
    design_dir_folder = os.path.split(design_dir)[1]
    probes = probes.drop_duplicates().reset_index().drop(['index'], axis = 1)
    probes_fasta_filename = '{}/{}_primerset_{}_barcode_selection_{}_full_length_probes_unique.fasta'.format(design_dir, design_dir_folder, primerset, barcode_selection)
    probes_list = [SeqRecord(Seq(probes['probe_full_seq'][i]), id = str(probes['probe_id'][i]), description = '') for i in range (0, probes.shape[0])]
    SeqIO.write(probes_list, probes_fasta_filename, 'fasta')
    return


def blast_probes(input_filename, blast_output_filename, blast_database, molecule):
    out_format = '6 qseqid sseqid pident qcovhsp length mismatch gapopen qstart qend sstart send evalue bitscore staxids qseq sseq'
    if molecule == 'RNA':
        # print('RNA')
        return_code = subprocess.call(['blastn', '-db', blast_database, '-query', input_filename, '-out', blast_output_filename, '-outfmt', out_format, '-task', 'blastn-short', '-max_hsps', '1', '-max_target_seqs', '100000', '-strand', 'minus', '-evalue', '100', '-num_threads', '1'])
    elif molecule == 'DNA':
        # print('DNA')
        return_code = subprocess.call(['blastn', '-db', blast_database, '-query', input_filename, '-out', blast_output_filename, '-outfmt', out_format, '-task', 'blastn-short', '-max_hsps', '1', '-max_target_seqs', '100000', '-strand', 'plus', '-evalue', '100', '-num_threads', '1'])
    return(return_code)


def blast_final_probes(design_dir, primerset, barcode_selection, blast_database, molecule):
# # Terminal
# blast_database = 'database/hao_zhou_mouse_gut/custom_blast_db.fasta'
# blast_database = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/databases/e_coli_mg1655/custom_blast_db.fasta'
# os.path.exists(blast_database)
    design_dir_folder = os.path.split(design_dir)[1]
    probes_fasta_filename = '{}/full_length_probe_selection.fasta'.format(design_dir)
    blast_output_filename = '{}/full_length_probe_selection.blast.out'.format(design_dir)
    out_format = '6 qseqid sseqid pident qcovhsp length mismatch gapopen qstart qend sstart send evalue bitscore staxids qseq sseq'
    if molecule == 'RNA':
        # print('molecule', molecule)
        return_code = subprocess.call(['blastn', '-db', blast_database, '-query', probes_fasta_filename, '-out', blast_output_filename, '-outfmt', out_format, '-task', 'blastn-short', '-max_hsps', '1', '-max_target_seqs', '100000', '-strand', 'minus', '-evalue', '100', '-num_threads', '1'])
    elif molecule == 'DNA':
        return_code = subprocess.call(['blastn', '-db', blast_database, '-query', probes_fasta_filename, '-out', blast_output_filename, '-outfmt', out_format, '-task', 'blastn-short', '-max_hsps', '1', '-max_target_seqs', '100000', '-strand', 'plus', '-evalue', '100', '-num_threads', '1'])
    # return_code = subprocess.call(['blastn', '-db', blast_database, '-query', probes_fasta_filename, '-out', blast_output_filename, '-outfmt', out_format, '-task', 'blastn-short', '-max_hsps', '1', '-max_target_seqs', '100000', '-strand', 'minus', '-evalue', '100', '-num_threads', '1'])
    # # non-full length
    # probes_fasta_filename = '{}/probe_selection.fasta'.format(design_dir)
    # blast_output_filename = '{}/probe_selection.blast.out'.format(design_dir)
    # out_format = '6 qseqid sseqid pident qcovhsp length mismatch gapopen qstart qend sstart send evalue bitscore staxids qseq sseq'
    # # return_code = subprocess.call(['blastn', '-db', blast_database, '-query', probes_fasta_filename, '-out', blast_output_filename, '-outfmt', out_format, '-task', 'blastn-short', '-max_hsps', '1', '-max_target_seqs', '100000', '-strand', 'plus', '-evalue', '100', '-num_threads', '1'])
    # return_code = subprocess.call(['blastn', '-db', blast_database, '-query', probes_fasta_filename, '-out', blast_output_filename, '-outfmt', out_format, '-task', 'blastn-short', '-max_hsps', '1', '-max_target_seqs', '100000', '-strand', 'minus', '-evalue', '100', '-num_threads', '1'])
# return_code = subprocess.check_call(['blastn', '-db', blast_database, '-query', infile, '-out', blast_output_filename, '-outfmt', out_format, '-task', 'blastn-short', '-max_hsps', '1', '-max_target_seqs', '100000', '-strand', 'plus', '-evalue', '100', '-num_threads', '1'])

    return(return_code)

def mch_test_v4(df):
    qseq_array = df.qseq.values
    sseq_array = df.sseq.values
    mch_array = np.zeros(len(qseq_array))
    for i in range(len(qseq_array)):
        qseq = qseq_array[i]
        sseq = sseq_array[i]
        if qseq != sseq:
            snp_indices = np.where(np.array(list(qseq)) != np.array(list(sseq)))[0]
            diffs = np.diff(snp_indices)
            mch_array[i] = np.max(np.append(diffs,[snp_indices[0], len(qseq) - 1 - snp_indices[-1]]))
        else:
            mch_array[i] = len(qseq)
    return(mch_array)

def calculate_tm(df, Na = 390, dnac1_oligo = 5):
    # get the query and subject sequences
    qseq_array = df.qseq.values
    sseq_array = df.sseq.values
    # Generate a new numpy array and iterate through the query values
    tm_array = np.zeros(len(qseq_array))
    for i in range(len(qseq_array)):
        qseq = qseq_array[i]
        # cseq = Seq(sseq_array[i]).complement()
        # Calculate the melting temp for the ith off target result using Bio.SeqUtils MeltingTemp
        tm_array[i] = mt.Tm_NN(qseq, Na = Na, saltcorr = 7, dnac1 = dnac1_oligo*15, dnac2 = 1)
    # Return the numpy array
    return(tm_array)

# Calculate the GC content of off target matches from blast results
def calculate_gc_count(df):
    # get the query and subject sequences
    qseq_array = df.qseq.values
    sseq_array = df.sseq.values
    # Generate a new numpy array and iterate through the query values
    gc_count_array = np.zeros(len(qseq_array), dtype = int)
    for i in range(len(qseq_array)):
        # Calculate the GC content using Bio.SeqUtils GC
        gc_count_array[i] = int(GC(qseq_array[i])*len(qseq_array[i])/100)
    # Return the numpy array
    return(gc_count_array)

# def summarize_final_probe_blast(df, mch):
#     df['mch'] = mch_test_v4(df)
#     df_filtered = df.loc[df.mch >= mch]
#     if df_filtered.shape[0] > 0:
#         probe_id = np.unique(df_filtered['probe_id'])[0]
#         probe_length = np.unique(df_filtered['probe_length'])[0]
#         check = np.sum(((df_filtered['probe_start'] >= 44) & (df_filtered['probe_end'] <= 44 + probe_length - 1) & (df_filtered['mch'] <= probe_length)) | (df_filtered['target_feature'] == df_filtered['feature']))/df_filtered.shape[0]
#     else:
#         print('{} has no blast hits...'.format(df.probe_id.values[0]))
#         check = 0
#     return(check)
def get_significant_blast_homologies(df, mch, mt_cutoff, ot_gc_cutoff, bitscore_thresh, dnaconc, sod):
    # # Terminal
    # df = probe_blast
    df['mch'] = mch_test_v4(df)
    df['melting_temp'] = calculate_tm(df, sod, dnaconc)
    df['GC_count'] = calculate_gc_count(df)
    df_filtered = df.loc[(df.mch > mch) | (df.melting_temp > mt_cutoff) | (df.GC_count > ot_gc_cutoff) | (df.bitscore > bitscore_thresh), :]
    if df_filtered.shape[0] == 0:
        print('{} has no blast hits...'.format(df.probe_id.values[0]))
    return(df_filtered)

# def check_final_probes_blast(design_dir, probes, mch, bot, primerset, barcode_selection):
#     # probes_filename = '{}/full_length_probes.csv'.format(design_dir)
#     design_dir_folder = os.path.split(design_dir)[1]
#     # Get the probe blast out filename
    # probes_blast_filename = '{}/full_length_probes_unique.blast.out'.format(design_dir)
    # os.path.exists(probes_blast_filename)
    # probes_blast = pd.read_table(probes_blast_filename, header = None)
    # probes_blast.columns = ['probe_id', 'feature', 'pid', 'qcovhsp', 'length', 'mismatch', 'gapopen', 'probe_start', 'probe_end', 'molecule_start', 'molecule_end', 'evalue', 'bitscore', 'staxids', 'qseq', 'sseq']
    # # probes = pd.read_csv(probes_filename, dtype = {'code': str})
    # probes['probe_length'] = probes['rna_seq'].apply(len) - 6
    # print(probes.shape)
    # probes_length_df = probes[['probe_id', 'target_feature', 'probe_length']]
    # probes_blast = probes_blast.merge(probes_length_df, on = 'probe_id', how = 'left')
    # probes_blast_summary = probes_blast.groupby('probe_id', axis = 0).apply(summarize_final_probe_blast, mch = mch)
    # probes_blast_summary = probes_blast_summary.reset_index()
    # probes_blast_summary.columns = ['probe_id', 'final_check']
    # probes_final_check = probes.merge(probes_blast_summary, on = 'probe_id', how = 'left')
    # print(probes_final_check.shape)
    # problematic_probes = probes_final_check.loc[(probes_final_check.final_check < bot)][['target_feature', 'probe_numeric_id', 'probe_id', 'final_check']].drop_duplicates()
    # print(problematic_probes.shape)
    # # for i in range(len(problematic_probes)):
    # #     probes_blast_problematic = probes_blast.loc[(probes_blast.probe_id.values == problematic_probes.probe_id.values[i]) & (probes_blast[target_rank].values != int(probes_blast.target_taxon.drop_duplicates()[0])),:]
    # #     probe_full_seq_modified = modify_problematic_probes(probes_blast_problematic, probes)
    # #     probes.loc[probes.probe_id.values == problematic_probes.probe_id.values[i], 'probe_full_seq'] = probe_full_seq_modified
    # probes_final_filter = probes_final_check.loc[~(probes_final_check.target_feature.isin(problematic_probes.target_feature.values) & probes_final_check.probe_numeric_id.isin(problematic_probes.probe_numeric_id.values))]
    # probes_final_filter.to_csv('{}/{}_primerset_{}_barcode_selection_{}_probes_final_check_filter_N.csv'.format(design_dir, design_dir_folder, primerset, barcode_selection))
    # print(probes_final_filter.shape)
    # probes_final_filter.loc[:,'probe_full_seq_length'] = probes_final_filter.probe_full_seq.apply(len)
    # probes_final_filter.to_csv('{}/{}_primerset_{}_barcode_selection_{}_full_length_probes.csv'.format(design_dir, design_dir_folder, primerset, barcode_selection), header = True, index = False)
    # probes_final_filter_summary = probes_final_filter['target_feature'].value_counts()
    # probes_final_filter_summary.to_csv('{}/{}_primerset_{}_barcode_selection_{}_probes_final_filter_summary.csv'.format(design_dir, design_dir_folder, primerset, barcode_selection))
    # probes_final_filter['probe_full_seq'].str.upper().to_csv('{}/{}_primerset_{}_barcode_selection_{}_full_length_probe_sequences.txt'.format(design_dir, design_dir_folder, primerset, barcode_selection), header = False, index = False)
    # probes_order_format = probes_final_filter[['probe_id', 'probe_full_seq']]
    # probes_order_format = probes_order_format.assign(Amount = '25nm', Purification = 'STD')
    # probes_order_format.loc[probes_order_format['probe_full_seq'].str.len() > 60]['Amount'] = '100nm'
    # probes_order_format.loc[probes_order_format['probe_full_seq'].str.len() > 60]['Purification'] = 'PAGE'
    # probes_order_format.to_excel('{}/{}_primerset_{}_barcode_selection_{}_full_length_probes_order_format.xlsx'.format(design_dir, design_dir_folder, primerset, barcode_selection))
#    return(probes_final_filter)

def get_blocking_probes_index(df):
    index = [*range(df.shape[0])]
    df['index'] = index
    return(df)


def check_final_probes_blast(design_dir, probe_selection, blast_output_filename, blast_database_filename, molecule,target_table_filename, total_transcript_table_filename, mch, bitscore_thresh, ot_gc_cutoff, mt_cutoff, bot, primerset, barcode_selection, dnaconc, sod, probe_selection_level):

# # Terminal
# design_dir = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/runs/2020_02_13_gfp_plasmid/simulation/DSGN_008'
# blast_database_filename = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/databases/e_coli_mg1655_coding/custom_blast_db.fasta'
# target_table_filename = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/runs/2020_02_13_gfp_plasmid/sample_008/inputs/target_table_008.csv'
# probe_selection = generate_full_probes(design_dir, target_table_filename, 0.99, plf = 0, primer = 0, primerset = 0, barcode_selection = 0)
# probe_selection[probe_selection['probe_id'] == 280].iloc[0,:]

    # probes_filename = '{}/full_length_probes.csv'.format(design_dir)
    design_dir_folder = os.path.split(design_dir)[1]

# # terminal
# probes = probe_selection
# mch = 14
# bitscore_thresh = 29
# ot_gc_cutoff = 7
# mt_cutoff = 60
# dnaconc = 5
# sod = 390
# bot = 0.99
# probe_selection_filename = '{}/'.format(design_dir)
# probe_selection[['seq', 'Sequence']]
# design_dir = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/runs/2020_02_13_gfp_plasmid/simulation/DSGN_003'

    # Get the probe blast out filename
    probes_blast_filename = blast_output_filename
    # print(probes_blast_filename)
    # probes_blast_filename = '{}/full_length_probe_selection.blast.out'.format(design_dir)
    # os.path.exists(probes_blast_filename)
    probes_blast = pd.read_csv(probes_blast_filename, header = None, delim_whitespace = True)
    # probes_blast = pd.read_csv(blast_output_filename, header = None, delim_whitespace = True)
    # probes_blast = pd.read_csv(probes_blast_filename, skiprows = 3, header = None, delim_whitespace = True)
    # probes_blast = pd.read_table(probes_blast_filename, skiprows = 3, header = None, delim_whitespace = True)
    probes_blast.columns = ['probe_id', 'feature', 'pid', 'qcovhsp', 'length', 'mismatch', 'gapopen', 'probe_start', 'probe_end', 'molecule_start', 'molecule_end', 'evalue', 'bitscore', 'staxids', 'qseq', 'sseq']
    # print('probes_blast',probes_blast.probe_id.unique())
    # print(probes_blast[probes_blast['probe_id'] == 681]['sseq'])
    # print('probes_blast',probes_blast)
    # probes_blast.loc[:,'random'] = 'gdts'
    # probes_blast.probe_id
    #
    # probe_id_to_feature = probe_selection[['Name','target_feature']]
    # probe_id_to_feature.loc[:,'Name'] = probe_id_to_feature.loc[:,'Name'].astype(str)
    # probe_id_to_feature.Name
    #
    # probe_blast = probe_blast.merge(probe_id_to_feature, left_on = 'probe_id', right_on = 'Name')
    # probes_blast_on_t = 0
    # probes_blast[probes_blast['feature'] == probe_id_to_feature['target_feature'][4]]
    # probes_blast_on_t.shape
    # probes_blast_on_t['target_feature']

    # probes = pd.read_csv(probes_filename, dtype = {'code': str})
    # probes['probe_length'] = probes['rna_seq'].apply(len) - 6
    # print(probes.shape)
    # probes_length_df = probes[['probe_id', 'target_feature', 'probe_length']]
    # probes_blast = probes_blast.merge(probes_length_df, on = 'probe_id', how = 'left')
        # Get all blast homologies that are above the thresholds
    # print(probes_blast[probes_blast['sseq'].apply(len) >= mch][['probe_id','sseq']])
    probes_blast_significant = get_significant_blast_homologies(probes_blast, mch = mch, mt_cutoff = mt_cutoff, ot_gc_cutoff = ot_gc_cutoff, bitscore_thresh = bitscore_thresh, dnaconc = dnaconc, sod = sod)
    # print(probes_blast_significant[['probe_id','sseq']])
    # print('probes_blast_significant',probes_blast_significant.probe_id.unique())
    # probe_blast_significant.shape
    probe_id_to_feature = probe_selection[['probe_id','target_feature']]
    # probe_id_to_feature.loc[:,'probe_id'] = probe_id_to_feature['probe_id'].astype(int)
    # probes_blast_significant.loc[:,'probe_id'] = probes_blast_significant['probe_id'].astype(int)
    #
    # probe_id_to_feature
    probes_blast_significant = probes_blast_significant.merge(probe_id_to_feature, on = 'probe_id')
    # print('probes_blast_significant',probes_blast_significant[['probe_id','feature','target_feature']])
    # probes_blast_significant_ot = pd.DataFrame()
    # for index, probe in probe_selection.iterrows():
    #     feature = probe['target_feature']
    #     target_table = pd.read_csv(target_table_filename)
    #     # probe_selection_level = target_table[target_table['target_feature'] == feature]['selection_level'].values[0]
    #     # print(index,'of',probe_selection.shape[0])
    #     if probe_selection_level == 'genome_bin':
    #         total_transcript_table = pd.read_csv(total_transcript_table_filename)
    #         feature_transcript_info = total_transcript_table[total_transcript_table['Gene'] == str(feature)]
    #         phylogenetic_marker_genes = ['dnaG', 'frr', 'infC', 'nusA', 'pgk', 'pyrG', 'rplA', 'rplB', 'rplC', 'rplD', 'rplE', 'rplF', 'rplK', 'rplL', 'rplM', 'rplN',
    #         'rplP', 'rplS', 'rplT', 'rpmA', 'rpoB', 'rpsB', 'rpsC', 'rpsE', 'rpsI', 'rpsJ', 'rpsK', 'rpsM', 'rpsS', 'smpB', 'tsf']
    #         feature_genome_bin = total_transcript_table[total_transcript_table['Genome_bin'] == feature_transcript_info['Genome_bin']]
    #         feature_genome_bin_marker = feature_genome_bin[feature_genome_bin['Preferred_name'].isin(phylogenetic_marker_genes)]
    #         genome_bin_feature_list = feature_genome_bin_marker['Gene'].values
    #         probes_blast_significant_ot_temp =  probe_blast[~probe_blast.feature.isin(genome_bin_feature_list)]
    #         probes_blast_significant_ot = probes_blast_significant_ot.append(probes_blast_significant_ot_temp)
    #     elif probe_selection_level == 'molecule':
    probes_blast_significant_ot = probes_blast_significant[probes_blast_significant['target_feature'] != probes_blast_significant['feature']]
    print('probe_blast_signficant on target final blast',probes_blast_significant[probes_blast_significant['target_feature'] == probes_blast_significant['feature']][['probe_id','sseq']])
    print('probes_blast_significant_ot',probes_blast_significant_ot[['probe_id','sseq']])
    # print('probes_blast_significant_ot',probes_blast_significant_ot.shape[0])
    #         probes_blast_significant_ot_temp = probes_blast_significant[probes_blast_significant['target_feature'] != probes_blast_significant['feature']]
    #         probes_blast_significant_ot = probes_blast_significant_ot.append(probes_blast_significant_ot_temp)
            # print(probe_blast_significant.shape[0])
    # probes_blast_significant_ot.iloc[0,:]
    # probes_blast_significant_ot[probes_blast_significant_ot['probe_id'] == 280].iloc[0,:]
    # probes_blast_significant_ot[probes_blast_significant_ot['probe_id'] == 280].shape
    # probes_blast_significant_ot.shape

    # probe_blast_significant.shape
    # probe_blast_significant_ot.shape
    try:
        blocking_probes_filename = glob.glob('{}/*_blocking_probes.csv'.format(design_dir))
        blocking_probes_list = [pd.read_csv(i) for i in blocking_probes_filename]
        blocking_probes = pd.concat(blocking_probes_list)
        blocking_probes = blocking_probes.merge(probe_id_to_feature, on = 'probe_id')
    except:
        blocking_probes = pd.DataFrame(columns = probes_blast_significant_ot.columns)
    # blocking_probes[blocking_probes['probe_id'] == 280].iloc[0,:]
    # blocking_probes[blocking_probes['probe_id'] == 280].shape
    # blocking_probes.shape
    # blocking_probes.shape
    # blocking_probes.columns
    # print('probes_blast_significant_ot',probes_blast_significant_ot.shape[0])
    probes_blast_significant_ot_unblocked = probes_blast_significant_ot[~probes_blast_significant_ot['sseq'].isin(blocking_probes['sseq'])]
    # print('blocking_probes',blocking_probes)
    probes_blast_significant_ot_unblocked[['probe_id','sseq']]
    print('There are {} probes not blocked.'.format(probes_blast_significant_ot_unblocked.shape[0]))
    # probes_blast_significant_ot_filename = '{}/full_length_off_target.csv'.format(design_dir)
    # probes_blast_significant_ot.to_csv(probes_blast_significant_ot_filename)
    # probes_blast_significant_ot_unblocked_filename = '{}/full_length_unblocked.csv'.format(design_dir)
    # probes_blast_significant_ot_unblocked.to_csv(probes_blast_significant_ot_unblocked_filename)
    # probe_blast_significant_ot_unblocked.shape

    # Concatenate the blocking probes with the unblocked off-target
    blocking_probes_all = pd.concat([blocking_probes, probes_blast_significant_ot_unblocked])
    target_table = pd.read_csv(target_table_filename)
    feature_to_name = target_table[['target_feature','target_name']]
    # Add main probe info
    blocking_probes_all = blocking_probes_all.merge(probe_selection[['p_start','probe_id','ln','seq']], on = 'probe_id', how = 'left')
    blocking_probes_all = blocking_probes_all.merge(feature_to_name, on = 'target_feature', how = 'left')
    blocking_probes_all = blocking_probes_all.groupby('probe_id', axis = 0).apply(get_blocking_probes_index)
    # Generate a name column that includes  target, and probe location on the target molecule info
    blocking_probes_all['Name'] = 'Block.' + blocking_probes_all['target_name'] + '_' + blocking_probes_all['p_start'].astype(int).astype(str)  + '_' + blocking_probes_all['index'].astype(int).astype(str)
    # print(blocking_probes_all[['Name','sseq']])
    # # Create a DataFrame for the flanking regions
    # flanking_regions = pd.DataFrame({'readout':['F3','F7'], 'readout_sequence':['GGGAGAATGAGGTGTAATGT','TTGGAGGTGTAGGGAGTAAA']})
    # # Merge the flanking regions with the probe_selection
    # blocking_probes_all = blocking_probes_all.merge(flanking_regions, left_on = 'readout_1', right_on = 'readout', how = 'left')
    # Create new columns for the IDT excel output
    # blocking_probes_all['Sequence'] = blocking_probes_all['sseq']
    # blocking_probes_all.loc[:, 'Scale'] = '25nm'
    blocking_probes_all.loc[:, 'Purification'] = 'STD'

    # blocking_probes_all = blocking_probes_all.merge(probe_selection[['probe_id','ln']], how = 'probe_id', on = 'left')
    idt_table = pd.DataFrame()
    blast_database = SeqIO.to_dict(SeqIO.parse(blast_database_filename, "fasta"))
# Check if the length is greater than 60
    for index, probe in blocking_probes_all.iterrows():
        # print(probe)
        # probe = blocking_probes_all[blocking_probes_all['sseq'] == seq]
        seq = probe.sseq
        record = re.sub('ref','',probe['feature'])
        if record[0] == '|':
            record = record[1:-1]
        off_target_full_sequence = blast_database[record]
        if len(seq) > 60:
            # blocking_probes_all.loc[blocking_probes_all.Sequence.values == seq, 'Scale'] = '100nm'
            probe['Scale'] = '100nm'
            probe['Sequence'] = seq
        elif len(seq) < 20:
            # off_target_full_sequence = blast_database[re.sub('[|ref]','',probe['feature'].values[0])]
            # # blocking_probes_all.loc[blocking_probes_all['sseq'] == seq, 'Sequence'] = off_target_full_sequence[probe['molecule_start'] : probe['molecule_start'] + 20]
            # print(probe.index.values[0])
            # print(blocking_probes_all.at[probe.index.values[0], 'Sequence'])
            # print(probe['molecule_start'].values[0])
            # blocking_probes_all.replace(probe['Sequence'].values[0], str(off_target_full_sequence[probe['molecule_start'][0] : probe['molecule_start'][0] + 20].seq))
            # blocking_probes_all.at[probe.index.values[0], 'Sequence'] = str(off_target_full_sequence[probe['molecule_start'].values[0] : probe['molecule_start'].values[0] + 20].seq)
            probe['Scale'] = '25nm'
            print('-\nmolecule\n',probe['molecule_start'] ,'\n', probe['molecule_end'],'\n',str(off_target_full_sequence[probe['molecule_end'] : probe['molecule_start']].seq.reverse_complement()))
            if probe['molecule_start'] < probe['molecule_end']:
                if molecule == 'RNA':
                    probe['Sequence'] = str(off_target_full_sequence[probe['molecule_start'] - 1 : probe['molecule_start'] + 19].seq.reverse_complement())
                elif molecule == 'DNA':
                    probe['Sequence'] = str(off_target_full_sequence[probe['molecule_start'] - 1 : probe['molecule_start'] + 19].seq)
                # print(1)
            else:
                if molecule == 'RNA':
                    probe['Sequence'] = str(off_target_full_sequence[probe['molecule_end'] - 1 : probe['molecule_end'] + 19].seq.reverse_complement())
                elif molecule == 'DNA':
                    probe['Sequence'] = str(off_target_full_sequence[probe['molecule_end'] - 1 : probe['molecule_end'] + 19].seq)
                # print(2)
            # print(probe[['sseq','Sequence']],  probe['molecule_start'], probe['molecule_end'])
            # print(probe['Name'], probe['molecule_start'], probe['molecule_end'], '\n', probe['sseq'], '\n', probe['Sequence'], '\n', probe['seq'])
        else:
            # off_target_full_sequence = blast_database[re.sub('ref|','',probe['feature'].values[0])]
            probe['Scale'] = '25nm'
            if probe['molecule_start'] < probe['molecule_end']:
                if molecule == 'RNA':
                    probe['Sequence'] = str(off_target_full_sequence[probe['molecule_start'] - 1 : probe['molecule_end'] + 1].seq.reverse_complement())
                elif molecule == 'DNA':
                    probe['Sequence'] = str(off_target_full_sequence[probe['molecule_start'] - 1 : probe['molecule_end'] + 1].seq)
                # print(3)
            else:
                if molecule == 'RNA':
                    probe['Sequence'] = str(off_target_full_sequence[probe['molecule_end'] - 1 : probe['molecule_start'] + 1].seq.reverse_complement())
                elif molecule == 'DNA':
                    probe['Sequence'] = str(off_target_full_sequence[probe['molecule_end'] - 1 : probe['molecule_start'] + 1].seq)
            # print(probe['Name'], probe['molecule_start'], probe['molecule_end'], '\n', probe['sseq'], '\n',probe['Sequence'], '\n', probe['seq'])
                # print(4)
            # probe['Sequence'] = str(off_target_full_sequence[probe['molecule_start'] : probe['molecule_end'] + 2].seq)
            # print(probe[['sseq','Sequence']], probe['molecule_start'], probe['molecule_end'])
            # blocking_probes_all.at[probe.index.values[0], 'Sequence'] = str(off_target_full_sequence[probe['molecule_start'].values[0] : probe['molecule_end'].values[0] + 2].seq)
        # temp['Name'] = probe['Name']
        print(probe['Sequence'])
        idt_table = idt_table.append(probe[['Name','Sequence','Scale']])
    blocking_probes_all = blocking_probes_all.merge(idt_table, on = 'Name', how = 'left')
    print('blocking_probes_all_dup', blocking_probes_all.shape[0])
    blocking_probes_all.drop_duplicates(subset = 'Sequence', keep = 'first', inplace = True)
    print('There are {} blocking probes.'.format(blocking_probes_all.shape[0]))
    blocking_probes_all_filename = '{}/blocking_probes_all.csv'.format(design_dir)
    blocking_probes_all.to_csv(blocking_probes_all_filename)

    #
    # # Short probes
    # short_probes_blast_filename = '{}/probe_selection.blast.out'.format(design_dir)
    # os.path.exists(short_probes_blast_filename)
    # short_probes_blast = pd.read_table(short_probes_blast_filename, header = None)
    # short_probes_blast.columns = ['probe_id', 'feature', 'pid', 'qcovhsp', 'length', 'mismatch', 'gapopen', 'probe_start', 'probe_end', 'molecule_start', 'molecule_end', 'evalue', 'bitscore', 'staxids', 'qseq', 'sseq']
    # short_probes_blast.shape
    # short_probe_blast_significant = get_significant_blast_homologies(short_probes_blast, mch = mch, mt_cutoff = mt_cutoff, ot_gc_cutoff = ot_gc_cutoff, bitscore_thresh = bitscore_thresh, dnaconc = dnaconc, sod = sod)
    # short_probe_blast_significant.shape
    # probe_id_to_feature = probe_selection[['Name','target_feature']]
    # probe_id_to_feature
    # short_probe_blast_significant = short_probe_blast_significant.merge(probe_id_to_feature, left_on = 'probe_id', right_on = 'Name')
    # short_probe_blast_significant_ot = short_probe_blast_significant[short_probe_blast_significant['target_feature'] != short_probe_blast_significant['feature']]
    # short_probe_blast_significant_on_t = short_probe_blast_significant[sum(short_probe_blast_significant['target_feature'] == short_probe_blast_significant['feature'])]
    # short_probe_blast_significant.shape
    # short_probe_blast_significant_ot.shape
    # short_probe_blast_significant_ot_unblocked = short_probe_blast_significant_ot[~short_probe_blast_significant_ot['sseq'].isin(blocking_probes['sseq'])]
    # short_probe_blast_significant_ot_unblocked.shape
    #
# blast_database = SeqIO.to_dict(SeqIO.parse(blast_database_filename, "fasta"))
# e_coli_genome = blast_database["NC_000913.3"]
# e_coli_genome = blast_database[re.sub('[|ref]','',blocking_probes_all[blocking_probes_all['sseq'] == seq]['feature'].values[0])]
# e_coli_genome.seq[2013741:2013752]

    # record_dict["NODE_3905_length_30719_cov_9.533818_22"]
    # record_dict["NODE_3527_length_32813_cov_15.070761_16"]


    # # probes_blast_summary = probes_blast.groupby('probe_id', axis = 0).apply(get_significant_blast_homologies, mch = mch, mt_cutoff = mt_cutoff, ot_gc_cutoff = ot_gc_cutoff, bitscore_thresh = bitscore_thresh, dnaconc = dnaconc, sod = sod)
    # probes_blast_summary.shape
    # probes_blast_summary = probes_blast_summary.reset_index()
    # probes_blast_summary.columns = ['probe_id', 'final_check']
    # probes_final_check = probes.merge(probes_blast_summary, on = 'probe_id', how = 'left')
    # print(probes_final_check.shape)
    # problematic_probes = probes_final_check.loc[(probes_final_check.final_check < bot)][['target_feature', 'probe_numeric_id', 'probe_id', 'final_check']].drop_duplicates()
    # print(problematic_probes.shape)
    # # for i in range(len(problematic_probes)):
    # #     probes_blast_problematic = probes_blast.loc[(probes_blast.probe_id.values == problematic_probes.probe_id.values[i]) & (probes_blast[target_rank].values != int(probes_blast.target_taxon.drop_duplicates()[0])),:]
    # #     probe_full_seq_modified = modify_problematic_probes(probes_blast_problematic, probes)
    # #     probes.loc[probes.probe_id.values == problematic_probes.probe_id.values[i], 'probe_full_seq'] = probe_full_seq_modified
    # probes_final_filter = probes_final_check.loc[~(probes_final_check.target_feature.isin(problematic_probes.target_feature.values) & probes_final_check.probe_numeric_id.isin(problematic_probes.probe_numeric_id.values))]
    # probes_final_filter.to_csv('{}/{}_primerset_{}_barcode_selection_{}_probes_final_check_filter_N.csv'.format(design_dir, design_dir_folder, primerset, barcode_selection))
    # print(probes_final_filter.shape)
    # probes_final_filter.loc[:,'probe_full_seq_length'] = probes_final_filter.probe_full_seq.apply(len)
    # probes_final_filter.to_csv('{}/{}_primerset_{}_barcode_selection_{}_full_length_probes.csv'.format(design_dir, design_dir_folder, primerset, barcode_selection), header = True, index = False)
    # probes_final_filter_summary = probes_final_filter['target_feature'].value_counts()
    # probes_final_filter_summary.to_csv('{}/{}_primerset_{}_barcode_selection_{}_probes_final_filter_summary.csv'.format(design_dir, design_dir_folder, primerset, barcode_selection))
    # probes_final_filter['probe_full_seq'].str.upper().to_csv('{}/{}_primerset_{}_barcode_selection_{}_full_length_probe_sequences.txt'.format(design_dir, design_dir_folder, primerset, barcode_selection), header = False, index = False)
    # probes_order_format = probes_final_filter[['probe_id', 'probe_full_seq']]
    # probes_order_format = probes_order_format.assign(Amount = '25nm', Purification = 'STD')
    # probes_order_format.loc[probes_order_format['probe_full_seq'].str.len() > 60]['Amount'] = '100nm'
    # probes_order_format.loc[probes_order_format['probe_full_seq'].str.len() > 60]['Purification'] = 'PAGE'
    # probes_order_format.to_excel('{}/{}_primerset_{}_barcode_selection_{}_full_length_probes_order_format.xlsx'.format(design_dir, design_dir_folder, primerset, barcode_selection))
    return(blocking_probes_all)



def generate_probe_statistics_plots(design_dir, probes_final_filter, plasmid_annotation_filtered, primerset, barcode_selection, theme_color):
    design_dir_folder = os.path.split(design_dir)[1]
    fig = plt.figure()
    fig.set_size_inches(cm_to_inches(5),cm_to_inches(4))
    plt.hist(probes_final_filter.Tm.values, bins = 20, color = (0,0.5,1), histtype = 'step')
    plt.xlabel(r'Melting Temperature [$^\circ$C]', fontsize = 8, color = theme_color)
    plt.ylabel('Frequency', fontsize = 8, color = theme_color)
    plt.tick_params(direction = 'in', length = 2, colors = theme_color)
    plt.axes().spines['left'].set_color(theme_color)
    plt.axes().spines['bottom'].set_color(theme_color)
    plt.axes().spines['top'].set_color(theme_color)
    plt.axes().spines['right'].set_color(theme_color)
    plt.subplots_adjust(left = 0.3, bottom = 0.22, top = 0.98, right = 0.98)
    probes_tm_histogram_filename = '{}/{}_primerset_{}_barcode_selection_{}_full_length_probes_tm_histogram.pdf'.format(design_dir, design_dir_folder, primerset, barcode_selection)
    fig.savefig(probes_tm_histogram_filename, dpi = 300, transparent = True)
    plt.close()
    fig = plt.figure()
    fig.set_size_inches(cm_to_inches(5),cm_to_inches(4))
    plt.hist(probes_final_filter.GC.values, bins = 20, color = (0,0.5,1), histtype = 'step')
    plt.xlabel('GC Content [%]', fontsize = 8, color = theme_color)
    plt.ylabel('Frequency', fontsize = 8, color = theme_color)
    plt.tick_params(direction = 'in', length = 2, colors = theme_color)
    plt.axes().spines['left'].set_color(theme_color)
    plt.axes().spines['bottom'].set_color(theme_color)
    plt.axes().spines['top'].set_color(theme_color)
    plt.axes().spines['right'].set_color(theme_color)
    plt.subplots_adjust(left = 0.3, bottom = 0.22, top = 0.98, right = 0.98)
    probes_gc_histogram_filename = '{}/{}_primerset_{}_barcode_selection_{}_full_length_probes_gc_histogram.pdf'.format(design_dir, design_dir_folder, primerset, barcode_selection)
    fig.savefig(probes_gc_histogram_filename, dpi = 300, transparent = True)
    plt.close()
    fig = plt.figure()
    fig.set_size_inches(cm_to_inches(5), cm_to_inches(4))
    plt.hist(probes_final_filter.probe_full_seq_length.values, color = (0,0.5,1), histtype = 'step')
    plt.xlabel('Probe Length [bp]', fontsize = 8, color = theme_color)
    plt.ylabel('Frequency', fontsize = 8, color = theme_color)
    plt.tick_params(direction = 'in', length = 2, colors = theme_color)
    plt.axes().spines['left'].set_color(theme_color)
    plt.axes().spines['bottom'].set_color(theme_color)
    plt.axes().spines['top'].set_color(theme_color)
    plt.axes().spines['right'].set_color(theme_color)
    plt.subplots_adjust(left = 0.3, bottom = 0.22, top = 0.98, right = 0.98)
    probe_length_histogram_filename = '{}/{}_primerset_{}_barcode_selection_{}_full_length_probes_length_histogram.pdf'.format(design_dir, design_dir_folder, primerset, barcode_selection)
    fig.savefig(probe_length_histogram_filename, dpi = 300, transparent = True)
    plt.close()
    probes_final_filter_summary = probes_final_filter['target_feature'].value_counts().reset_index()
    probes_final_filter_summary.columns = ['FEATURE', 'PROBE_COUNTS']
    plasmid_annotation_filtered['FEATURE'] = ['feature_{}'.format(f) for f in plasmid_annotation_filtered.Feature.values]
    probes_final_filter_summary = probes_final_filter_summary.merge(plasmid_annotation_filtered, on = 'FEATURE', how = 'left')
    fig = plt.figure()
    fig.set_size_inches(cm_to_inches(5),cm_to_inches(4))
    plt.plot(probes_final_filter_summary.Length.values, probes_final_filter_summary.PROBE_COUNTS.values, 'o', color = (0,0.5,1), alpha = 0.8)
    plt.xlabel('Feature Length [bp]', fontsize = 8, color = theme_color)
    plt.ylabel('Probe Counts', fontsize = 8, color = theme_color)
    plt.tick_params(direction = 'in', length = 2, colors = theme_color)
    plt.axes().spines['left'].set_color(theme_color)
    plt.axes().spines['bottom'].set_color(theme_color)
    plt.axes().spines['top'].set_color(theme_color)
    plt.axes().spines['right'].set_color(theme_color)
    plt.subplots_adjust(left = 0.3, bottom = 0.22, top = 0.98, right = 0.98)
    counts_v_length_filename = '{}/{}_primerset_{}_barcode_selection_{}_probe_counts_vs_feature_length.pdf'.format(design_dir, design_dir_folder, primerset, barcode_selection)
    plt.savefig(counts_v_length_filename, dpi = 300, transparent = True)
    return

def generate_probe_summary_file(design_dir, probes_final_filter, primerset, barcode_selection):
    design_dir_folder = os.path.split(design_dir)[1]
    probes_summary_filename = '{}/full_length_probes_summary.txt'.format(design_dir)
    # probes_summary_filename = '{}/{}_full_length_probes_summary.txt'.format(design_dir, design_dir_folder, primerset, barcode_selection)
    with open(probes_summary_filename, 'w') as tf:
        tf.write('Probe design complete.')
    return

def get_helper_probes(design_dir, probe_selection, target_table_filename):
    helper_probes_filename = glob.glob('{}/*_helper_probes.csv'.format(design_dir))
    helper_probes_list = [pd.read_csv(i) for i in helper_probes_filename]
    helper_probes = pd.concat(helper_probes_list)
    target_table = pd.read_csv(target_table_filename)
    helper_probes = helper_probes.merge(target_table, on = 'target_feature', how = 'left')
    helper_probes['Name'] = 'Help.' + helper_probes['target_name'] + '_' + helper_probes['p_start'].astype(int).astype(str)
    helper_probes['Sequence'] = helper_probes['seq']
    helper_probes.loc[:, 'Scale'] = '25nm'
    helper_probes.loc[:, 'Purification'] = 'STD'
    # Check if the length is greater than 60
    for seq in helper_probes.Sequence.values:
        if len(seq) > 60:
            helper_probes.loc[helper_probes.Sequence.values == seq, 'Scale'] = '100nm'
    return(helper_probes)


###############################################################################################################
# main function
###############################################################################################################

def main():
    parser = argparse.ArgumentParser('Design FISH probes for a complex microbial community')
    # input blast filename
    parser.add_argument('design_dir', type = str, help = 'Input file containing blast results')

    parser.add_argument('sim_primer3_dir', type = str, help = 'Input file containing blast results')

    parser.add_argument('target_table_filename', type = str, help = 'Directory containing inputs')

    parser.add_argument('blast_database', type = str, help = 'Input file containing blast lineage')

    parser.add_argument('bot', type = float, help = 'Input file containing blast lineage')

    parser.add_argument('mch', type = int, help = 'Input file containing blast lineage')

    parser.add_argument('bplc', type = int, help = 'Input file containing blast lineage')

    parser.add_argument('-ps', '--primerset', dest = 'primerset', default = 'B', type = str, help = 'Input file containing blast lineage')
    parser.add_argument('-plf', '--plf', dest = 'plf', default = 'T', type = str, help = 'Input file containing blast lineage')
    parser.add_argument('-p', '--primer', dest = 'primer', default = 'T', type = str, help = 'Input file containing blast lineage')
    parser.add_argument('-bs', '--barcode_selection', dest = 'barcode_selection', default = 'MostSimple', type = str, help = 'Input file containing blast lineage')
    parser.add_argument('-tc', '--theme_color', dest = 'theme_color', default = 'black', help = 'Input file containing blast lineage')

    parser.add_argument('-bt', '--bitscore_thresh', dest = 'bitscore_thresh', type = float, default = 27, help = 'Number of top probes to keep')

    parser.add_argument('-sod', '--sod', dest = 'sod', type = float, default = 390, help = 'sodium concentration in nM')

    parser.add_argument('-dnaconc', '--dnaconc', dest = 'dnaconc', type = float, default = 5, help = 'oligo concentration in nM')

    parser.add_argument('-mt', '--mt_cutoff', dest = 'mt_cutoff', type = float, default = 60, help = 'oligo concentration in nM')

    parser.add_argument('-otgc', '--ot_gc_cutoff', dest = 'ot_gc_cutoff', type = float, default = 7, help = 'oligo concentration in nM')

    parser.add_argument('-ttt', '--total_transcript_table_filename', dest = 'total_transcript_table_filename', type = str, default = '', help = 'oligo concentration in nM')

    parser.add_argument('-mol', '--molecule', dest = 'molecule', type = str, default = 'RNA', help = 'Is the target DNA or RNA (i.e. are we hybridizing to non-coding or coding)?')

    parser.add_argument('-psl', '--probe_selection_level', dest = 'probe_selection_level', type = str, default = 'molecule', help = 'what level are we selecting probes on?')


    args = parser.parse_args()

    print('Generating full probes for design {}...'.format(os.path.basename(args.design_dir)))
    # Get all the feature best probes filenames
    feature_best_probes = glob.glob('{}/*_probe_selection.csv'.format(args.design_dir))
    # # Iterate through the filenames
    # for file in feature_best_probes:
    #     feature_best_probes_sa_filename = re.sub('_probe_selection.csv', '_probe_selection_sa.csv', file)
    #     add_spacer(file, args.sim_primer3_dir, feature_best_probes_sa_filename)
    # plasmid_annotation_filtered = pd.read_csv('{}/plasmid_annotation_filtered.csv'.format(args.sim_primer3_dir))
    # if not os.path.exists
    probe_selection = generate_full_probes(args.design_dir, args.target_table_filename, args.blast_database, args.mch, args.bitscore_thresh, args.ot_gc_cutoff, args.mt_cutoff, args.dnaconc, args.sod, args.molecule, args.probe_selection_level)
    # print('probe_selection\n',probe_selection[['Name','Sequence']])
    probes_fasta_filename = '{}/full_length_probe_selection.fasta'.format(args.design_dir)
    # write_final_probes_fasta(probe_selection, args.design_dir, args.primerset, args.barcode_selection)
    write_final_probes_fasta(probe_selection, probes_fasta_filename)
# write_final_probes_fasta(probe_selection, design_dir, 0,0)
    # write_final_unique_probes_fasta(oligo_df, args.design_dir, args.primerset, args.barcode_selection)
    blast_output_filename =  re.sub('.fasta','.blast.out',probes_fasta_filename)
    blast_probes(probes_fasta_filename, blast_output_filename, args.blast_database, args.molecule)
    # blast_final_probes(args.design_dir, args.primerset, args.barcode_selection, args.blast_database, args.molecule)
# blast_final_probes(design_dir, '','',blast_database)
    blocking_probes = check_final_probes_blast(args.design_dir, probe_selection, blast_output_filename, args.blast_database, args.molecule, args.target_table_filename, args.total_transcript_table_filename,args.mch, args.bitscore_thresh, args.ot_gc_cutoff, args.mt_cutoff, args.bot, args.primerset, args.barcode_selection, args.dnaconc, args.sod, args.probe_selection_level)
    helper_probes = get_helper_probes(args.design_dir, probe_selection, args.target_table_filename)
    # generate_blocking_probes(args.design_dir, args.bplc, plf = args.plf, primer = args.primer, primerset = args.primerset, barcode_selection = args.barcode_selection)
    # # generate_probe_statistics_plots(args.design_dir, probes_final_filter, plasmid_annotation_filtered, args.primerset, args.barcode_selection, args.theme_color)
    generate_probe_summary_file(args.design_dir, probe_selection, args.primerset, args.barcode_selection)
    probe_selection_order_format = pd.concat([probe_selection[['Name','Sequence','Scale','Purification']], blocking_probes[['Name','Sequence','Scale','Purification']], helper_probes[['Name','Sequence','Scale','Purification']]])
    design_dir_folder = os.path.split(args.design_dir)[1]
    probe_selection_all = probe_selection.to_csv('{}/{}_probe_selection_all.csv'.format(args.design_dir, design_dir_folder))
    probe_selection_order_format.to_excel('{}/{}_full_length_probe_selection_order_format.xlsx'.format(args.design_dir, design_dir_folder), index = False)
    return

if __name__ == '__main__':
    main()
