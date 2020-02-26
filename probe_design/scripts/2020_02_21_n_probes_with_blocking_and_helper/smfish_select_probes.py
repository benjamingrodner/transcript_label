import argparse
import threading
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
import tables
pd.options.mode.chained_assignment = None

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

###############################################################################################################
# HiPR-FISH Probe Design Pipeline
###############################################################################################################

###############################################################################################################
# Workflow functions
###############################################################################################################
# def select_all_specific_probes(probe_summary_info, bot, bitscore_thresh):
def select_all_specific_probes(probe_summary_info, bot, bitscore_thresh, probes_merge_filename):
    # Read csv for all probes summary (for debugging)
    probe_summary_info = pd.read_csv(probes_merge_filename)
    if np.max(probe_summary_info['blast_on_target_rate']) > bot:
        best_probes = probe_summary_info[(probe_summary_info['blast_on_target_rate'] > bot) & (probe_summary_info['off_target_max_bitscore'] < bitscore_thresh)]
        best_probes['selection_method'] = 'AllSpecific'
        print('There are {} best probes'.format(best_probes.shape[0]))
    else:
        best_probes = probe_summary_info.iloc[[0],:]
        best_probes.loc[:,'selection_method'] = 'AllSpecificSingleBest'
        print('There are no probes that fit the parameters. Selecting the single best.')
    return(best_probes)

# def select_all_specific_p_start_group_probes(probe_summary_info, group_distance, bot, bitscore_thresh):
#     best_probes_group = pd.DataFrame()
#     if np.max(probe_summary_info['blast_on_target_rate']) > bot:
#         best_probes = probe_summary_info.loc[(probe_summary_info['blast_on_target_rate'] > bot) & (probe_summary_info['off_target_max_bitscore'] < bitscore_thresh)]
#         for group in range(int(np.floor(1500/group_distance))):
#             best_probes_temp = best_probes.loc[best_probes.p_start_group.values == group]
#             best_probes_temp.sort_values(['taxon_coverage', 'off_target_full_qcovhsp_fraction', 'off_target_max_mch', 'off_target_max_bitscore', 'on_target_full_match', 'quality'], ascending = [False, True, True, True, False, True], inplace = True)
#             if not best_probes_temp.empty:
#                 best_probes_group = best_probes_group.append(best_probes_temp.iloc[[0],:], sort = True)
#                 best_probes_group.loc[:,'selection_method'] = 'AllSpecificPStartGroup'
#     else:
#         probe_summary_info.sort_values(['taxon_coverage', 'off_target_full_qcovhsp_fraction', 'off_target_max_mch', 'off_target_max_bitscore', 'on_target_full_match', 'quality'], ascending = [False, True, True, True, False, True], inplace = True)
#         best_probes_group = probe_summary_info.iloc[[0],:]
#         best_probes_group.loc[:,'selection_method'] = 'AllSpecificPStartGroupSingleBest'
#     return(best_probes_group)
def select_all_specific_p_start_group_probes(probe_summary_info, group_distance, bot, bitscore_thresh, mt_cutoff, ot_gc_cutoff):
    best_probes_group = pd.DataFrame()
    best_probes = probe_summary_info.loc[(probe_summary_info['blast_on_target_rate'] > bot) & (probe_summary_info['off_target_max_bitscore'] < bitscore_thresh) & (probe_summary_info['off_target_max_tm'] < mt_cutoff) & (probe_summary_info['off_target_max_gc'] < ot_gc_cutoff)]
    if not best_probes.empty:
        for group in range(int(np.floor(1500/group_distance))):
            best_probes_temp = best_probes.loc[best_probes.mean_probe_start_group.values == group]
            best_probes_temp.sort_values(['taxon_coverage', 'off_target_full_qcovhsp_fraction', 'off_target_max_mch', 'off_target_max_bitscore', 'off_target_max_tm', 'off_target_max_gc', 'on_target_full_match', 'quality'], ascending = [False, True, True, True, True, True, False, True], inplace = True)
            if not best_probes_temp.empty:
                best_probes_group = best_probes_group.append(best_probes_temp.iloc[[0],:], sort = False)
                best_probes_group.loc[:,'selection_method'] = 'AllSpecificPStartGroup'
    else:
        probe_summary_info.sort_values(['blast_on_target_rate', 'taxon_coverage', 'off_target_full_qcovhsp_fraction', 'off_target_max_mch', 'off_target_max_bitscore', 'off_target_max_tm', 'off_target_max_gc', 'on_target_full_match', 'quality'], ascending = [False, False, True, True, True, True, True, False, True], inplace = True)
        best_probes_group = probe_summary_info.iloc[[0],:]
        best_probes_group.loc[:,'selection_method'] = 'AllSpecificPStartGroupSingleBest'
    return(best_probes_group)

def select_min_overlap_probes(probe_summary_info, taxon_fasta_filename, probe_evaluation_filename, max_continuous_homology):
    probe_summary_info.sort_values(['blast_on_target_rate', 'taxon_coverage', 'quality'], ascending = [False, False, True], inplace = True)
    if (probe_summary_info['blast_on_target_rate'][0] > 0.9999 and probe_summary_info['taxon_coverage'][0] > 0.9999):
        best_probes = probe_summary_info.iloc[[0], :]
        best_probes.loc[:,'selection_method'] = 'MinOverlapSingleBest'
    elif probe_summary_info['blast_on_target_rate'][0] > 0.9999:
        taxon_molecules = [record.id for record in SeqIO.parse(taxon_fasta_filename, 'fasta')]
        taxon_molecules_set = [sub_slash(mol) for mol in set(taxon_molecules)]
        probe_summary_filtered = probe_summary_info[probe_summary_info['blast_on_target_rate'] > 0.9999]
        probe_ids = probe_summary_filtered['probe_id'].unique()
        cover_matrix = np.zeros((len(taxon_molecules), len(probe_ids)), dtype = int)
        cost = np.ones(len(probe_ids), dtype = int)
        for i in range(probe_summary_filtered.shape[0]):
            probe_idx = probe_summary_filtered.probe_id.values[i]
            probe_info = probe_summary_filtered[probe_summary_filtered.probe_id == probe_idx]
            probe_name = 'probe_' + str(probe_idx)
            probe_blast = pd.read_hdf(probe_evaluation_filename, probe_name)
            probe_blast = probe_blast[probe_blast['mch'] >= max_continuous_homology]
            blasted_molecules = list(probe_blast['molecule_id'])
            indices = [i for i, e in enumerate(blasted_molecules) if e in taxon_molecules_set]
            cover_matrix[indices,i] = 1
        cover_matrix_filt = np.array(np.delete(cover_matrix, np.where(np.sum(cover_matrix, axis = 1) == 0), axis = 0), dtype = bool)
        g = setcover.SetCover(cover_matrix_filt, cost)
        g.SolveSCP()
        set_cover_indices = np.flatnonzero(g.s*1 == 1)
        best_probes = probe_summary_filtered.iloc[set_cover_indices,:]
        if set_cover_indices.shape[0] > 1:
            best_probes.loc[:,'selection_method'] = 'MinOverlap'
        else:
            best_probes.loc[:,'selection_method'] = 'MinOverlapSingleBest'
    else:
        best_probes = probe_summary_info.iloc[[0],:]
        best_probes.loc[:,'selection_method'] = 'MinOverlapSingleBest'
    return(best_probes)

def select_top_probes(probe_summary_info, tpn, bot):
    if np.max(probe_summary_info['blast_on_target_rate']) > bot:
        best_probes_all = probe_summary_info[probe_summary_info['blast_on_target_rate'] > bot]
        if best_probes_all.shape[0] > nprobes:
            best_probes = best_probes_all.iloc[0:nprobes,:]
            best_probes.loc[:,'selection_method'] = 'Top' + str(nprobes)
        else:
            best_probes = best_probes_all
            best_probes.loc[:,'selection_method'] = 'AllTop'
    else:
        best_probes = probe_summary_info.iloc[[0],:]
        best_probes.loc[:,'selection_method'] = 'AllTopSingleBest'
    return(best_probes)

# Tm, GC calculation
def calculate_tm(df, Na = 390, dnac1_oligo = 5):
    qseq_array = df.qseq.values
    sseq_array = df.sseq.values
    tm_array = np.zeros(len(qseq_array))
    for i in range(len(qseq_array)):
        qseq = qseq_array[i]
        cseq = Seq(sseq_array[i]).complement()
        tm_array[i] = mt.Tm_NN(qseq, Na = Na, saltcorr = 7, dnac1 = dnac1_oligo*15, dnac2 = 1)
    return(tm_array)

def calculate_gc_count(df):
    qseq_array = df.qseq.values
    sseq_array = df.sseq.values
    gc_count_array = np.zeros(len(qseq_array), dtype = int)
    for i in range(len(qseq_array)):
        gc_count_array[i] = int(GC(qseq_array[i])*len(qseq_array[i])/100)
    return(gc_count_array)


def probe_blast_summarize(probe_blast, max_continuous_homology):
    probe_blast_filtered = probe_blast[(probe_blast['mch'] >= max_continuous_homology) & (probe_blast['length'] >= max_continuous_homology)]
    # print(probe_blast_filtered.shape)
    if probe_blast_filtered.shape[0] > 0:
        blast_on_target_rate = probe_blast_filtered.loc[:,'target_feature_hit'].sum()/probe_blast_filtered.shape[0]
        on_target_full_match = probe_blast.target_feature_hit_full_match.sum()/probe_blast_filtered.shape[0]
        # if blast_on_target_rate > 0.1 or on_target_full_match > 0.1:
        #     print(blast_on_target_rate, on_target_full_match)
    else:
        blast_on_target_rate = 0
        on_target_full_match = 0
    probe_blast_below_mch = probe_blast[probe_blast['mch'] < max_continuous_homology]
    if probe_blast_below_mch.shape[0] > 0:
        off_target_min_evalue = probe_blast_below_mch.evalue.min()
        off_target_max_bitscore = probe_blast_below_mch.bitscore.max()
        off_target_max_mch = probe_blast_below_mch.mch.max()
        probe_blast_below_mch.sort_values(['mch', 'qcovhsp'], ascending = [False, False], inplace = True)
        off_target_max_mch_sort = probe_blast_below_mch.mch.values[ 0]
        off_target_max_mch_qcovhsp = probe_blast_below_mch.qcovhsp.values[0]
        off_target_full_qcovhsp_fraction = np.sum(probe_blast_below_mch.qcovhsp.values > 99.9)/probe_blast_below_mch.shape[0]
    else:
        off_target_min_evalue = 100
        off_target_max_bitscore = 0
        off_target_max_mch = 0
        off_target_max_mch_sort = 0
        off_target_max_mch_qcovhsp = 0
        off_target_full_qcovhsp_fraction = 0
    return(pd.Series({'probe_id': probe_blast.probe_id.values[0],
                      'blast_on_target_rate': blast_on_target_rate,
                      'on_target_full_match': on_target_full_match,
                      'off_target_min_evalue': off_target_min_evalue,
                      'off_target_max_bitscore': off_target_max_bitscore,
                      'off_target_max_mch': off_target_max_mch,
                      'off_target_max_mch_sort': off_target_max_mch_sort,
                      'off_target_max_mch_qcovhsp': off_target_max_mch_qcovhsp,
                      'off_target_full_qcovhsp_fraction': off_target_full_qcovhsp_fraction}))

def get_probes(probe_evaluation_filename, input_probe_directory):
    target_feature = re.sub('_probe_evaluation.h5', '', os.path.basename(probe_evaluation_filename))
    probes = pd.read_table('{}/{}.int'.format(input_probe_directory, target_feature), skiprows = 3, header = None, delim_whitespace = True)
    probes.columns = ['probe_id', 'seq', 'p_start', 'ln', 'N', 'GC', 'Tm', 'self_any_th', 'self_end_th', 'hairpin', 'quality']
    probes.loc[:,'target_feature'] = target_feature
    return(target_feature, probes)


def summarize_probes(probe_evaluation_filename, probe_summary_info_filename, min_tm, max_tm, gc_cutoff, max_continuous_homology, probe_selection_method, tpn, bot, bitscore_thresh):
    sim_blast_dir, feature_evaluaton_filename = os.path.split(probe_evaluation_filename)
    sim_dir = os.path.split(sim_blast_dir)[0]
    sim_primer3_dir = '{}/primer3'.format(sim_dir)
    target_feature, probes = get_probes(probe_evaluation_filename, sim_primer3_dir)
    # print(target_feature, probes)
    feature_fasta_filename = '{}/{}.fasta'.format(sim_primer3_dir, target_feature)
    probe_summary = pd.DataFrame()
    print('There are {} probes'.format(probes.shape[0]))
    for probe_id_idx in probes.probe_id:
    # for probe_id_idx in range(1):
        probe_name = 'probe_' + str(probe_id_idx)
        # print(probe_name)
        try:
            probe_blast = pd.read_hdf(probe_evaluation_filename, probe_name)
            probe_blast['melting_temp'] = calculate_tm(probe_blast, Na, dnac1oligo)
            probe_blast['GC_count'] = calculate_gc_count(probe_blast)
            probe_blast.loc[:,'target_feature_hit'] = (probe_blast.feature.values.astype(str) == str(target_feature))
            probe_blast.loc[:,'target_feature_hit_full_match'] = (probe_blast.feature.values.astype(str) == str(target_feature))*(probe_blast.pid.values >= 99.9)*(probe_blast.qcovhsp.values >= 99.9)
            # print(probe_blast.shape)
            if probe_blast.shape[0] > 0:
                probe_summary = probe_summary.append(probe_blast_summarize(probe_blast, max_continuous_homology = max_continuous_homology), ignore_index = True, sort = True)
        except KeyError:
            pass
    if not probe_summary.empty:
        probe_summary.loc[:, 'probe_id'] = probe_summary.probe_id.astype(int)
    probes_merge = probes.merge(probe_summary, on = 'probe_id', how = 'left', sort = True)
    probes_merge = probes_merge[(probes_merge['Tm'] > min_tm) & (probes_merge['Tm'] < max_tm) & (probes_merge['GC'] > gc_cutoff)]
    probes_merge_filename = re.sub('_probe_selection.csv', '_all_probes.csv', probe_summary_info_filename)
    probes_merge.to_csv(probes_merge_filename)
    print('Selecting probes now...')
    if probes.shape[0] > 0:
        if probe_selection_method == 'SingleBestProbe':
            best_probes = pd.DataFrame(probes_merge.iloc[[0],:])
            best_probes.loc[:,'selection_method'] = 'SingleBest'
        elif probe_selection_method == 'AllSpecific':
            best_probes = select_all_specific_probes(probes_merge, bot, bitscore_thresh, probes_merge_filename)
        elif probe_selection_method == 'MinOverlap':
            best_probes = select_min_overlap_probes(probes_merge, taxon_fasta_filename, probe_evaluation_filename, max_continuous_homology)
        elif probe_selection_method == 'TopN':
            best_probes = select_top_probes(probes_merge, tpn)
    probe_blast_off_target = pd.DataFrame()
    for probe_idx in best_probes.probe_id:
        probe_info = best_probes.loc[best_probes.probe_id == probe_idx]
        probe_name = 'probe_' + str(probe_idx)
        probe_blast = pd.read_hdf(probe_evaluation_filename, probe_name)
        probe_blast_off_target = probe_blast_off_target.append(probe_blast.loc[(probe_blast.mch.values < max_continuous_homology) & (probe_blast.gapopen.values == 0)])
        probe_blast = probe_blast[probe_blast['mch'] >= max_continuous_homology]
        if not probe_blast.empty:
            probe_blast.loc[:,'target_feature_hit'] = probe_blast.feature.values.astype(str) == str(target_feature)
    if not best_probes.empty:
        best_probes.to_csv(probe_summary_info_filename, index = False)
        probe_off_target_summary_filename = re.sub('_probe_selection.csv', '_off_target_summary_info.csv', probe_summary_info_filename)
        probe_blast_off_target.to_csv(probe_off_target_summary_filename, index = False)
    else:
        best_probes = pd.DataFrame(columns = ['GC', 'N', 'Tm', 'blast_on_target_rate', 'hairpin', 'ln', 'off_target_full_qcovhsp_fraction',
                               'off_target_max_mch', 'off_target_max_mch_qcovhsp', 'off_target_max_mch_sort', 'off_target_min_evalue',
                               'on_target_full_match', 'probe_id', 'quality', 'self_any_th', 'self_end_th',
                               'seq', 'target_feature', 'selection_method'])
        best_probes.to_csv(probe_summary_info_filename, index = False)
        probe_blast_off_target = pd.DataFrame(columns = ['probe_id', 'molecule_id', 'pid', 'qcovhsp', 'length', 'mismatch', 'gapopen', 'probe_start', 'probe_end', 'molecule_start',
                                          'molecule_end', 'evalue', 'bitscore', 'staxids', 'qseq', 'sseq', 'mch'])
        probe_off_target_summary_filename = re.sub('_probe_selection.csv', '_off_target_summary_info.csv', probe_summary_info_filename)
        probe_blast_off_target.to_csv(probe_off_target_summary_filename, index = False)
    # if probes.shape[0] > 0:
    #     best_probes = pd.DataFrame()
    #     group_distance = 120
    #     while (best_probes.shape[0] <= 15) & (group_distance > 20):
    #         group_distance -= 20
    #         probes_merge.loc[:,'mean_probe_start_group'] = np.floor(probes_merge.mean_probe_start.values/group_distance).astype(int)
    #         if probes_merge.empty:
    #             print('{}, {}'.format(probe_evaluation_filename, probe_summary.off_target_max_gc.min()))
    #         best_probes = select_all_specific_p_start_group_probes(probes_merge, group_distance, bot, bitscore_thresh, mt_cutoff, ot_gc_cutoff)
    #     if best_probes.shape[0] > 15:
    #         group_distance += 20
    #         probes_merge['mean_probe_start_group'] = np.floor(probes_merge.mean_probe_start.values/group_distance).astype(int)
    #         best_probes = select_all_specific_p_start_group_probes(probes_merge, group_distance, bot, bitscore_thresh, mt_cutoff, ot_gc_cutoff)
    #     elif best_probes.shape[0] == 0:
    #         probes_sorted = probes_merge.sort_values(by = 'blast_on_target_rate', ascending = False)
    #         best_probes = pd.DataFrame(probes_sorted.iloc[[0],:])
    #         best_probes.loc[:,'selection_method'] = 'SingleBest'
    #     helper_probes_all = pd.DataFrame()
    #     for probe_idx in best_probes.probe_id:
    #         mean_probe_start = best_probes.loc[best_probes.probe_id == probe_idx, 'mean_probe_start'].values[0]
    #         mean_probe_end = best_probes.loc[best_probes.probe_id == probe_idx, 'mean_probe_end'].values[0]
    #         five_prime_helpers = probes_merge.loc[(probes_merge.mean_probe_start.values > mean_probe_start - 103) & (probes_merge.mean_probe_end.values < mean_probe_start - 3)]
    #         three_prime_helpers = probes_merge.loc[(probes_merge.mean_probe_start.values > mean_probe_end + 3) & (probes_merge.mean_probe_end.values < mean_probe_end + 103)]
    #         helper_probes = pd.concat([five_prime_helpers,three_prime_helpers])
    #         helper_probes['helper_group'] = ((helper_probes.mean_probe_start.values - mean_probe_start)/20).astype(int)
    #         helper_probes['helper_source_probe'] = probe_idx
    #         helper_probes_all = helper_probes_all.append(helper_probes)
    #     probe_off_target_summary = pd.DataFrame()
    #     probe_blast_off_target = pd.DataFrame()
    #     probe_blast_cumulative = pd.DataFrame(columns = ['molecule_id'])
    #     for probe_idx in best_probes.probe_id:
    #         probe_info = best_probes.loc[best_probes.probe_id == probe_idx]
    #         probe_blast = pd.read_hdf(probe_evaluation_filename, probe_idx)
    #         probe_blast_off_target = probe_blast_off_target.append(probe_blast.loc[(probe_blast.mch.values < max_continuous_homology) & (probe_blast.gapopen.values == 0)])
    #         probe_blast = probe_blast[probe_blast['mch'] >= max_continuous_homology]
    #         probe_blast.loc[:,'target_taxon_hit'] = probe_blast[target_rank].values.astype(str) == str(target_taxon)
    #         if probe_info.selection_method.values[0] != 'SingleBest':
    #             probe_blast_cumulative = probe_blast_cumulative.merge(probe_blast.loc[probe_blast.target_taxon_hit == True].loc[:,['molecule_id']], on = 'molecule_id', how = 'outer', sort = False)
    #         probe_off_target_summary = probe_off_target_summary.append(get_probe_blast_off_target_ranks(probe_blast, probe_info, target_rank, target_taxon), ignore_index = True, sort = False)
    #     if not probe_off_target_summary.empty:
    #         best_probes = best_probes.merge(probe_off_target_summary, on = 'probe_id', how = 'left', sort = False)
    #     if not best_probes.empty:
    #         probe_off_target_summary_filename = re.sub('_probe_selection.csv', '_off_target_summary_info.csv', probe_summary_info_filename)
    #         helper_probes_filename = re.sub('_probe_selection.csv', '_helper_probes.csv', probe_summary_info_filename)
    #         prove_coverage_filename = re.sub('_probe_selection.csv', '_probe_coverage.csv', probe_summary_info_filename)
    #         probe_blast_off_target.to_csv(probe_off_target_summary_filename, index = False)
    #         helper_probes_all.to_csv(helper_probes_filename, index = False)
    #         best_probes.to_csv(probe_summary_info_filename, index = False)
    #         probe_blast_cumulative.to_csv(prove_coverage_filename, index = False)
    #     else:
    #         best_probes = pd.DataFrame(columns = ['GC', 'N', 'Tm', 'blast_on_target_rate', 'hairpin', 'length', 'off_target_full_qcovhsp_fraction',
    #                                'off_target_max_mch', 'off_target_max_mch_qcovhsp', 'off_target_max_mch_sort', 'off_target_min_evalue',
    #                                'on_target_full_match', 'mean_probe_start', 'mean_probe_start_group', 'probe_id', 'quality', 'self_any_th', 'self_end_th',
    #                                'seq', 'target_taxon', 'target_taxon_full', 'taxon_abundance', 'taxon_coverage', 'selection_method',
    #                                'off_target_class', 'off_target_family', 'off_target_genus', 'off_target_norank', 'off_target_order',
    #                                'off_target_phylum', 'off_target_species', 'off_target_superkingdom'])
    #         best_probes.to_csv(probe_summary_info_filename, index = False)
    #         probe_blast_off_target = pd.DataFrame(columns = ['probe_id', 'molecule_id', 'pid', 'qcovhsp', 'length', 'mismatch', 'gapopen', 'probe_start', 'probe_end', 'molecule_start',
    #         	                              'molecule_end', 'evalue', 'bitscore', 'staxids', 'qseq', 'sseq', 'mch', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'])
    #         probe_off_target_summary_filename = re.sub('_probe_selection.csv', '_off_target_summary_info.csv', probe_summary_info_filename)
    #         probe_blast_off_target.to_csv(probe_off_target_summary_filename, index = False)
    #         probe_blast_cumulative = pd.DataFrame(columns = ['molecule_id'])
    #         prove_coverage_filename = re.sub('_probe_selection.csv', '_probe_coverage.csv', probe_summary_info_filename)
    #         probe_blast_cumulative.to_csv(prove_coverage_filename)
    return

###############################################################################################################
# main function
###############################################################################################################


def main():
    parser = argparse.ArgumentParser('Design FISH probes for a complex microbial community')
    # input blast filename
    parser.add_argument('probe_evaluation_complete_filename', type = str, help = 'Input file containing blast results')

    parser.add_argument('design_id', type = str, help = 'Input file containing blast lineage')

    parser.add_argument('probe_summary_info_filename', type = str, help = 'Input file containing blast lineage')

    parser.add_argument('-tmin', '--min_tm', dest = 'min_tm', type = float, default = 55.0, help = 'Boolean to indicate whether to group sequences by similarity instead of taxonomic info')

    parser.add_argument('-tmax', '--max_tm', dest = 'max_tm', type = float, default = 55.0, help = 'Boolean to indicate whether to group sequences by similarity instead of taxonomic info')

    parser.add_argument('-m', '--mch', dest = 'mch', type = int, default = 14, help = 'Boolean to indicate whether to group sequences by similarity instead of taxonomic info')

    parser.add_argument('-tpn', '--top_n_probes', dest = 'tpn', type = int, default = 1, help = 'Number of top probes to keep')

    parser.add_argument('-gc', '--gc', dest = 'gc', type = float, default = 60.0, help = 'Number of top probes to keep')

    parser.add_argument('-bot', '--bot', dest = 'bot', type = float, default = 0.0, help = 'Number of top probes to keep')

    parser.add_argument('-bt', '--bitscore_thresh', dest = 'bitscore_thresh', type = float, default = 25, help = 'Number of top probes to keep')

    parser.add_argument('-c', '--probe_selection_method', dest = 'probe_selection_method', type = str, default = 'AllSpecific', help = 'Probe selection method. AllSpecific (default) | Single Best | MinOverlap | Top N')

    args = parser.parse_args()
    sim_blast_dir, probe_evaluation_complete_filename = os.path.split(args.probe_evaluation_complete_filename)
    sim_dir = os.path.split(sim_blast_dir)[0]
    feature = re.sub('_probe_evaluation_complete.txt', '', probe_evaluation_complete_filename)
    print('Selecting probes in {} for {}'.format(feature, args.design_id))
    probe_evaluation_filename = re.sub('.complete.txt', '.h5', args.probe_evaluation_complete_filename)
    summarize_probes(probe_evaluation_filename, args.probe_summary_info_filename, args.min_tm, args.max_tm, args.gc, args.mch, args.probe_selection_method, args.tpn, args.bot, args.bitscore_thresh)
    return

if __name__ == '__main__':
    main()
