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

def select_n_specific_probes(probe_summary_info, n, bot, bitscore_thresh, mt_cutoff, ot_gc_cutoff):
    best_probes = pd.DataFrame()
    blocking_probes = pd.DataFrame()
    # Add a column for the location of the end of the probes
    probe_summary_info.loc[:,'p_end'] = probe_summary_info.p_start + probe_summary_info.ln
    # # Filter the probes to a temporary DataFrame using blast on target rate and off-target {bitscore, max melting temperature, and GC content}
    # best_probes_temp = probe_summary_info.loc[(probe_summary_info['blast_on_target_rate'] > bot) & (probe_summary_info['off_target_max_bitscore'] < bitscore_thresh) & (probe_summary_info['off_target_max_tm'] < mt_cutoff) & (probe_summary_info['off_target_max_gc'] < ot_gc_cutoff)]
    # Filter probes based on number of strong off-target hybridizations
    best_probes_temp = probe_summary_info.loc[(probe_summary_info['off_target_worst_count'] == 0)]
# # For Terminal debugging
# n = 10
# mt_cutoff = 60
# ot_gc_cutoff = 7
# bot = 0.99
# probe_summary_info = probes_merge
# probe_summary_info.columns
# best_probes_temp = probes_merge.iloc[0:10,:]
# print(best_probes_temp)
# best_probes = pd.DataFrame()
    # Check if the filtered dataframe has any probes in it
    if not best_probes_temp.shape[0] == 0:
        # If it does, sort the probes based on various off target binding characteristics
        # best_probes_temp.sort_values(['blast_on_target_rate', 'taxon_coverage', 'off_target_full_qcovhsp_fraction', 'off_target_max_mch', 'off_target_max_bitscore', 'off_target_max_tm', 'off_target_max_gc', 'on_target_full_match', 'quality'], ascending = [False, False, True, True, True, True, True, False, True], inplace = True)
        best_probes_temp.sort_values(['blast_on_target_rate', 'off_target_full_qcovhsp_fraction', 'off_target_max_mch', 'off_target_max_bitscore', 'off_target_max_tm', 'off_target_max_gc', 'on_target_full_match', 'quality'], ascending = [False, True, True, True, True, True, False, True], inplace = True)
        # Add a column for the probes with no blocking probes
        best_probes_temp.loc[:,'blocking'] = False
        # # Add a column for the location of the end of the probes
        # best_probes_temp.loc[:,'p_end'] = best_probes_temp.p_start + best_probes_temp.ln
        # Append the best probe to the output DataFrame
        best_probes = best_probes.append(best_probes_temp.iloc[[0],:], sort = True)
        # Start counter
        i = 1
        # Go through the filtered probes in the temporary DataFrame until the output DataFrame has enough probes, or until you've gone through all the probes in the DataFrame
        while (best_probes.shape[0] < n) & (best_probes_temp.shape[0] > i):
            # print(i)
            # # Terminal
            # i = 1
            # Get the start, end, and middle locations for the next best probe
            p_start = best_probes_temp.p_start.values[i]
            p_end = best_probes_temp.p_end.values[i]
            p_middle = int(np.floor((p_start + p_end)/2))
            # print(p_middle)
            # print(best_probes_temp.iloc[[0], :]['p_end'])
            # print(p_start, p_end)
            # Determine if the probe overlaps with any other probes in the output DataFrame
            bool_start = sum(((best_probes.p_start < p_start) & (best_probes.p_end > p_start)))
            bool_end = sum(((best_probes.p_start < p_end) & (best_probes.p_end > p_end)))
            bool_middle = sum(((best_probes.p_start < p_middle) & (best_probes.p_end > p_middle)))
            if (bool_start == 0) & (bool_end == 0) & (bool_middle == 0):
                # Append the non-overlapping probe to the output DataFrame
                best_probes = best_probes.append(best_probes_temp.iloc[i,:], sort = True)
            # Move to the next probe
            i += 1
    else:
        pass
    print('-\nThere are {} probes that fit the specificity conditions.\n-'.format(best_probes.shape[0]))
    # Check if the output DataFrame has enough probes
    if best_probes.shape[0] < n:
        # Check if there are probes already designed
        if not best_probes.empty:
            # Remove the chosen probes from the full probe DataFrame
            remaining_probes = probe_summary_info[~probe_summary_info['probe_id'].isin(best_probes.probe_id.astype(int))]
            # # Sort the remaining probes by off-target binding parameters
            # remaining_probes.sort_values(['blast_on_target_rate', 'off_target_full_qcovhsp_fraction', 'off_target_max_mch', 'off_target_max_bitscore', 'off_target_max_tm', 'off_target_max_gc', 'on_target_full_match', 'quality'], ascending = [False, True, True, True, True, True, False, True], inplace = True)
            # Sort the remaining probes by number of strong off-target interactions
            remaining_probes.sort_values(['off_target_worst_count'], ascending = [True], inplace = True)
            # Add a column for the probes with blocking probes
            remaining_probes.loc[:,'blocking'] = True
            # # Add a column for the location of the end of the probes
            # remaining_probes.loc[:,'p_end'] = remaining_probes.p_start + remaining_probes.ln
            # Start counter
            i = 0
        else:
            # If there were no probes from the filter, use the full probe dataframe
            remaining_probes = probe_summary_info
            # # Sort the remaining probes by off-target binding parameters
            # remaining_probes.sort_values(['blast_on_target_rate', 'off_target_full_qcovhsp_fraction', 'off_target_max_mch', 'off_target_max_bitscore', 'off_target_max_tm', 'off_target_max_gc', 'on_target_full_match', 'quality'], ascending = [False, True, True, True, True, True, False, True], inplace = True)
            # Sort the remaining probes by number of strong off-target interactions
            remaining_probes.sort_values(['off_target_worst_count'], ascending = [True], inplace = True)
            # Add a column for the probes with blocking probes
            remaining_probes.loc[:,'blocking'] = True
            # # Add a column for the location of the end of the probes
            # remaining_probes.loc[:,'p_end'] = remaining_probes.p_start + remaining_probes.ln
            # Append the best probe to the output DataFrame
            best_probes = best_probes.append(remaining_probes.iloc[0,:], sort = True)
            # Start counter
            i = 1
        # # Add a column for the probes with blocking probes
        # remaining_probes.loc[:,'blocking'] = True
        # Go through the remaining probes until the output DataFrame has enough probes, or until you've gone through all the probes in the DataFrame
        while (best_probes.shape[0] < n) & (remaining_probes.shape[0] > i):
            # print(i)
            # # Terminal
            # i = 1
            # Get the start, end, and middle locations for the next best probe
            p_start = remaining_probes.p_start.values[i]
            p_end = remaining_probes.p_end.values[i]
            p_middle = int(np.floor((p_start + p_end)/2))
            # print(p_middle)
            # print(best_probes_temp.iloc[[0], :]['p_end'])
            # print(p_start, p_end)
            # Determine if the probe overlaps with any other probes in the output DataFrame
            bool_start = sum(((best_probes.p_start < p_start) & (best_probes.p_end > p_start)))
            bool_end = sum(((best_probes.p_start < p_end) & (best_probes.p_end > p_end)))
            bool_middle = sum(((best_probes.p_start < p_middle) & (best_probes.p_end > p_middle)))
            if (bool_start == 0) & (bool_end == 0) & (bool_middle == 0):
                # Append the non-overlapping probe to the output DataFrame
                best_probes = best_probes.append(remaining_probes.iloc[i,:], sort = True)
            # Move to the next probe
            i += 1
    else:
        pass
    best_probes.loc[:,'selection_method'] = 'nSpecificProbes'
    # print('Best Probes \n', best_probes)
    # print(best_probes[['p_start','p_end']])
    # print(sum(best_probes['blocking']))

    # if best_probes.shape[0] >= n:
    #     best_probes = best_probes.iloc[[0:n-1],:]
    #     best_probes.loc[:,'selection_method'] = 'nSpecificNoBlocking'
    # else:
    #     # probe_summary_info.sort_values(['blast_on_target_rate', 'taxon_coverage', 'off_target_full_qcovhsp_fraction', 'off_target_max_mch', 'off_target_max_bitscore', 'off_target_max_tm', 'off_target_max_gc', 'on_target_full_match', 'quality'], ascending = [False, False, True, True, True, True, True, False, True], inplace = True)
    #     probe_summary_info.sort_values(['blast_on_target_rate', 'off_target_full_qcovhsp_fraction', 'off_target_max_mch', 'off_target_max_bitscore', 'off_target_max_tm', 'off_target_max_gc', 'on_target_full_match', 'quality'], ascending = [False, True, True, True, True, True, False, True], inplace = True)
    #     best_probes = best_probes.append(probe_summary_info.iloc[0,:])
    #     best_probes['p_end'] = probe_summary_info.p_start[0] + probe_summary_info.ln[0]
    #     i = 1
    #     while best_probes.shape[0] < n:
    #         p_start = probe_summary_info.p_start[i]
    #         p_end = p_start + probe_summary_info.ln[i]
    #         bool_start = sum(((best_probes.p_start < p_start) & (best_probes.p_end > p_start)))
    #         bool_end = sum(((best_probes.p_end < p_end) & (best_probes.p_end > p_end)))
    #         if (bool_start == 0) & (bool_end = 0):
    #             best_probes = best_probes.append(probe_summary_info.iloc[i,:])
    #         i = i + 1
    #     best_probes.loc[:,'selection_method'] = 'nSpecificWithBlocking'

        # probe_summary_info.sort_values(['blast_on_target_rate', 'taxon_coverage', 'off_target_full_qcovhsp_fraction', 'off_target_max_mch', 'off_target_max_bitscore', 'off_target_max_tm', 'off_target_max_gc', 'on_target_full_match', 'quality'], ascending = [False, False, True, True, True, True, True, False, True], inplace = True)
        # best_probes = probe_summary_info.iloc[[0],:]
        # best_probes.loc[:,'selection_method'] = 'AllSpecificSingleBest'
    return(best_probes)


# Select best probe for a p start group
def select_all_specific_p_start_group_probes(probe_summary_info, group_distance, bot, bitscore_thresh, mt_cutoff, ot_gc_cutoff):
    # Start a dataframe to receive the best probes
    best_probes_group = pd.DataFrame()
    # Filter the probes by blast on target rate and off-target: bitscore, max melting temperature, and GC content
    best_probes = probe_summary_info.loc[(probe_summary_info['blast_on_target_rate'] > bot) & (probe_summary_info['off_target_max_bitscore'] < bitscore_thresh) & (probe_summary_info['off_target_max_tm'] < mt_cutoff) & (probe_summary_info['off_target_max_gc'] < ot_gc_cutoff)]
    # Check if there are any probes that fit the conditions
    if not best_probes.empty:
        # Iterate through all groups
        for group in range(int(np.floor(1500/group_distance))):
            print("Group\n", group, "\n")
            # best_probes_temp = best_probes.loc[best_probes.mean_probe_start_group.values == group]
            # Find the best probes for the group
            best_probes_temp = best_probes.loc[best_probes.probe_start_group.values == group]
            print("Best probes for Group \n", best_probes_temp, "\n")
            # Sort the best probes for the group to minimize off-target hybridization measurements
            best_probes_temp.sort_values(['taxon_coverage', 'off_target_full_qcovhsp_fraction', 'off_target_max_mch', 'off_target_max_bitscore', 'off_target_max_tm', 'off_target_max_gc', 'on_target_full_match', 'quality'], ascending = [False, True, True, True, True, True, False, True], inplace = True)
            # If the group has any best probes, add the single best probe for the group to the the output DataFrame
            if not best_probes_temp.empty:
                best_probes_group = best_probes_group.append(best_probes_temp.iloc[[0],:], sort = False)
                best_probes_group.loc[:,'selection_method'] = 'AllSpecificPStartGroup'
    # If there are no probes that fit the conditions, pick the best single probe
    else:
        # probe_summary_info.sort_values(['blast_on_target_rate', 'taxon_coverage', 'off_target_full_qcovhsp_fraction', 'off_target_max_mch', 'off_target_max_bitscore', 'off_target_max_tm', 'off_target_max_gc', 'on_target_full_match', 'quality'], ascending = [False, False, True, True, True, True, True, False, True], inplace = True)
        probe_summary_info.sort_values(['blast_on_target_rate', 'off_target_full_qcovhsp_fraction', 'off_target_max_mch', 'off_target_max_bitscore', 'off_target_max_tm', 'off_target_max_gc', 'on_target_full_match', 'quality'], ascending = [False, True, True, True, True, True, False, True], inplace = True)
        best_probes_group = probe_summary_info.iloc[[0],:]
        print("Single Best probe for Group \n", best_probes_group, "\n")
        best_probes_group.loc[:,'selection_method'] = 'AllSpecificPStartGroupSingleBest'
    return(best_probes_group)

def get_off_target_last_common_taxon_rank(df, target_rank, target_taxon):
    ncbi = NCBITaxa()
    if (target_taxon != 0) & (df.loc[target_rank] != 0):
        if not pd.isnull(df.loc[target_rank]):
            last_common_taxon = ncbi.get_topology([df[target_rank], target_taxon])
            last_common_taxon_rank = last_common_taxon.rank
            if last_common_taxon_rank != 'no rank':
                lineage = ncbi.get_lineage(last_common_taxon.taxid)
                last_common_taxon_rank = ncbi.get_rank([lineage[-1]])[lineage[-1]]
            else:
                last_common_taxon_rank = 'no rank'
        else:
            last_common_taxon_rank = 'no rank'
    else:
        last_common_taxon_rank = 'no rank'
    return(last_common_taxon_rank)


# def probe_blast_summarize(probe_blast, max_continuous_homology):
# def probe_blast_summarize(probe_blast, max_continuous_homology, taxon_abundance, target_rank):
#     probe_blast_filtered = probe_blast[(probe_blast['mch'] >= max_continuous_homology) | (probe_blast['length'] >= max_continuous_homology)]
#     if probe_blast_filtered.shape[0] > 0:
#         blast_on_target_rate = probe_blast_filtered.loc[:,'target_taxon_hit'].sum()/probe_blast_filtered.shape[0]
#         taxon_coverage = probe_blast_filtered.loc[:,'target_taxon_hit'].sum()/taxon_abundance
#         on_target_full_match = probe_blast.target_taxon_hit_full_match.sum()/taxon_abundance
#         min_probe_start = probe_blast_filtered.molecule_end.min()
#         max_probe_start = probe_blast_filtered.molecule_end.max()
#         mean_probe_start = probe_blast_filtered.molecule_end.mean()
#
#         # on_target_full_match = probe_blast.target_feature_hit_full_match.sum()/probe_blast_filtered.shape[0]
#
#     else:
#         blast_on_target_rate = 0
#         taxon_coverage = 0
#         on_target_full_match = 0
#         min_probe_start = 0
#         max_probe_start = 0
#         mean_probe_start = 0
#     probe_blast_below_mch = probe_blast[probe_blast['mch'] < max_continuous_homology]
#     if probe_blast_below_mch.shape[0] > 0:
#         off_target_min_evalue = probe_blast_below_mch.evalue.min()
#         off_target_max_bitscore = probe_blast_below_mch.bitscore.max()
#         off_target_max_mch = probe_blast_below_mch.mch.max()
#         probe_blast_below_mch.sort_values(['mch', 'qcovhsp'], ascending = [False, False], inplace = True)
#         off_target_max_mch_sort = probe_blast_below_mch.mch.values[0]
#         off_target_max_mch_qcovhsp = probe_blast_below_mch.qcovhsp.values[0]
#         off_target_full_qcovhsp_fraction = np.sum(probe_blast_below_mch.qcovhsp.values > 99.9)/probe_blast_below_mch.shape[0]
#         off_target_max_tm = probe_blast_below_mch.melting_temp.max()
#         off_target_max_gc = probe_blast_below_mch.GC_count.max()
#     else:
#         off_target_min_evalue = 100
#         off_target_max_bitscore = 0
#         off_target_max_mch = 0
#         off_target_max_mch_sort = 0
#         off_target_max_mch_qcovhsp = 0
#         off_target_full_qcovhsp_fraction = 0
#         off_target_max_tm = 0
#         off_target_max_gc = 0
#     return(pd.Series({'probe_id': probe_blast.probe_id.values[0],
#                       'min_probe_start': min_probe_start,
#                       'max_probe_start': max_probe_start,
#                       'mean_probe_start': mean_probe_start,
#                       'blast_on_target_rate': blast_on_target_rate,
#                       'taxon_coverage': taxon_coverage,
#                       'on_target_full_match': on_target_full_match,
#                       'off_target_min_evalue': off_target_min_evalue,
#                       'off_target_max_bitscore': off_target_max_bitscore,
#                       'off_target_max_mch': off_target_max_mch,
#                       'off_target_max_mch_sort': off_target_max_mch_sort,
#                       'off_target_max_mch_qcovhsp': off_target_max_mch_qcovhsp,
#                       'off_target_full_qcovhsp_fraction': off_target_full_qcovhsp_fraction,
#                       'off_target_max_tm': off_target_max_tm,
#                       'off_target_max_gc': off_target_max_gc,
#                       'taxon_abundance': taxon_abundance}))
# Append the analysis of blast results to the summary of all probes
# def probe_blast_summarize(probe_blast, target_feature, max_continuous_homology, bot, bitscore_thresh, ot_gc_cutoff):
def probe_blast_summarize(probe_blast, max_continuous_homology, probe_blast_off_target_worst):
    # # Sort the off target probes by off-target binding parameters
    # probe_blast_off_target.sort_values(['mch', 'bitscore', 'melting_temp', 'GC_count'], ascending = [False, False, False, False], inplace = True)
    # Filter the blast results based on those with significant binding potential
    probe_blast_filtered = probe_blast[(probe_blast['mch'] >= max_continuous_homology) & (probe_blast['length'] >= max_continuous_homology)]
    # # Create a list of the off-target bindings
    # probe_blast_off_target = probe_blast_filtered[probe_blast_filtered.feature.values.astype(str) != str(target_feature)]
    # # Filter the probes to a temporary DataFrame using blast on target rate and off-target {bitscore, max melting temperature, and GC content}
    # probe_blast_off_target_worst = probe_blast_off_target.loc[(probe_blast_off_target['blast_on_target_rate'] > bot) & (probe_blast_off_target['off_target_max_bitscore'] < bitscore_thresh) & (probe_blast_off_target['off_target_max_tm'] < mt_cutoff) & (probe_blast_off_target['off_target_max_gc'] < ot_gc_cutoff)]
    # print(probe_blast_filtered.shape)
    # if the filter did not remove all probes, calculate the rate at which probes hybridize on target
    if probe_blast_filtered.shape[0] > 0:
        blast_on_target_rate = probe_blast_filtered.loc[:,'target_feature_hit'].sum()/probe_blast_filtered.shape[0]
        on_target_full_match = probe_blast.target_feature_hit_full_match.sum()/probe_blast_filtered.shape[0]
        # Determine the most common blast site
        probe_start = probe_blast_filtered.molecule_end.mode()[0]
        # Get the count of off-target probes that have significant homology
        off_target_worst_count = probe_blast_off_target_worst.shape[0]
        # if blast_on_target_rate > 0.1 or on_target_full_match > 0.1:
        #     print(blast_on_target_rate, on_target_full_match)
    else:
        blast_on_target_rate = 0
        on_target_full_match = 0
        probe_start = probe_blast.molecule_end.mode()[0]
    # Filter the blast results based on those with short homology regions
    probe_blast_below_mch = probe_blast[probe_blast['mch'] < max_continuous_homology]
    # If the filter did not remove all probes, calculate the minimums and maximums for various readouts of off-target binding
    if probe_blast_below_mch.shape[0] > 0:
        off_target_min_evalue = probe_blast_below_mch.evalue.min()
        off_target_max_bitscore = probe_blast_below_mch.bitscore.max()
        off_target_max_mch = probe_blast_below_mch.mch.max()
        probe_blast_below_mch.sort_values(['mch', 'qcovhsp'], ascending = [False, False], inplace = True)
        off_target_max_mch_sort = probe_blast_below_mch.mch.values[ 0]
        off_target_max_mch_qcovhsp = probe_blast_below_mch.qcovhsp.values[0]
        off_target_full_qcovhsp_fraction = np.sum(probe_blast_below_mch.qcovhsp.values > 99.9)/probe_blast_below_mch.shape[0]
        off_target_max_tm = probe_blast_below_mch.melting_temp.max()
        off_target_max_gc = probe_blast_below_mch.GC_count.max()
        # off_target_worst_count = probe_blast_off_target_worst.shape[0]
    else:
        off_target_min_evalue = 100
        off_target_max_bitscore = 0
        off_target_max_mch = 0
        off_target_max_mch_sort = 0
        off_target_max_mch_qcovhsp = 0
        off_target_full_qcovhsp_fraction = 0
        off_target_max_tm = 0
        off_target_max_gc = 0
        # off_target_worst_count = probe_blast_off_target_worst.shape[0]
    # Return the results of the various analyses as a dictionary
    return(pd.Series({'probe_id': probe_blast.probe_id.values[0],
                      'blast_on_target_rate': blast_on_target_rate,
                      'on_target_full_match': on_target_full_match,
                      'probe_start': probe_start,
                      'off_target_min_evalue': off_target_min_evalue,
                      'off_target_max_bitscore': off_target_max_bitscore,
                      'off_target_max_mch': off_target_max_mch,
                      'off_target_max_mch_sort': off_target_max_mch_sort,
                      'off_target_max_mch_qcovhsp': off_target_max_mch_qcovhsp,
                      'off_target_full_qcovhsp_fraction': off_target_full_qcovhsp_fraction,
                      'off_target_max_tm': off_target_max_tm,
                      'off_target_max_gc': off_target_max_gc,
                      'off_target_worst_count': off_target_worst_count
                      }))

def get_probe_blast_off_target_ranks(probe_blast, probe_info, target_rank, target_taxon):
    ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'no rank']
    if probe_info.blast_on_target_rate.values[0] < 1-1e-8:
        probe_blast_filtered = probe_blast[probe_blast['target_taxon_hit'] == False]
        off_target_last_common_taxon_rank = probe_blast_filtered.apply(get_off_target_last_common_taxon_rank, target_rank = target_rank, target_taxon = target_taxon, axis = 1)
        n_total = off_target_last_common_taxon_rank.shape[0]
        if n_total > 0:
            rank_fractions = [np.sum(off_target_last_common_taxon_rank == rank)/n_total for rank in ranks]
        else:
            rank_fractions = [np.nan for rank in ranks]
    else:
        rank_fractions = [0 for rank in ranks]
    probe_off_target_summary = pd.DataFrame([probe_blast.probe_id.values[0]] + rank_fractions).transpose()
    probe_off_target_summary.columns = ['probe_id'] + ['off_target_' + s for s in ranks]
    return(probe_off_target_summary)

def sub_slash(str):
    return(re.sub('/', '_', str))

# calculate the melting temperature of off target matches from blast results
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

def get_blast_lineage(blast_lineage_filename):
    blast_lineage_df = pd.read_table(blast_lineage_filename, dtype = str)
    lineage_columns = ['molecule_id', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    blast_lineage_slim = blast_lineage_df[lineage_columns]
    blast_lineage_slim.loc[:,'molecule_id'] = blast_lineage_slim.molecule_id.apply(sub_slash)
    return(blast_lineage_slim)

def get_taxon_abundance(blast_lineage, target_taxon):
    taxon_abundance = blast_lineage
    return(taxon_abundance)

# def get_probes(probe_evaluation_filename, probes_summary_filename, target_rank):
#     target_taxon = re.sub('_probe_evaluation.h5', '', os.path.basename(probe_evaluation_filename))
#     probes = pd.read_hdf(probes_summary_filename, key = '{}/{}_{}'.format(target_rank, target_rank, target_taxon), mode = 'r')
#     probes.loc[:,'target_taxon'] = target_taxon
#     return(target_taxon, probes)
# Get probes and feature name from evaluation file
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

def get_full_probe_id(pid, target_rank, target_taxon):
    return('{}_{}_{}'.format(target_rank, target_taxon, pid))

# def summarize_probes(probe_evaluation_filename, blast_lineage_filename, probe_summary_info_filename, sample_directory, probes_summary_filename, target_rank, min_tm, max_tm, gc_cutoff, max_continuous_homology, bot, bitscore_thresh, Na, dnac1oligo, mt_cutoff, ot_gc_cutoff):
def summarize_probes(probe_evaluation_filename, probe_summary_info_filename, sample_directory, probes_summary_filename, probe_selection_level, total_transcript_table_filename, n, min_tm, max_tm, gc_cutoff, max_continuous_homology, bot, bitscore_thresh, Na, dnac1oligo, mt_cutoff, ot_gc_cutoff):
    # # evaluation_dir, taxon_evaluaton_filename = os.path.split(probe_evaluation_filename)
    # # blast_lineage = get_blast_lineage(blast_lineage_filename)
    # # taxon_abundance = blast_lineage.groupby(target_rank).molecule_id.count().reset_index()
    # # taxon_abundance.loc[:, target_rank] = taxon_abundance.loc[:, target_rank].astype(str)
    # # target_taxon, probes = get_probes(probe_evaluation_filename, probes_summary_filename, target_rank)
    # # target_taxon_abundance = taxon_abundance.loc[taxon_abundance.loc[:,target_rank].values == target_taxon, 'molecule_id'].values[0]
    # # probe_summary = pd.DataFrame()
# # # Terminal editing
# probe_evaluation_filename = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/runs/2020_02_13_gfp_plasmid/sample_002/blast/GFP_probe_evaluation.h5'
# sim_primer3_dir = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/runs/2020_02_13_gfp_plasmid/sample_002/primer3'
    # Get info from eval filename for future filenames
    sim_blast_dir, feature_evaluaton_filename = os.path.split(probe_evaluation_filename)
    sim_dir = os.path.split(sim_blast_dir)[0]
    sim_primer3_dir = '{}/primer3'.format(sim_dir)
    # Get probes and feature name from evaluation file
    target_feature, probes = get_probes(probe_evaluation_filename, sim_primer3_dir)
    # print(target_feature, probes)
    # build a fasta for the feature
    feature_fasta_filename = '{}/{}.fasta'.format(sim_primer3_dir, target_feature)
    # initiate a dataframe to summarize all probes
    probe_summary = pd.DataFrame()
    # iterate through the probe ids
    for probe_idx in probes.probe_id:
        # probe_name = '{}_{}/{}_{}_{}'.format(target_rank, target_taxon, target_rank, target_taxon, probe_idx)
        # generate a new probe name
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
        if probe_selection_level == 'genome_bin':
            total_transcript_table = pd.read_csv(total_transcript_table_filename)
            feature_transcript_info = total_transcript_table[total_transcript_table['Gene'] == str(target_feature)]
            phylogenetic_marker_genes = ['dnaG', 'frr', 'infC', 'nusA', 'pgk', 'pyrG', 'rplA', 'rplB', 'rplC', 'rplD', 'rplE', 'rplF', 'rplK', 'rplL', 'rplM', 'rplN',
            'rplP', 'rplS', 'rplT', 'rpmA', 'rpoB', 'rpsB', 'rpsC', 'rpsE', 'rpsI', 'rpsJ', 'rpsK', 'rpsM', 'rpsS', 'smpB', 'tsf']
            feature_genome_bin = total_transcript_table[total_transcript_table['Genome_bin'] == feature_transcript_info['Genome_bin']]
            feature_genome_bin_marker = feature_genome_bin[feature_genome_bin['Preferred_name'].isin(phylogenetic_marker_genes)]
            genome_bin_feature_list = feature_genome_bin_marker['Gene'].values
            probe_blast_off_target =  probe_blast[~probe_blast.feature.isin(genome_bin_feature_list)]
        elif probe_selection_level == 'molecule':
                probe_blast_off_target = probe_blast[probe_blast.feature.values.astype(str) != str(target_feature)]
        # Filter the probes to a temporary DataFrame using blast on target rate and off-target {bitscore, max melting temperature, and GC content}
        probe_blast_off_target_worst = probe_blast_off_target.loc[(probe_blast_off_target['mch'] > max_continuous_homology) |
                                                                (probe_blast_off_target['bitscore'] > bitscore_thresh) |
                                                                (probe_blast_off_target['melting_temp'] > mt_cutoff) |
                                                                (probe_blast_off_target['GC_count'] > ot_gc_cutoff)]
        # Append the probe blast analysis to the summary of all probes
        probe_summary = probe_summary.append(probe_blast_summarize(probe_blast, max_continuous_homology, probe_blast_off_target_worst), ignore_index = True, sort = True)
    # probe_summary = probe_summary.append(probe_blast_summarize(probe_blast, target_feature, max_continuous_homology, bot, bitscore_thresh, ot_gc_cutoff), ignore_index = True, sort = True)

    # probes.loc[:,'probe_id'] = probes.probe_id.apply(get_full_probe_id, args = (target_rank, target_taxon))
    # Generate a column with the probe names as integers
    probe_summary.loc[:, 'probe_id'] = probe_summary.probe_id.astype(int)
    # Merge the summary of blast results with the probe design file
    probes_merge = probes.merge(probe_summary, on = 'probe_id', how = 'left', sort = False)
    # Filter the probes based on those that satisfy the melting temperature and GC content requirements
    probes_merge = probes_merge.loc[(probes_merge['Tm'] >= min_tm) & (probes_merge['Tm'] <= max_tm) & (probes_merge['GC'] >= gc_cutoff),:]
    # probes_merge['mean_probe_end'] = probes_merge.mean_probe_start.values + probes_merge.length.values - 1
    # probes_merge['mean_probe_start_group'] = np.floor(probes_merge.mean_probe_start.values/20).astype(int)
    # # Create groups for probes
    # probes_merge['probe_start_group'] = np.floor(probes_merge.molecule_end.values/20).astype(int)
# # # Input for terminal editing
# probe_summary_info_filename = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/runs/2020_02_13_gfp_plasmid/simulation/DSGN_002/GFP_probe_selection.csv'
    # Set up a csv for the probe blast analysis results and write the merged probe info to the file
    # return the total number of probes
    # print('There are {} probes'.format(probes_merge.shape[0]))
    probes_merge_filename = re.sub('_probe_selection.csv', '_all_probes.csv', probe_summary_info_filename)
    probes_merge.to_csv(probes_merge_filename)

# # Editing: Read in written csv for code troubleshooting
# probes_merge_filename = re.sub('_probe_selection.csv', '_all_probes.csv', probe_summary_info_filename)
# probes_merge = pd.read_csv(probes_merge_filename)
# probes_merge[probes_merge['probe_id'] == 554]
    # probes_merge.shape
    # probes_merge.columns
    # probes_merge.probe_start[0:10]
    # probes_merge.p_start[0:10]

    # Check if any probes exist for this feature
    if probes.shape[0] > 0:
        # Start a DataFrame for the selected probes
        best_probes = pd.DataFrame()
        # # Select the number of probes to design
        # n = 2
# # # Terminal Editing
# # mt_cutoff = 60
# # ot_gc_cutoff = 7
# # bot = 0.99

        # Select all specific probes
        best_probes = select_n_specific_probes(probes_merge, n, bot, bitscore_thresh, mt_cutoff, ot_gc_cutoff)
        best_probes.to_csv(probe_summary_info_filename, index = False)

        # # Editing
        # best_probes = pd.read_csv(probe_summary_info_filename)

        # Set a spacing of probes along the feature
        # group_distance = 120
        # Select p start group probes for the first 15 best probes
        # while (best_probes.shape[0] <= 15) & (group_distance > 20):
        # # while (best_probes.shape[0] <= 2) & (group_distance > 20):
        #     # Subtract 20 from the group distance for every best probe added
        #     group_distance -= 20
        #     print("Group Distance\n", group_distance, "\n")
        #     # probes_merge.loc[:,'mean_probe_start_group'] = np.floor(probes_merge.mean_probe_start.values/group_distance).astype(int)
        #     probes_merge.loc[:,'probe_start_group'] = np.floor(probes_merge.probe_start.values/group_distance).astype(int)
        #     print("Probe Start Group\n", probes_merge.probe_start_group[0:10], "\n")
    #         if probes_merge.empty:
    #             print('{}, {}'.format(probe_evaluation_filename, probe_summary.off_target_max_gc.min()))
            # Select the best probes for a start group
            # best_probes = select_all_specific_p_start_group_probes(probes_merge, group_distance, bot, bitscore_thresh, mt_cutoff, ot_gc_cutoff)
    #     if best_probes.shape[0] > 15:
    #         group_distance += 20
    #         probes_merge['mean_probe_start_group'] = np.floor(probes_merge.mean_probe_start.values/group_distance).astype(int)
    #         best_probes = select_all_specific_p_start_group_probes(probes_merge, group_distance, bot, bitscore_thresh, mt_cutoff, ot_gc_cutoff)
    #     elif best_probes.shape[0] == 0:
    #         probes_sorted = probes_merge.sort_values(by = 'blast_on_target_rate', ascending = False)
    #         best_probes = pd.DataFrame(probes_sorted.iloc[[0],:])
    #         best_probes.loc[:,'selection_method'] = 'SingleBest'
        # Set up a DataFrame to recieve helper probes
        print('Selected {} probes'.format(best_probes.shape[0]))

        helper_probes = pd.DataFrame(columns = probes_merge.columns)
        helper_probes['p_end'] = []
        # Iterate through best probes
        for probe_idx in best_probes.probe_id:
            # mean_probe_start = best_probes.loc[best_probes.probe_id == probe_idx, 'mean_probe_start'].values[0]
            # mean_probe_end = best_probes.loc[best_probes.probe_id == probe_idx, 'mean_probe_end'].values[0]
            # five_prime_helpers = probes_merge.loc[(probes_merge.mean_probe_start.values > mean_probe_start - 103) & (probes_merge.mean_probe_end.values < mean_probe_start - 3)]
            # three_prime_helpers = probes_merge.loc[(probes_merge.mean_probe_start.values > mean_probe_end + 3) & (probes_merge.mean_probe_end.values < mean_probe_end + 103)]
            # Get the start and end of the best probe
            p_start = best_probes.loc[best_probes.probe_id == probe_idx, 'p_start'].values[0]
            p_end = best_probes.loc[best_probes.probe_id == probe_idx, 'p_end'].values[0]
            # Add a column to all summarized probes for the location of the end of the probes
            probes.loc[:,'p_end'] = probes.p_start + probes.ln
            # probes_merge.loc[:,'p_end'] = probes_merge.p_start + probes_merge.ln
            # Design helper probes that go on either end of the best probe
            # five_prime_helpers = probes_merge.loc[(probes_merge.p_start.values > p_start - 103) & (probes_merge.p_end.values < p_start - 3), :]
            # three_prime_helpers = probes_merge.loc[(probes_merge.p_start.values > p_end + 3) & (probes_merge.p_end.values < p_end + 103), :]
            five_prime_helpers = probes.loc[(probes.p_start.values > p_start - 33) & (probes.p_end.values < p_start - 3), :]
            three_prime_helpers = probes.loc[(probes.p_start.values > p_end + 3) & (probes.p_end.values < p_end + 33), :]
            # print(p_start, p_end, probes_merge.shape,five_prime_helpers.shape,three_prime_helpers.shape)
            # Concatenate the helpers on each end together
            helper_probes_temp = pd.concat([five_prime_helpers,three_prime_helpers])
            # print('There are {} potential helper probes'.format(helper_probes_temp.shape[0]))
            for helper in helper_probes_temp.probe_id:
                # Get the start and end positions for the helper probe
                h_start = helper_probes_temp.loc[helper_probes_temp.probe_id == helper, 'p_start'].values[0]
                h_end = helper_probes_temp.loc[helper_probes_temp.probe_id == helper, 'p_end'].values[0]
                h_middle = int(np.floor((h_start + h_end)/2))
                # print(h_start, h_middle, h_end)
                # Determine if the helper overlaps with any other probes
                bool_start_p = sum((best_probes.p_start < h_start) & (best_probes.p_end > h_start))
                bool_start_h = sum((helper_probes.p_start < h_start) & (helper_probes.p_end > h_start))
                bool_end_p = sum((best_probes.p_start < h_end) & (best_probes.p_end > h_end))
                bool_end_h = sum((helper_probes.p_start < h_end) & (helper_probes.p_end > h_end))
                bool_middle_p = sum((best_probes.p_start < h_middle) & (best_probes.p_end > h_middle))
                bool_middle_h = sum((helper_probes.p_start < h_middle) & (helper_probes.p_end > h_middle))
                # print(bool_start_p, bool_middle_p, bool_end_p,bool_start_h, bool_middle_h, bool_end_h)
                if (bool_start_p == 0) & (bool_end_p == 0) & (bool_middle_p == 0) & (bool_start_h == 0) & (bool_end_h == 0) & (bool_middle_h == 0):
                    # Append the non-overlapping helper to the output DataFrame
                    helper_probes = helper_probes.append(helper_probes_temp.loc[helper_probes_temp.probe_id == helper, :], sort = True)
            # helper_probes['helper_group'] = ((helper_probes.mean_probe_start.values - mean_probe_start)/20).astype(int)
            # helper_probes['helper_group'] = ((helper_probes.mean_probe_start.values - mean_probe_start)/20).astype(int)
            # Add a column for the best probe id
            helper_probes['helper_source_probe'] = probe_idx
            # print('There are {} helper probes'.format(helper_probes.shape[0]))
            # # Add the helpers to the output DataFrame
            # helper_probes_all = helper_probes.append(helper_probes)
# best_probes[['probe_id','p_start','p_end']]
# helper_probes[['p_start','p_end']]
# # Terminal Editing
# probe_summary_info_filename = '/workdir/bmg224/hiprfish/transcript_label/probe_design/runs/2020_02_05_hao_zhou_mouse_gut/simulation/DSGN_002/NODE_42025_length_3778_cov_100.983615_3_probe_selection.csv'
# best_probes = pd.read_csv(probe_summary_info_filename)
# best_probes[['probe_id','off_target_max_mch', 'off_target_max_tm']]
# min(probes_merge.off_target_max_tm.values)
# np.mean(probes_merge.off_target_max_tm.values)
# min(probes_merge.off_target_max_gc.values)
# probes_merge.off_target_max_tm[probes_merge.off_target_max_tm.values <= 60]
# probes_merge.off_target_max_gc[probes_merge.off_target_max_gc.values <= 10]
# Na = 390
# dnac1oligo = 5
# max_continuous_homology = 14
# bitscore_thresh = 29
        # Initiate DataFrames to output the off target binding of selected probes
        probe_off_target_summary = pd.DataFrame()
        probe_blast_off_target_all = pd.DataFrame()
        probe_blast_cumulative = pd.DataFrame(columns = ['molecule_id'])
        blocking_probes = pd.DataFrame()
        # Iterate through all selected probes
        for probe_idx in best_probes.probe_id:
            # ensure id is an integer
            probe_idx = int(probe_idx)
            # Get probe names
            probe_info = best_probes.loc[best_probes.probe_id == probe_idx]
            print(probe_info[['seq','off_target_worst_count']])
            # print('-\nProbe {}\n-\n{}'.format(probe_idx, probe_info[['seq','off_target_worst_count']]))
            # Get the hdf probe name
            probe_name = 'probe_' + str(probe_idx)
            # Read in the probe blast results
            probe_blast = pd.read_hdf(probe_evaluation_filename, probe_name)
            # # write the melting temp and gc content of each off target match
            probe_blast['melting_temp'] = calculate_tm(probe_blast, Na, dnac1oligo)
            probe_blast['GC_count'] = calculate_gc_count(probe_blast)
            if probe_selection_level == 'genome_bin':
                total_transcript_table = pd.read_csv(total_transcript_table_filename)
                feature_transcript_info = total_transcript_table[total_transcript_table['Gene'] == str(target_feature)]
                phylogenetic_marker_genes = ['dnaG', 'frr', 'infC', 'nusA', 'pgk', 'pyrG', 'rplA', 'rplB', 'rplC', 'rplD', 'rplE', 'rplF', 'rplK', 'rplL', 'rplM', 'rplN',
                'rplP', 'rplS', 'rplT', 'rpmA', 'rpoB', 'rpsB', 'rpsC', 'rpsE', 'rpsI', 'rpsJ', 'rpsK', 'rpsM', 'rpsS', 'smpB', 'tsf']
                feature_genome_bin = total_transcript_table[total_transcript_table['Genome_bin'] == feature_transcript_info['Genome_bin']]
                feature_genome_bin_marker = feature_genome_bin[feature_genome_bin['Preferred_name'].isin(phylogenetic_marker_genes)]
                genome_bin_feature_list = feature_genome_bin_marker['Gene'].values
                probe_blast_off_target =  probe_blast[~probe_blast.feature.isin(genome_bin_feature_list)]
            elif probe_selection_level == 'molecule':
                probe_blast_off_target = probe_blast[probe_blast.feature.values.astype(str) != str(target_feature)]
            probe_blast_off_target.sort_values(['qcovhsp', 'mch', 'bitscore', 'melting_temp', 'GC_count'], ascending = [False, False, False, False, False], inplace = True)
            probe_blast_off_target_all = probe_blast_off_target_all.append(probe_blast_off_target)
            # print(probe_info)
            # # Read in the probe blast results
            # probe_blast = pd.read_hdf(probe_evaluation_filename, probe_idx)
            # # write the melting temp and gc content of each off target match
            # probe_blast['melting_temp'] = calculate_tm(probe_blast, Na, dnac1oligo)
            # probe_blast['GC_count'] = calculate_gc_count(probe_blast)
            # # Append blast results for off target sites where there are no gaps in the binding region
            # probe_blast_off_target = probe_blast.loc[(probe_blast.mch.values < max_continuous_homology) & (probe_blast.gapopen.values == 0)]
            # # Reduce the blast results to those with long sequences
            # probe_blast_long = probe_blast[probe_blast['mch'] >= max_continuous_homology]
            # # Determine which probes blast only to their target with long homology
            # probe_blast.loc[:,'target_feature_hit'] = (probe_blast.feature.values.astype(str) == str(target_feature))
            # probe_blast.loc[:,'target_taxon_hit'] = probe_blast[target_rank].values.astype(str) == str(target_taxon)
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
                # print('There are {} blocking probes for probe {}'.format(probe_blast_off_target_worst.shape[0],probe_idx))
                blocking_probes = blocking_probes.append(probe_blast_off_target_worst)
                # probe_blast_off_target_all = probe_blast_off_target_all.append(probe_blast_off_target)

        #
        #     # Append blast results for off target sites where there are no gaps in the binding region
        #     probe_blast_off_target_all = probe_blast_off_target_all.append(probe_blast.loc[(probe_blast.mch.values < max_continuous_homology) & (probe_blast.gapopen.values == 0)])

            # if probe_info.selection_method.values[0] != 'SingleBest':
            #     probe_blast_cumulative = probe_blast_cumulative.merge(probe_blast.loc[probe_blast.target_feature_hit == True].loc[:,['molecule_id']], on = 'molecule_id', how = 'outer', sort = False)
                # probe_blast_cumulative = probe_blast_cumulative.merge(probe_blast.loc[probe_blast.target_taxon_hit == True].loc[:,['molecule_id']], on = 'molecule_id', how = 'outer', sort = False)
    #         probe_off_target_summary = probe_off_target_summary.append(get_probe_blast_off_target_ranks(probe_blast, probe_info, target_rank, target_taxon), ignore_index = True, sort = False)
        # if not probe_off_target_summary.empty:
        #     best_probes = best_probes.merge(probe_off_target_summary, on = 'probe_id', how = 'left', sort = False)
        if not best_probes.empty:
            probe_off_target_summary_filename = re.sub('_probe_selection.csv', '_off_target_summary_info.csv', probe_summary_info_filename)
            print('There are {} blocking probes.'.format(blocking_probes.shape[0]))
            blocking_probes_filename = re.sub('_probe_selection.csv', '_blocking_probes.csv', probe_summary_info_filename)
            print('There are {} helper probes'.format(helper_probes.shape[0]))
            helper_probes_filename = re.sub('_probe_selection.csv', '_helper_probes.csv', probe_summary_info_filename)
            # prove_coverage_filename = re.sub('_probe_selection.csv', '_probe_coverage.csv', probe_summary_info_filename)
            probe_blast_off_target_all.to_csv(probe_off_target_summary_filename, index = False)
            blocking_probes.to_csv(blocking_probes_filename, index = False)
            helper_probes.to_csv(helper_probes_filename, index = False)
            best_probes.to_csv(probe_summary_info_filename, index = False)
            # probe_blast_cumulative.to_csv(prove_coverage_filename, index = False)
        # else:
        #     pass
        else:
            best_probes = pd.DataFrame(columns = ['GC', 'N', 'Tm', 'blast_on_target_rate', 'hairpin', 'length', 'off_target_full_qcovhsp_fraction',
                                   'off_target_max_mch', 'off_target_max_mch_qcovhsp', 'off_target_max_mch_sort', 'off_target_min_evalue',
                                   'on_target_full_match', 'mean_probe_start', 'mean_probe_start_group', 'probe_id', 'quality', 'self_any_th', 'self_end_th',
                                   'seq', 'target_taxon', 'target_taxon_full', 'taxon_abundance', 'taxon_coverage', 'selection_method',
                                   'off_target_class', 'off_target_family', 'off_target_genus', 'off_target_norank', 'off_target_order',
                                   'off_target_phylum', 'off_target_species', 'off_target_superkingdom'])
            best_probes.to_csv(probe_summary_info_filename, index = False)
            probe_blast_off_target_all = pd.DataFrame(columns = ['probe_id', 'molecule_id', 'pid', 'qcovhsp', 'length', 'mismatch', 'gapopen', 'probe_start', 'probe_end', 'molecule_start',
            	                              'molecule_end', 'evalue', 'bitscore', 'staxids', 'qseq', 'sseq', 'mch', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'])
            probe_off_target_summary_filename = re.sub('_probe_selection.csv', '_off_target_summary_info.csv', probe_summary_info_filename)
            probe_blast_off_target_all.to_csv(probe_off_target_summary_filename, index = False)
            probe_blast_cumulative = pd.DataFrame(columns = ['molecule_id'])
            prove_coverage_filename = re.sub('_probe_selection.csv', '_probe_coverage.csv', probe_summary_info_filename)
            probe_blast_cumulative.to_csv(prove_coverage_filename)
    return


###############################################################################################################
# main function
###############################################################################################################

def main():
    parser = argparse.ArgumentParser('Design FISH probes for a complex microbial community')

    parser.add_argument('probe_evaluation_complete_filename', type = str, help = 'Input file containing blast results')

    parser.add_argument('design_id', type = str, help = 'Input file containing blast lineage')

    parser.add_argument('probe_summary_info_filename', type = str, help = 'Input file containing blast lineage')

    # parser.add_argument('-t', '--target_rank', dest = 'target_rank', type = str, default = 'phylum', help = 'Input file containing blast lineage')

    parser.add_argument('-n', '--number_of_probes', dest = 'number_of_probes', type = int, default = 2, help = 'Indicate the number of probes to design for the target')

    parser.add_argument('-tmin', '--min_tm', dest = 'min_tm', type = float, default = 55.0, help = 'Boolean to indicate whether to group sequences by similarity instead of taxonomic info')

    parser.add_argument('-tmax', '--max_tm', dest = 'max_tm', type = float, default = 55.0, help = 'Boolean to indicate whether to group sequences by similarity instead of taxonomic info')

    parser.add_argument('-m', '--mch', dest = 'mch', type = int, default = 14, help = 'Boolean to indicate whether to group sequences by similarity instead of taxonomic info')

    parser.add_argument('-gc', '--gc', dest = 'gc', type = float, default = 40.0, help = 'Number of top probes to keep')

    parser.add_argument('-bot', '--bot', dest = 'bot', type = float, default = 0.0, help = 'Number of top probes to keep')

    parser.add_argument('-bt', '--bitscore_thresh', dest = 'bitscore_thresh', type = float, default = 27, help = 'Number of top probes to keep')

    parser.add_argument('-sod', '--sod', dest = 'sod', type = float, default = 390, help = 'sodium concentration in nM')

    parser.add_argument('-dnaconc', '--dnaconc', dest = 'dnaconc', type = float, default = 5, help = 'oligo concentration in nM')

    parser.add_argument('-mt', '--mt_cutoff', dest = 'mt_cutoff', type = float, default = 60, help = 'oligo concentration in nM')

    parser.add_argument('-otgc', '--ot_gc_cutoff', dest = 'ot_gc_cutoff', type = float, default = 7, help = 'oligo concentration in nM')

    parser.add_argument('-psl', '--probe_selection_level', dest = 'probe_selection_level', type = str, default = 'molecule', help = 'how general is the probe selection')

    parser.add_argument('-ttt', '--total_transcript_table_filename', dest = 'total_transcript_table_filename', type = str, default = 'molecule', help = 'how general is the probe selection')

    args = parser.parse_args()
    design_level_directory = os.path.split(args.probe_evaluation_complete_filename)[0]
    evaluation_dir = os.path.split(design_level_directory)[0]
    sample_dir = os.path.split(evaluation_dir)[0]
    probe_evaluation_filename = re.sub('_complete.txt', '.h5', args.probe_evaluation_complete_filename)
    probes_summary_filename = '{}/probes_summary/probes_summary.h5'.format(sample_dir)
    blast_lineage_filename = '{}/utilities/blast_lineage.tab'.format(sample_dir)
    # summarize_probes(probe_evaluation_filename, blast_lineage_filename, args.probe_summary_info_filename, sample_dir, probes_summary_filename, args.target_rank, args.min_tm, args.max_tm, args.gc, args.mch, args.bot, args.bitscore_thresh, args.sod, args.dnaconc, args.mt_cutoff, args.ot_gc_cutoff)
    # summarize_probes(probe_evaluation_filename, blast_lineage_filename, args.probe_summary_info_filename, sample_dir, probes_summary_filename, args.min_tm, args.max_tm, args.gc, args.mch, args.bot, args.bitscore_thresh, args.sod, args.dnaconc, args.mt_cutoff, args.ot_gc_cutoff)
    summarize_probes(probe_evaluation_filename, args.probe_summary_info_filename, sample_dir, probes_summary_filename, args.probe_selection_level, args.total_transcript_table_filename, args.number_of_probes, args.min_tm, args.max_tm, args.gc, args.mch, args.bot, args.bitscore_thresh, args.sod, args.dnaconc, args.mt_cutoff, args.ot_gc_cutoff)
    return

if __name__ == '__main__':
    main()
