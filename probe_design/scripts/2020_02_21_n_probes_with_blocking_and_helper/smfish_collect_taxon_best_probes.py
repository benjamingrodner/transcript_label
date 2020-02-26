"""
Collect HiPRFISH probe design results
Hao Shi 2017
"""

import argparse
import os
import re
import glob
import pandas as pd
import numpy as np
from Bio import SeqIO

###############################################################################################################
# HiPR-FISH : collect probe design results
###############################################################################################################

def collect_taxon_best_probes(design_directory, sim_input_filename, feature_best_probes_filename, feature_best_probes_filtered_filename, output_probes_summary_filename, bot):
    # Pull out the design id for the sample
    simulation_directory, design_id = os.path.split(design_directory)
    # Pull the simulation run directory
    data_dir = os.path.split(simulation_directory)[0]
    # Read the simulatino table
    sim_tab = pd.read_csv(sim_input_filename)
    # Read the sample value from the simulation table
    sample = sim_tab[sim_tab.DESIGN_ID == design_id].SAMPLE.values[0]
    # Create an evaluation directory in the sample directory
    feature_evaluation_directory = '{}/{}/evaluation/'.format(data_dir, sample)
    # Get all of the filenames of feature probe selections
    feature_evaluation_filename_list = glob.glob('{}/*_probe_selection.csv'.format(design_directory))
    # Read in a list of tables containing probe design results
    best_probes_list = [pd.read_csv(filename) for filename in feature_evaluation_filename_list]
    # Sort the tables by ascending quality, as assigned on probe design
    best_probes_quality_sorted_list = [df.sort_values(by = ['quality'], ascending = True) for df in best_probes_list]
    # Combine all the tables together
    best_probes_df = pd.concat(best_probes_quality_sorted_list)
    # Get the counts of probes for each feature
    best_probes_summary = best_probes_df['target_feature'].value_counts().reset_index()
    best_probes_summary.columns = ['FEATURE', 'PROBE_COUNTS']
    # filter the best probes by blast on target rate
    best_probes_filtered = best_probes_df[best_probes_df['blast_on_target_rate'] > bot]
    # Get the filtered counts of probes for each feature
    best_probes_filtered_summary = best_probes_filtered['target_feature'].value_counts().reset_index()
    best_probes_filtered_summary.columns = ['FEATURE', 'PROBE_COUNTS']
    # Write the quality sorted list, the filtered list, and the summary info to files 
    best_probes_df.to_csv(feature_best_probes_filename, index = False)
    best_probes_filtered.to_csv(feature_best_probes_filtered_filename, index = False)
    best_probes_filtered_summary.to_csv(output_probes_summary_filename)
    return

###############################################################################################################
# main function
###############################################################################################################

def main():
    parser = argparse.ArgumentParser('Collect summary statistics of HiPRFISH probes for a complex microbial community')

    # data directory
    parser.add_argument('design_directory', type = str, help = 'Directory of the data files')

    # input simulation table
    parser.add_argument('sim_input_filename', type = str, help = 'Input csv table containing simulation information')

    # output simulation results table
    parser.add_argument('feature_best_probes_filename', type = str, help = 'Output csv table containing simulation results')

    # output simulation results table
    parser.add_argument('feature_best_probes_filtered_filename', type = str, help = 'Output csv table containing simulation results')

    # output simulation results table
    parser.add_argument('output_probes_summary_filename', type = str, help = 'Output csv table containing simulation results')

    parser.add_argument('bot', type = float, help = 'Output csv table containing simulation results')

    args = parser.parse_args()

    collect_taxon_best_probes(args.design_directory, args.sim_input_filename, args.feature_best_probes_filename, args.feature_best_probes_filtered_filename, args.output_probes_summary_filename, args.bot)

if __name__ == '__main__':
    main()
