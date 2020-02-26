
"""
Collect HiPRFISH probe design results
Hao Shi 2017
"""

import os
import argparse
import pandas as pd
import numpy as np
# from Bio import SeqIO


###############################################################################################################
# HiPR-FISH : collect probe design results
###############################################################################################################

def collect_probe_coverage_results(data_dir, simulation_table, pipeline_version, output_filename):
    print('Loading samples table: {}'.format(simulation_table))
    sim_tab = pd.read_csv(simulation_table)
    sim_tab['COVERED_FEATURE_RICHNESS'] = np.nan
    sim_tab['COVERED_FEATURE_RICHNESS_FRACTION'] = np.nan
    print('Loading result files:')
    for i in range(0, sim_tab.shape[0]):
        sample = sim_tab.SAMPLE[i]
        feature_best_probes_filename = '{}/simulation/{}/feature_best_probes.csv'.format(data_dir, sim_tab.DESIGN_ID[i])
        feature_best_probes_filtered_filename = '{}/simulation/{}/feature_best_probes_filtered.csv'.format(data_dir, sim_tab.DESIGN_ID[i])
        if os.path.exists(feature_best_probes_filename):
            feature_best_probes = pd.read_csv(feature_best_probes_filename)
            feature_best_probes_filtered = pd.read_csv(feature_best_probes_filtered_filename)
            sim_tab.loc[i, 'COVERED_FEATURE_RICHNESS'] = feature_best_probes_filtered.target_feature.drop_duplicates().shape[0]
            sim_tab.loc[i, 'COVERED_FEATURE_RICHNESS_FRACTION'] = feature_best_probes_filtered.target_feature.drop_duplicates().shape[0]/feature_best_probes.target_feature.drop_duplicates().shape[0]
            print('Saving collected results to %s...' % (output_filename))
        else:
            print('Sample result file {} does not exist'.format(taxon_best_probes_filename))
        sim_tab.loc[i, 'PIPELINE_VERSION'] = pipeline_version
        sim_tab.to_csv(output_filename, index = False, header = True)
    return(sim_tab)

###############################################################################################################
# main function
###############################################################################################################

def main():
    parser = argparse.ArgumentParser('Collect summary statistics of HiPRFISH probes for a complex microbial community')

    # data directory
    parser.add_argument('data_dir', type = str, help = 'Directory of the data files')

    # input simulation table
    parser.add_argument('simulation_table', type = str, help = 'Input csv table containing simulation information')

    parser.add_argument('pipeline_version', type = str, help = 'Input csv table containing simulation information')

    # output simulation results table
    parser.add_argument('simulation_results', type = str, help = 'Output csv table containing simulation results')

    args = parser.parse_args()

    sim_tab = collect_probe_coverage_results(args.data_dir, args.simulation_table, args.pipeline_version, args.simulation_results)

if __name__ == '__main__':
    main()
