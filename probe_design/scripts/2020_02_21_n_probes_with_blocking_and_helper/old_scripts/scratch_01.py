import pandas as pd
import os
from snapgene_reader import snapgene_file_to_dict
import string
from matplotlib import pyplot as plt

# gfp probe check
def parse_snapgene_file(snapgene_filename):
    plasmid = snapgene_file_to_dict(snapgene_filename)
    plasmid_seq = plasmid['seq']
    plasmid_annotation = pd.DataFrame(columns = ['Feature', 'Name', 'Start', 'End', 'Length', 'Sequence'])
    n_features = len(plasmid['features'])
    for i in range(n_features):
        plasmid_annotation.loc[i, 'Feature'] = i
        plasmid_annotation.loc[i, 'Name'] = plasmid['features'][i]['name']
        start = plasmid['features'][i]['start']
        end = plasmid['features'][i]['end']
        plasmid_annotation.loc[i, 'Start'] = start
        plasmid_annotation.loc[i, 'End'] = end
        plasmid_annotation.loc[i, 'Length'] = end - start
        plasmid_annotation.loc[i, 'Sequence'] = plasmid['seq'][start:end]
    return(plasmid_annotation)


snapgene_filename = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/runs/2019_12_03_gfp_plasmid_dh5alpha/input/pJKR-H-tetR_Addgene_62561.dna'
os.path.exists(snapgene_filename)
plasmid_annotation = parse_snapgene_file(snapgene_filename)
plasmid_annotation.columns
plasmid_annotation['Name']
plasmid_annotation['Feature']


gfp_best_probes_off_target_summary_filename = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/runs/simulation/DSGN_001/feature_8_off_target_summary_info.csv'
gfp_best_probes_off_target_summary = pd.read_csv(gfp_best_probes_off_target_summary_filename)
gfp_best_probes_filtered_filename = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/runs/simulation/DSGN_001/feature_8_probe_selection.csv'
gfp_best_probes_filtered = pd.read_csv(gfp_best_probes_filtered_filename)

probes_used_list = ['caccggtgttgttccgatcctgg','cgtggtgaaggtgaaggtgacgc','ctgcaccaccggtaaactgccgg','ccgtggccgaccctggtt','ccctgacctacggtgttcagtgc',
                    'cccggaccacatgaaacagcacg','cggttctgttcagctggctgacc','cgatcggtgacggtccggttctg','ggtgacggtccggttctgt','cgttaccgctgctggtatcaccc']
probes_used_list = [i.upper() for i in probes_used_list]
probes_used_list[0]
probes_used_list[0][::-1]
def rc(string):
    tab = string.maketrans("ACTG", "TGAC")
    return(string.translate(tab)[::-1])

rc(probes_used_list[0])
probes_used_list = [rc(i) for i in probes_used_list]

# probes_used_ids = [gfp_best_probes_filtered[gfp_best_probes_filtered['seq'] == i]['probe_id'] for i in probes_used_list]

probes_used_ids = gfp_best_probes_filtered[gfp_best_probes_filtered['seq'].isin(probes_used_list)]['probe_id']

# gfp_best_probes_filtered[gfp_best_probes_filtered['probe_id'] == 221]['seq']
# probes_used_list[0]

# gfp_best_probes_off_target_summary.columns
# gfp_used_probes_off_target_summary = pd.DataFrame(columns = gfp_best_probes_off_target_summary.columns)
# for i in probes_used_ids:
#     new_vals = gfp_best_probes_off_target_summary[gfp_best_probes_off_target_summary['probe_id'] == i]
#     gfp_used_probes_off_target_summary = pd.concat([gfp_used_probes_off_target_summary, new_vals])

gfp_used_probes_off_target_summary = gfp_best_probes_off_target_summary[gfp_best_probes_off_target_summary['probe_id'].isin(probes_used_ids)]

gfp_used_probes_off_target_summary.shape
gfp_used_probes_off_target_summary.columns
# gfp_used_probes_off_target_summary['sseq'].iloc[0].count('G') + gfp_used_probes_off_target_summary['sseq'].iloc[0].count('C')
gc_content = [gfp_used_probes_off_target_summary['sseq'].iloc[i].count('G')
                                                    + gfp_used_probes_off_target_summary['sseq'].iloc[i].count('C')
                                                    for i in range(gfp_used_probes_off_target_summary.shape[0])]
# gc_content = [str(gfp_used_probes_off_target_summary.loc[i, 'sseq']).count('G') + gfp_used_probes_off_target_summary.loc[i, 'sseq'].count('c')]
# gfp_used_probes_off_target_summary.loc[1, 'sseq']
plt.xlabel('Off-target GC count for GFP plasmid probes')
plt.ylabel('Frequency')
plt.hist(gc_content, bins = max(gc_content) - min(gc_content), color = 'blue')
hist_filename = '/workdir/bmg224/hiprfish/mobile_elements/experiments/2020_01_30_mixed_plasmid/figures/off_target_gc_histogram.png'
plt.savefig(hist_filename, format = 'png')
# plt.show()
plt.close()
