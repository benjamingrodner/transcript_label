# Load best probes
best_probes_filename = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/runs/2020_02_13_gfp_plasmid/simulation/DSGN_002/GFP_probe_selection.csv'
best_probes = pd.read_csv(best_probes_filename)
best_probes.shape
# Get all probes filtered by tm and gc
probes_merge_filename = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/runs/2020_02_13_gfp_plasmid/simulation/DSGN_002/GFP_all_probes.csv'
probes_merge = pd.read_csv(probes_merge_filename)
probes_merge.shape
# Get all probes
probes_filename = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/runs/2020_02_13_gfp_plasmid/sample_002/primer3/GFP.int'
probes = pd.read_csv(probes_filename, skiprows = 3, header = None, delim_whitespace = True)
target_feature = 'GFP'
probes.columns = ['probe_id', 'seq', 'p_start', 'ln', 'N', 'GC', 'Tm', 'self_any_th', 'self_end_th', 'hairpin', 'quality']
probes.loc[:,'target_feature'] = target_feature
# Get the start and end of the best probe
p_start = best_probes.loc[best_probes.probe_id == probe_idx, 'p_start'].values[0]
p_end = best_probes.loc[best_probes.probe_id == probe_idx, 'p_end'].values[0]
# Add a column to all summarized probes for the location of the end of the probes
probes.loc[:,'p_end'] = probes.p_start + probes.ln
# probes_merge.loc[:,'p_end'] = probes_merge.p_start + probes_merge.ln
# Design helper probes that go on either end of the best probe
five_prime_helpers = probes.loc[(probes.p_start.values > p_start - 103) & (probes.p_end.values < p_start - 3), :]
three_prime_helpers = probes.loc[(probes.p_start.values > p_end + 3) & (probes.p_end.values < p_end + 103), :]
# five_prime_helpers = probes_merge.loc[(probes_merge.p_start.values > p_start - 103) & (probes_merge.p_end.values < p_start - 3), :]
# three_prime_helpers = probes_merge.loc[(probes_merge.p_start.values > p_end + 3) & (probes_merge.p_end.values < p_end + 103), :]
helper_probes_temp = pd.concat([five_prime_helpers,three_prime_helpers])
# Pick a probe
probe_idx = 0
# Get helper probes
helper_probes = pd.DataFrame(columns = probes_merge.columns)
helper_probes['p_end'] = []
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
# HDF stuff
probe_evaluation_filename = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/runs/2020_02_13_gfp_plasmid/sample_003/blast/GFP_probe_evaluation.h5'
probe_name = 'probe_' + str(probe_idx)
probe_blast = pd.read_hdf(probe_evaluation_filename, probe_name)
probe_blast
probe_blast['melting_temp'] = calculate_tm(probe_blast, Na, dnac1oligo)
probe_blast['GC_count'] = calculate_gc_count(probe_blast)
probe_blast.loc[:,'target_feature_hit'] = (probe_blast.feature.values.astype(str) == str(target_feature))
# Create a list of the off-target bindings
probe_blast_off_target = probe_blast[probe_blast.feature.values.astype(str) != str(target_feature)]
max_continuous_homology = 10
bitscore_thresh = 25
mt_cutoff = 60
ot_gc_cutoff = 7
# Filter the probes to a temporary DataFrame using blast on target rate and off-target {bitscore, max melting temperature, and GC content}
probe_blast_off_target_worst = probe_blast_off_target.loc[(probe_blast_off_target['mch'] > max_continuous_homology) |
                                                        (probe_blast_off_target['bitscore'] > bitscore_thresh) |
                                                        (probe_blast_off_target['melting_temp'] > mt_cutoff) |
                                                        (probe_blast_off_target['GC_count'] > ot_gc_cutoff)]
