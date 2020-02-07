

# Integrate bot, bitscore, off target max tm, and off target gc content into selection
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

######
# bigger defs
######

# Add tm and gc count to probe_blast dataframe
for probe_idx in probes.probe_id:
    probe_name = '{}_{}/{}_{}_{}'.format(target_rank, target_taxon, target_rank, target_taxon, probe_idx)
    probe_blast = pd.read_hdf(probe_evaluation_filename, probe_name)
    probe_blast['melting_temp'] = calculate_tm(probe_blast, Na, dnac1oligo)
    probe_blast['GC_count'] = calculate_gc_count(probe_blast)
    probe_blast.loc[:,'target_taxon_hit'] = (probe_blast[target_rank].values.astype(str) == str(target_taxon))
    probe_blast.loc[:,'target_taxon_hit_full_match'] = (probe_blast[target_rank].values.astype(str) == str(target_taxon))*(probe_blast.pid.values >= 99.9)*(probe_blast.qcovhsp.values >= 99.9)
    probe_summary = probe_summary.append(probe_blast_summarize(probe_blast, max_continuous_homology = max_continuous_homology, taxon_abundance = target_taxon_abundance, target_rank = target_rank), ignore_index = True, sort = False)

# Space the probes out nicely
probes_merge['mean_probe_end'] = probes_merge.mean_probe_start.values + probes_merge.length.values - 1
probes_merge['mean_probe_start_group'] = np.floor(probes_merge.mean_probe_start.values/20).astype(int)
if probes.shape[0] > 0:
    best_probes = pd.DataFrame()
    group_distance = 120
    while (best_probes.shape[0] <= 15) & (group_distance > 20):
        group_distance -= 20
        probes_merge.loc[:,'mean_probe_start_group'] = np.floor(probes_merge.mean_probe_start.values/group_distance).astype(int)
        if probes_merge.empty:
            print('{}, {}'.format(probe_evaluation_filename, probe_summary.off_target_max_gc.min()))
        best_probes = select_all_specific_p_start_group_probes(probes_merge, group_distance, bot, bitscore_thresh, mt_cutoff, ot_gc_cutoff)
    if best_probes.shape[0] > 15:
        group_distance += 20
        probes_merge['mean_probe_start_group'] = np.floor(probes_merge.mean_probe_start.values/group_distance).astype(int)
        best_probes = select_all_specific_p_start_group_probes(probes_merge, group_distance, bot, bitscore_thresh, mt_cutoff, ot_gc_cutoff)
    elif best_probes.shape[0] == 0:
        probes_sorted = probes_merge.sort_values(by = 'blast_on_target_rate', ascending = False)
        best_probes = pd.DataFrame(probes_sorted.iloc[[0],:])
        best_probes.loc[:,'selection_method'] = 'SingleBest'

# Blocking probe stuff
def generate_blocking_probes(design_dir, bplc, target_rank, plf = 'T', primer = 'T', primerset = 'B', barcode_selection = 'MostSimple'):
    design_dir_folder = os.path.split(design_dir)[1]
    blocking_probes_filenames = glob.glob('{}/*_off_target_summary_info.csv'.format(design_dir))
    probes_filename = '{}/{}_primerset_{}_barcode_selection_{}_full_length_probes.csv'.format(design_dir, design_dir_folder, primerset, barcode_selection)
    probes = pd.read_csv(probes_filename)
    if '{}/0_off_target_summary_info.csv'.format(design_dir) in blocking_probes_filenames:
        blocking_probes_filenames.remove('{}/0_off_target_summary_info.csv'.format(design_dir))
    blocking_probes_list = []
    for f in blocking_probes_filenames:
        target_taxon = re.sub('_off_target_summary_info.csv', '', os.path.basename(f))
        if not probes.loc[probes.target_taxon.values.astype(str) == target_taxon, :].empty:
            off_target_summary_full = pd.read_csv(f)
            if not off_target_summary_full.empty:
                off_target_summary_full.loc[:, 'max_average_encoding_interference_fraction_0bp'] = 0.0
                off_target_summary_full.loc[:, 'max_average_encoding_interference_fraction_5bp'] = 0.0
                off_target_summary_full.loc[:, 'max_average_encoding_interference_fraction_10bp'] = 0.0
                for bp in off_target_summary_full.probe_id.drop_duplicates().values:
                    off_target_summary = off_target_summary_full.loc[off_target_summary_full.probe_id.values == bp, :]
                    off_target_summary.loc[:, 'average_encoding_interference_fraction_0bp'] = 0.0
                    off_target_summary.loc[:, 'average_encoding_interference_fraction_5bp'] = 0.0
                    off_target_summary.loc[:, 'average_encoding_interference_fraction_10bp'] = 0.0
                    blocked_taxa = off_target_summary[target_rank].drop_duplicates()
                    for tt in blocked_taxa.values:
                        off_target_summary_blocked_taxon = off_target_summary.loc[off_target_summary[target_rank].values == tt, :]
                        off_target_summary_p_start = off_target_summary_blocked_taxon.loc[:, ['molecule_start', 'molecule_end']].drop_duplicates()
                        taxon_probes = probes.loc[probes.target_taxon.values == tt, :]
                        encoding_p_start = taxon_probes.loc[:,['mean_probe_start', 'length']].drop_duplicates()
                        if not encoding_p_start.empty:
                            significant_overlap_fraction = np.zeros((encoding_p_start.shape[0], 3))
                            for i in range(encoding_p_start.shape[0]):
                                min_end = off_target_summary_p_start.molecule_start.apply(min, args = (encoding_p_start.mean_probe_start.values[i] + encoding_p_start.length.values[i], ))
                                max_start = off_target_summary_p_start.molecule_end.apply(max, args = (encoding_p_start.mean_probe_start.values[i], ))
                                overlap = min_end - max_start
                                significant_overlap_fraction[i, 0] = np.sum(overlap > 0)/overlap.shape[0]
                                significant_overlap_fraction[i, 1] = np.sum(overlap > 5)/overlap.shape[0]
                                significant_overlap_fraction[i, 2] = np.sum(overlap > 10)/overlap.shape[0]
                            off_target_summary.loc[off_target_summary[target_rank].values == tt, 'average_encoding_interference_fraction_0bp'] = np.average(significant_overlap_fraction[:,0])
                            off_target_summary.loc[off_target_summary[target_rank].values == tt, 'average_encoding_interference_fraction_5bp'] = np.average(significant_overlap_fraction[:,1])
                            off_target_summary.loc[off_target_summary[target_rank].values == tt, 'average_encoding_interference_fraction_10bp'] = np.average(significant_overlap_fraction[:,2])
                    off_target_summary_full.loc[off_target_summary_full.probe_id.values == bp, 'max_average_encoding_interference_fraction_0bp'] = off_target_summary.average_encoding_interference_fraction_0bp.max()
                    off_target_summary_full.loc[off_target_summary_full.probe_id.values == bp, 'max_average_encoding_interference_fraction_5bp'] = off_target_summary.average_encoding_interference_fraction_5bp.max()
                    off_target_summary_full.loc[off_target_summary_full.probe_id.values == bp, 'max_average_encoding_interference_fraction_10bp'] = off_target_summary.average_encoding_interference_fraction_10bp.max()
            blocking_probes_list.append(off_target_summary_full)
    blocking_probes = pd.concat(blocking_probes_list, sort = False)
    blocking_probes = blocking_probes.loc[blocking_probes.length.values >= bplc, :]
    blocking_probes.to_csv('{}/{}_primerset_{}_barcode_selection_{}_full_length_blocking_probes.csv'.format(design_dir, design_dir_folder, primerset, barcode_selection), index = False)
    blocking_probes_order_format = blocking_probes[['sseq', 'length']]
    blocking_probes_order_format.columns = ['Sequence', 'Length']
    blocking_probes_order_format = blocking_probes_order_format.drop_duplicates()
    blocking_probes_order_format = blocking_probes_order_format.sort_values(by = 'Length', ascending = False)
    # blocking_probes_order_format = blocking_probes_order_format.loc[blocking_probes_order_format.length.values >= bplc]
    blocking_probes_order_format.insert(0, 'Blocking_Probe_Name', ['BP_{}'.format(i) for i in range(blocking_probes_order_format.shape[0])])
    blocking_probes_order_format = blocking_probes_order_format.assign(Amount = '25nm', Purification = 'STD')
    blocking_probes_order_format.loc[blocking_probes_order_format.Length.values < 15, 'Amount'] = '100nm'
    blocking_probes_order_format = blocking_probes_order_format[['Blocking_Probe_Name', 'Sequence', 'Amount', 'Purification', 'Length']]
    blocking_probes_order_format.to_excel('{}/{}_primerset_{}_barcode_selection_{}_full_length_blocking_probes_order_format.xlsx'.format(design_dir, design_dir_folder, primerset, barcode_selection))
