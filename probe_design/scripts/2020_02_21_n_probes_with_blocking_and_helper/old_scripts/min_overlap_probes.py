def main():
    parser = argparse.ArgumentParser('Design FISH probes for a complex microbial community')

    parser.add_argument('probe_evaluation_complete_filename', type = str, help = 'Input file containing blast results')

    parser.add_argument('design_id', type = str, help = 'Input file containing blast lineage')

    parser.add_argument('probe_summary_info_filename', type = str, help = 'Input file containing blast lineage')

    parser.add_argument('-tmin', '--min_tm', dest = 'min_tm', type = float, default = 55.0, help = 'Boolean to indicate whether to group sequences by similarity instead of taxonomic info')

    parser.add_argument('-tmax', '--max_tm', dest = 'max_tm', type = float, default = 55.0, help = 'Boolean to indicate whether to group sequences by similarity instead of taxonomic info')

    parser.add_argument('-m', '--mch', dest = 'mch', type = int, default = 14, help = 'Boolean to indicate whether to group sequences by similarity instead of taxonomic info')

    parser.add_argument('-tpn', '--top_n_probes', dest = 'tpn', type = int, default = 1, help = 'Number of top probes to keep')

    parser.add_argument('-gc', '--gc', dest = 'gc', type = float, default = 60.0, help = 'Number of top probes to keep')

    parser.add_argument('-bot', '--bot', dest = 'bot', type = float, default = 0.0, help = 'Number of top probes to keep')

    parser.add_argument('-bt', '--bitscore_thresh', dest = 'bitscore_thresh', type = float, default = 0.001, help = 'Number of top probes to keep')

    parser.add_argument('-ommcht', '--off_target_max_mch_threshold', dest = 'off_target_max_mch_threshold', type = int, default = 9, help = 'Number of top probes to keep')

    parser.add_argument('-c', '--probe_selection_method', dest = 'probe_selection_method', type = str, default = 'AllSpecific', help = 'Probe selection method. AllSpecific (default) | Single Best | MinOverlap | Top N')

    args = parser.parse_args()

if __name__ == '__main__':
    main()
