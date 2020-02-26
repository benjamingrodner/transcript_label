import argparse
import pandas as pd
import subprocess
import os
import multiprocessing
import glob
import re
import itertools
import numpy as np
import random
from ete3 import NCBITaxa
# from SetCoverPy import setcover
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, generic_dna
from Bio.Blast.Applications import NcbiblastnCommandline
# from joblib import Parallel, delayed
# from snapgene_reader import snapgene_file_to_dict

###############################################################################################################
# HiPR-FISH : design probes
###############################################################################################################

###############################################################################################################
# Workflow functions
###############################################################################################################

# parse plasmid annotation files
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

def parse_genomic_information(gtf_filename, genome_fasta_filename):
    gtf = pr.read_gtf(gtf_filename)
    chromosomes = list(SeqIO.parse(genome_fasta_filename), 'fasta')
    return

# write individual feature fasta files
def write_individual_feature_fasta(plasmid_annotation, length_threshold, output_dir):
    plasmid_annotation_filtered = plasmid_annotation.loc[plasmid_annotation.Length.values >= length_threshold, :]
    plasmid_annotation_filtered = plasmid_annotation_filtered.reset_index().drop(columns = ['index'])
    plasmid_annotation_filtered.to_csv('{}/plasmid_annotation_filtered.csv'.format(output_dir), index = None)
    for i in range(plasmid_annotation_filtered.shape[0]):
        output_filename = '{}/feature_{}.fasta'.format(output_dir, plasmid_annotation_filtered.loc[i, 'Feature'])
        with open(output_filename, 'w') as output_handle:
            feature_sequence = SeqRecord(Seq(plasmid_annotation_filtered.loc[i, 'Sequence']).reverse_complement(), id = 'feature_{}'.format(plasmid_annotation_filtered.loc[i, 'Feature']), description = '')
            SeqIO.write(feature_sequence, output_handle, 'fasta')
    return

# write commbined features fasta files
def write_feature_fasta(plasmid_annotation, length_threshold, output_filename):
    try:
        os.remove(output_filename)
    except:
        pass
    plasmid_annotation_filtered = plasmid_annotation.loc[plasmid_annotation.Length.values >= length_threshold, :]
    plasmid_annotation_filtered = plasmid_annotation_filtered.reset_index().drop(columns = ['index'])
    plasmid_seqs = [SeqRecord(Seq(plasmid_annotation_filtered.loc[i, 'Sequence']).reverse_complement(), id = 'feature_{}'.format(plasmid_annotation_filtered.loc[i, 'Feature']), description = '') for i in range(plasmid_annotation_filtered.shape[0])]
    with open(output_filename, 'w') as output_handle:
        SeqIO.write(plasmid_seqs, output_handle, 'fasta')
    return

def combine_feature_and_custom_blast_database(feature_fasta, custom_blast_fasta):
    combined_seqs = []
    feature_seqs = SeqIO.parse(feature_fasta, 'fasta')
    for record in feature_seqs:
        combined_seqs.append(record)
    custom_blast_seqs = SeqIO.parse(custom_blast_fasta, 'fasta')
    for record in custom_blast_seqs:
        combined_seqs.append(record)
    dir = os.path.split(custom_blast_fasta)[0]
    output_filename = '{}/custom_blast_db.fasta'.format(dir)
    with open(output_filename, 'w') as output_handle:
        SeqIO.write(combined_seqs, output_handle, 'fasta')
    return(output_filename)

def combine_database_fastas(database_fastas, database_dir):
    combined_seqs = []
    # database_fastas = glob.glob('{}/*.fasta'.format(database_dir))
    # database_fastas = database_fastas + glob.glob('{}/*.fna'.format(database_dir))
    for database in database_fastas:
        database_seqs = SeqIO.parse(database, 'fasta')
        for record in database_seqs:
            combined_seqs.append(record)
    # for record in feature_seqs:
    #     combined_seqs.append(record)
    # custom_blast_seqs = SeqIO.parse(custom_blast_fasta, 'fasta')
    # for record in custom_blast_seqs:
    #     combined_seqs.append(record)
    # dir = os.path.split(custom_blast_fasta)[0]
    output_filename = '{}/custom_blast_db.fasta'.format(database_dir)
    with open(output_filename, 'w') as output_handle:
        SeqIO.write(combined_seqs, output_handle, 'fasta')
    return(output_filename)

def make_blast_db(input_fasta_file):
    subprocess.call(['makeblastdb','-in', input_fasta_file, '-parse_seqids', '-dbtype', 'nucl'])
    return


def make_blast_db_from_csv(blast_db_csv_filename):
    # Read in the csv file to a pandas table
    blast_db_table = pd.read_csv(blast_db_csv_filename)
    # Drop all duplicate nodes
    blast_db_table_reduced = blast_db_table.drop_duplicates(subset = 'Gene')
    # Pull the sequence node names and sequences
    blast_db_fasta = [SeqRecord(Seq(str(blast_db_table_reduced.iloc[i, :]['Sequence'])), id = blast_db_table_reduced.iloc[i, :]['Gene'], description = '') for i in range(blast_db_table_reduced.shape[0])]
    # Get directory for blast database
    blast_db_dir = os.path.split(blast_db_csv_filename)[0]
    # Make a filename for the custom blast db
    blast_db_fasta_filename = '{}/custom_blast_db.fasta'.format(blast_db_dir)
    # Write sequences and names to file
    with open(blast_db_fasta_filename, 'w') as output_handle:
        SeqIO.write(blast_db_fasta, output_handle, 'fasta')
    # Make a blast database in the file
    # make_blast_db(blast_db_fasta_filename)
    return(blast_db_fasta_filename)


def generate_reverse_complement_fasta(gene_sequences_filename, gene_sequences_rc_filename):
    # Parse the target gene fasta
    gene_sequences = SeqIO.parse(gene_sequences_filename, 'fasta')
    # Generate a reverse compliment with the same id
    gene_sequences_rc = [SeqRecord(record.seq.reverse_complement(), id = record.id, description = '') for record in gene_sequences]
    # Save the new file
    with open(gene_sequences_rc_filename, 'w') as output_handle:
        SeqIO.write(gene_sequences_rc, output_handle, 'fasta')
    return


# probe design
def probe_design(feature_sequences_filename, outdir, include_start, include_end):
    # running the function
    settings_filename = outdir + '/primer3_settings.txt'
    primer3_input_filename = outdir + '/features_primer3_input.txt'
    output_filename= outdir + '/feature_probes.fasta'

    # write primer3 settings file
    primer3_settings = ['Primer3 File - http://primer3.sourceforge.net',
                        'P3_FILE_TYPE=settings',
                        '',
                        'P3_FILE_ID=FISH probe design',
                        'P3_FILE_FLAG=1',
                        'PRIMER_FIRST_BASE_INDEX=1',
                        'PRIMER_TASK=generic',
                        'PRIMER_EXPLAIN_FLAG=1',
                        'PRIMER_NUM_RETURN=10000',
                        'PRIMER_PICK_LEFT_PRIMER=0',
                        'PRIMER_PICK_INTERNAL_OLIGO=1',
                        'PRIMER_PICK_RIGHT_PRIMER=0',
                        'PRIMER_INTERNAL_OPT_SIZE=20',
                        'PRIMER_INTERNAL_MIN_SIZE=18',
                        'PRIMER_INTERNAL_MAX_SIZE=23',
                        'PRIMER_INTERNAL_MIN_TM=' + str(40),
                        'PRIMER_INTERNAL_MAX_SELF_ANY_TH=1000.00',
                        'PRIMER_INTERNAL_MAX_HAIRPIN_TH=1000.0',
                        'PRIMER_INTERNAL_MAX_NS_ACCEPTED=0',
                        'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1',
                        'PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=0',
                        'PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/programs/primer3-2.3.5/src/primer3_config/',
                        'PRIMER_LOWERCASE_MASKING=0',
                        'PRIMER_PICK_ANYWAY=1',
                        '=']
    with open(settings_filename, 'w') as handle:
        for l in primer3_settings:
            _ = handle.write('{}\n'.format(l))
    # write primer3 input file
    feature_seqs = SeqIO.parse(feature_sequences_filename, 'fasta')
    for record in feature_seqs:
        primer3_record = ['SEQUENCE_ID=' + str(record.id), 'SEQUENCE_TEMPLATE=' + str(record.seq).upper(), 'SEQUENCE_INCLUDED_REGION=' + str(include_start) + ',' + str(len(record.seq) - include_end - include_start), 'P3_FILE_FLAG=1', 'PRIMER_EXPLAIN_FLAG=1', '=']
        pd.DataFrame(primer3_record).to_csv(primer3_input_filename, header = None, index = False, mode = 'a', sep = ' ')

    os.chdir(outdir)
    if os.path.exists(output_filename):
        print('No probe design was rerun')
    else:
        print('Probe design was run')
        # subprocess.check_call(['primer3_core', '-p3_settings_file', settings_filename, '-output', output_filename, '-format_output', primer3_input_filename])
        subprocess.check_call(['/programs/primer3-2.3.5/src/primer3_core', '-p3_settings_file', settings_filename, '-output', output_filename, '-format_output', primer3_input_filename])
    return

def split_feature_probe_file(probe_filename):
# def split_feature_probe_file(probe_filename, molecule):
    probe_dir, feature_filename = os.path.split(probe_filename)
    feature = re.sub('.int', '', feature_filename)
    feature_directory = '{}/{}'.format(probe_dir, feature)
    if not os.path.exists(feature_directory):
        os.makedirs(feature_directory)
    probes = pd.read_table(probe_filename, skiprows = 3, header = None, delim_whitespace = True)
    probes.columns = ['probe_num', 'seq', 'start', 'length', 'N', 'GC', 'Tm', 'self_any_th', 'self_end_th', 'hair-pin', 'quality']
    # if molecule == 'RNA':
    #     probes_list = [SeqRecord(Seq(probes['seq'][i]).reverse_complement(), id = str(probes['probe_num'][i]), description = '') for i in range (0, probes.shape[0])]
    # elif molecule == 'DNA':
    #     probes_list = [SeqRecord(Seq(probes['seq'][i]), id = str(probes['probe_num'][i]), description = '') for i in range (0, probes.shape[0])]
    probes_list = [SeqRecord(Seq(probes['seq'][i]), id = str(probes['probe_num'][i]), description = '') for i in range (0, probes.shape[0])]
    probes_fasta_filenames = ['{}/{}/{}_probe_{}.fasta'.format(probe_dir, feature, feature, probe_num) for probe_num in range(probes.shape[0])]
    for i in range(probes.shape[0]):
        SeqIO.write(probes_list[i], probes_fasta_filenames[i], 'fasta')

def rename_fasta_ids(fasta_filename, target, substitute):
    fasta = SeqIO.parse(fasta_filename, "fasta")
    fasta_new = []
    for record in fasta:
        # print(re.sub('\|','_',record.description))
        # print(record.id)
        f = SeqRecord(record.seq)
        f.id = re.sub(target, substitute,record.id)
        f.description = re.sub(target, substitute,record.description)
        fasta_new.append(f)
    # combined_seqs
    dir, name = os.path.split(fasta_filename)
    name = re.sub('.fasta|.fna','',name)
    # name = re.sub('.fna','',name)
    output_filename = '{}/{}_renamed.fasta'.format(dir, name)
    # combined_filename
    with open(output_filename, 'w') as output_handle:
        SeqIO.write(fasta_new, output_handle, 'fasta')
    return(output_filename)

def get_gene_sequences_filename(data_dir):
    gene_sequences_filename = glob.glob('{}/inputs/*.fasta'.format(data_dir))[0]
    return(gene_sequences_filename)

def get_target_table_filename(data_dir):
    target_table_filename = glob.glob('{}/inputs/target_table*.csv'.format(data_dir))[0]
    return(target_table_filename)

def get_specific_gene_id(fasta_filename, target_gene_names):
    fasta = SeqIO.parse(fasta_filename, "fasta")
    target_gene_ids = pd.DataFrame()
    for record in fasta:
        description = record.description
        gene_name = description.split(' ')[1]
        gene_name = re.sub(r"gene=|\[|\]", '', gene_name)
        if gene_name in target_gene_names:
            target_gene_ids = target_gene_ids.append({'name': gene_name, 'id': record.id}, ignore_index = True)
    return(target_gene_ids)


###############################################################################################################
# main function
###############################################################################################################


def main():
    parser = argparse.ArgumentParser('Design FISH probes for a complex microbial community')
    # Target sequences
    parser.add_argument('sample_dir', type = str, help = 'sample dir containing inputs etc')
    # parser.add_argument('gene_sequences_filename', type = str, help = 'Input FASTA file containing target sequences')
    # parser.add_argument('snapgene_filename', type = str, help = 'Input FASTA file containing 16S sequences')
    # parser.add_argument('custom_blast_db', type = str, help = 'Blast database')
    # Csv containing information from hao Zhou
    parser.add_argument('custom_blast_db', type = str, help = 'Blast database')
    # parser.add_argument('blast_db_csv_filename', type = str, help = 'Blast database')
    parser.add_argument('-l', '--length_threshold', dest = 'length_threshold', type = int, default = 50, help = 'Length threshold')
    parser.add_argument('-b', '--include_start', dest = 'include_start', type = int, default = 1, help = 'Starting position of included region for probe design, measured from the beginning of 16S rRNA molecules')
    parser.add_argument('-e', '--include_end', dest = 'include_end', type = int, default = 1, help = 'Ending position of included region for probe design, measured from the end of 16S rRNA molecules')
    parser.add_argument('-fnid', '--feature_in_db', dest = 'feature_in_db', type = bool, default = False, help = 'Is the target feature in the database?')
    parser.add_argument('-mol', '--molecule', dest = 'molecule', type = str, default = 'RNA', help = 'Is the target DNA or RNA (i.e. are we hybridizing to non-coding or coding)?')
    parser.add_argument('-sm', '--select_markers', dest = 'select_markers', type = str, default = 'False', help = 'Do we want to select probes for marker genes ')

    args = parser.parse_args()

    # sim_input_dir = os.path.split(args.gene_sequences_filename)[0]
    # sim_input_dir = os.path.split(args.snapgene_filename)[0]
    # sim_dir = os.path.split(args.sample_dir)[0]
    # sim_primer3_dir = '{}/primer3'.format(sim_dir)
    sim_primer3_dir = '{}/primer3'.format(args.sample_dir)

    # plasmid_annotation = parse_snapgene_file(args.snapgene_filename)
    # feature_sequences_filename = '{}/features.fasta'.format(sim_primer3_dir)
    # write_individual_feature_fasta(plasmid_annotation, args.length_threshold, sim_primer3_dir)
    # write_feature_fasta(plasmid_annotation, args.length_threshold, feature_sequences_filename)
    if os.path.isdir(args.custom_blast_db):
        blast_db_fasta_filename = '{}/custom_blast_db.fasta'.format(args.custom_blast_db)
        if not os.path.exists(blast_db_fasta_filename):
            database_fastas = glob.glob('{}/*.fasta'.format(args.custom_blast_db))
            database_fastas = database_fastas + glob.glob('{}/*.fna'.format(args.custom_blast_db))
            custom_blast_db_renamed_list = []
            for database in database_fastas:
                database_name = re.sub('.fasta|.fna', '',os.path.split(database)[1])
                custom_blast_db_renamed = rename_fasta_ids(database, '\|', '_')
                custom_blast_db_renamed_list.append(custom_blast_db_renamed)
            blast_db_fasta_filename = combine_database_fastas(custom_blast_db_renamed_list, args.custom_blast_db)
            # print(args.feature_in_db)
# # Temrinal
# custom_blast_db = '/workdir/bmg224/hiprfish/transcript_label/probe_design/database/cg_ec_coding'
# database_fastas = glob.glob('{}/*.fasta'.format(custom_blast_db))
# database_fastas
# custom_blast_db_renamed_list = []
# for database in database_fastas:
#     # database_name = re.sub('.fasta|.fna', '',os.path.split(database)[1])
#     custom_blast_db_renamed = rename_fasta_ids(database, '\|', '_')
#     custom_blast_db_renamed_list.append(custom_blast_db_renamed)
# combined_blast_db_filename = combine_database_fastas(custom_blast_db_renamed_list, custom_blast_db)

    # combined_blast_db_filename = combine_feature_and_custom_blast_database(args.gene_sequences_filename, args.custom_blast_db)
            if not args.feature_in_db:
                # custom_blast_db_renamed = rename_fasta_ids(combined_blast_db_filename, '\|', '_')
                gene_sequences_filename = get_gene_sequences_filename(args.sample_dir)
                blast_db_fasta_filename = combine_feature_and_custom_blast_database(gene_sequences_filename, blast_db_fasta_filename)
                # make_blast_db(combined_blast_db_filename)
            # else:
                # custom_blast_db_final = rename_fasta_ids(combined_blast_db_filename, '\|', '_')
                # make_blast_db(custom_blast_db_renamed)
        if not os.path.exists(blast_db_fasta_filename + '.nhr'):
            make_blast_db(blast_db_fasta_filename)
    else:
        blast_db_dir = os.path.split(args.custom_blast_db)[0]
        blast_db_fasta_filename = '{}/custom_blast_db.fasta'.format(blast_db_dir)
        if not os.path.exists(blast_db_fasta_filename):
            if '.csv' in args.custom_blast_db:
                # blast_db_dir = os.path.split(args.custom_blast_db)[0]
                # blast_db_fasta_filename = '{}/custom_blast_db.fasta'.format(blast_db_dir)
                blast_db_fasta_filename = make_blast_db_from_csv(args.custom_blast_db)
                if not args.feature_in_db:
                    gene_sequences_filename = get_gene_sequences_filename(args.sample_dir)
                    blast_db_fasta_filename = combine_feature_and_custom_blast_database(gene_sequences_filename, custom_blast_db_from_csv)
            elif '.fasta' in args.custom_blast_db or '.fna' in args.custom_blast_db:
                if not args.feature_in_db:
                    custom_blast_db_renamed = rename_fasta_ids(args.custom_blast_db, '\|', '_')
                    gene_sequences_filename = get_gene_sequences_filename(args.sample_dir)
                    blast_db_fasta_filename = combine_feature_and_custom_blast_database(gene_sequences_filename, custom_blast_db_renamed)
                    # make_blast_db(combined_blast_db_filename)
                else:
                    blast_db_fasta_filename = rename_fasta_ids(args.custom_blast_db, '\|', '_')
                    # make_blast_db(custom_blast_db_renamed)
        if not os.path.exists(blast_db_fasta_filename + '.nhr'):
            make_blast_db(blast_db_fasta_filename)

    # # Get directory for blast database
    # blast_db_dir = os.path.split(args.blast_db_csv_filename)[0]
    # # Make a filename for the custom blast db
    # blast_db_fasta_filename = '{}/custom_blast_db.fasta'.format(blast_db_dir)
    # # Pull the sequence and node name from the csv of total transcripts, put in a fasta file, and make a blast database
    # if not os.path.exists(blast_db_fasta_filename):
    #     make_blast_db_from_csv(args.blast_db_csv_filename)

    # Get the input_dir and file end
    # input_dir, filename_end = os.path.split(args.gene_sequences_filename)
    # # Create a filename for the reverse complement probe design fasta
    # gene_sequences_rc_filename = input_dir + '/rc/' + re.sub('.fasta', '_rc.fasta', filename_end)
    # # make a reverse complement fasta for the target genes
    # if not os.path.exists(gene_sequences_rc_filename):
    #     os.makedirs(input_dir + '/rc')
    #     generate_reverse_complement_fasta(args.gene_sequences_filename, gene_sequences_rc_filename)

    # Design probes on the reverse complement sequences
    sample = os.path.split(args.sample_dir)[1]
    gene_sequences_filename = '{}/inputs/{}.fasta'.format(args.sample_dir, sample)
    if not args.feature_in_db:
        gene_sequences_filename = get_gene_sequences_filename(args.sample_dir)
    elif not os.path.exists(gene_sequences_filename):
        dict = SeqIO.to_dict(SeqIO.parse(blast_db_fasta_filename, "fasta"))
        target_table_filename = get_target_table_filename(args.sample_dir)
        target_table = pd.read_csv(target_table_filename)
        gene_sequences = []
        for index, target in target_table.iterrows():
            target_name = target.target_name
            gene_id = get_specific_gene_id(blast_db_fasta_filename, [target_name]).id.values[0]
            target_table.loc[index, 'target_feature'] = gene_id
            sequence = dict[gene_id].seq
            record = SeqRecord(sequence, id = gene_id)
            gene_sequences.append(record)
        target_table.to_csv(target_table_filename)
        # sample = os.path.split(args.sample_dir)[1]
        # gene_sequences_filename = '{}/inputs/{}.fasta'.format(args.sample_dir, sample)
        with open(gene_sequences_filename, 'w') as output_handle:
            SeqIO.write(gene_sequences, output_handle, 'fasta')
    if args.molecule == 'RNA':
        gene_sequences_rc_filename = '{}_rc.fasta'.format(gene_sequences_filename)
        if not os.path.exists(gene_sequences_rc_filename):
            generate_reverse_complement_fasta(gene_sequences_filename, gene_sequences_rc_filename)
        probe_design(gene_sequences_rc_filename, sim_primer3_dir, args.include_start, args.include_end)
        feature_probes_filename = glob.glob('{}/*.int'.format(sim_primer3_dir))
        for f in feature_probes_filename:
            split_feature_probe_file(f)
            # if not os.path.exists(sim_dir + '/blast'):
            #     os.makedirs(sim_dir + '/blast')
    elif args.molecule == 'DNA':
        probe_design(gene_sequences_filename, sim_primer3_dir, args.include_start, args.include_end)
        feature_probes_filename = glob.glob('{}/*.int'.format(sim_primer3_dir))
        for f in feature_probes_filename:
            split_feature_probe_file(f)
        # if not os.path.exists(sim_dir + '/blast'):
        #     os.makedirs(sim_dir + '/blast')




    # Write marker fasta and list of marker names if you want to design marker probes
    # if args.select_markers:
    #     total_transcript_table = pd.read_csv(args.custom_blast_db)
    #     phylogenetic_marker_genes = ['dnaG', 'frr', 'infC', 'nusA', 'pgk', 'pyrG', 'rplA', 'rplB', 'rplC', 'rplD', 'rplE', 'rplF', 'rplK', 'rplL', 'rplM', 'rplN',
    #                                 'rplP', 'rplS', 'rplT', 'rpmA', 'rpoB', 'rpsB', 'rpsC', 'rpsE', 'rpsI', 'rpsJ', 'rpsK', 'rpsM', 'rpsS', 'smpB', 'tsf']
    #     target_table_filename = get_target_table_filename(args.sample_dir)
    #     target_table = pd.read_csv(target_table_filename)
    #     input_dir = os.path.split(target_table_filename)
    #     for target_feature in target_table.target_feature:
    #         feature_marker_dir = '{}/{}'.format(input_dir, target_feature)
    #         if not os.path.exists(feature_marker_dir):
    #             os.mkdir(feature_marker_dir)
    #         feature_marker_fasta_filename = '{}/markers.fasta'.format(feature_marker_dir)
    #         marker_id_list_filename = '{}/marker_ids.csv'.format(feature_marker_dir)
    #         if not os.path.exists(feature_marker_fasta_filename) or not os.path.exists(marker_id_list_filename):
    #             target_genome_bin_classification = total_transcript_table[total_transcript_table['Gene'] == target_feature].Genome_bin_classification
    #             transcripts_target_genome_bin = total_transcript_table[total_transcript_table['Genome_bin_classification'] == target_genome_bin_classification]
    #             markers_target_genome_bin = transcripts_target_genome_bin[transcripts_target_genome_bin['Preferred_name'].isin(phylogenetic_marker_genes)]
    #             # Save marker id list for feature
    #             marker_id_list = markers_target_genome_bin['Gene']
    #             # feature_marker_dir = '{}/{}'.format(input_dir, target_feature)
    #             # if not os.path.exists(feature_marker_dir):
    #             #     os.mkdir(feature_marker_dir)
    #             # marker_id_list_filename = '{}/marker_ids.csv'.format(feature_marker_dir)
    #             marker_id_list.to_csv(marker_id_list_filename)
    #             # Write fasta for marker genes
    #             # feature_marker_fasta_filename = '{}/markers.fasta'.format(feature_marker_dir)
    #             feature_marker_fasta = []
    #             for index, marker in markers_target_genome_bin:
    #                 sequence = marker.Sequence
    #                 id = marker.Gene
    #                 record = SeqRecord(sequence, id = id)
    #                 feature_marker_fasta.append(record)
    #             with open(feature_marker_fasta_filename, 'w') as output_handle:
    #                 SeqIO.write(feature_marker_fasta, output_handle, 'fasta')
    #                 # marker_probe_design_outdir = '{}/{}_markers'.format(sim_primer3_dir, target_feature)
    #         marker_probe_design_outdir = '{}/{}_markers'.format(sim_primer3_dir, target_feature)
    #         if args.molecule == 'RNA':
    #             feature_marker_fasta_rc_filename = '{}_rc.fasta'.format(feature_marker_fasta_filename)
    #             if not os.path.exists(feature_marker_fasta_rc_filename):
    #                 generate_reverse_complement_fasta(feature_marker_fasta_filename, feature_marker_fasta_rc_filename)
    #         # marker_probe_design_outdir = '{}/{}_markers'.format(sim_primer3_dir, target_feature)
    #             probe_design(feature_marker_fasta_rc_filename, marker_probe_design_outdir, args.include_start, args.include_end)
    #             feature_marker_probes_filename = glob.glob('{}/*.int'.format(marker_probe_design_outdir))
    #             for f in feature_marker_probes_filename:
    #                 split_feature_probe_file(f)
    #         if args.molecule == 'DNA':
    #         # marker_probe_design_outdir = '{}/{}_markers'.format(sim_primer3_dir, target_feature)
    #             probe_design(feature_marker_fasta_filename, marker_probe_design_outdir, args.include_start, args.include_end)
    #             feature_marker_probes_filename = glob.glob('{}/*.int'.format(marker_probe_design_outdir))
    #             for f in feature_marker_probes_filename:
    #                 split_feature_probe_file(f)



    # probe_design(gene_sequences_filename, sim_primer3_dir, args.include_start, args.include_end)
    # feature_probes_filename = glob.glob('{}/*.int'.format(sim_primer3_dir))
    # for f in feature_probes_filename:
    #     split_feature_probe_file(f, args.molecule)
    # # if not os.path.exists(sim_dir + '/blast'):
    # #     os.makedirs(sim_dir + '/blast')
    if not os.path.exists(args.sample_dir + '/blast'):
        os.makedirs(args.sample_dir + '/blast')
    return

if __name__ == '__main__':
    main()
