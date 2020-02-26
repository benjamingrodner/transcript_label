from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
import glob
import re
import argparse
import random
import os

###############################################################################################################
###############################################################################################################
# Workflow functions
###############################################################################################################
###############################################################################################################


################################################################################
### Randomly reduce the number of taxa in the PacBio dictionary and output a new
### list of seqences
################################################################################
def reduce_taxon_fasta_random(all_taxids, reduction, blast_lineage, seq_dict, target_rank):
    ## How many taxids to remove
    n_reduction = np.round(all_taxids.shape[0]*reduction).astype(int)
    ## Random reduced list of taxonomic ids
    reduced_taxids = random.sample(list(all_taxids), n_reduction)
    ## Group molecular ids by taxonimic ids and filter using the reduced list
    blast_lineage_reduced = blast_lineage.groupby([target_rank]).filter(lambda x: x[target_rank].unique()[0] in reduced_taxids)
    ## New reduced library
    seq_dict_reduced = [seq_dict[mid] for mid in blast_lineage_reduced['molecule_id']]
    return(seq_dict_reduced, blast_lineage_reduced)


################################################################################
### Keep some percentage of the most abundant taxa
################################################################################
def reduce_taxon_fasta_top(all_taxids, reduction, blast_lineage, seq_dict):
    ## Order taxa by abundance
    ## Make a list of reduced taxids

    ## Group molecular ids by taxonimic ids and filter using the reduced list
    blast_lineage_reduced = blast_lineage.groupby([target_rank]).filter(lambda x: x[target_rank].unique()[0] in reduced_taxids)
    ## New reduced library
    seq_dict_reduced = [seq_dict[mid] for mid in blast_lineage_reduced['molecule_id']]
    return(seq_dict_reduced)


################################################################################
### Keep some percentage of the least abundant taxa
################################################################################
def reduce_taxon_fasta_bottom(all_taxids, reduction, blast_lineage, seq_dict):
    ## Order taxa by abundance
    ## Make a list of reduced taxids

    ## Group molecular ids by taxonimic ids and filter using the reduced list
    blast_lineage_reduced = blast_lineage.groupby([target_rank]).filter(lambda x: x[target_rank].unique()[0] in reduced_taxids)
    ## New reduced library
    seq_dict_reduced = [seq_dict[mid] for mid in blast_lineage_reduced['molecule_id']]
    return(seq_dict_reduced, blast_lineage_reduced)


################################################################################
### Create appropriate folders and save appropriate files
################################################################################
def save_reduced_things(seq_dict_reduced, blast_lineage_reduced, sample, method, r, data_dir):
    ## Save library as fasta in new sample name
    sample_reduced = sample + '_' + method + '_' + str(r) + '_reduction'
    sample_reduced_directory = data_dir + sample_reduced + '/'
    if not os.path.exists(sample_reduced_directory):
        os.makedirs(sample_reduced_directory)
    reduced_pacbio_fasta_directory = sample_reduced_directory + 'input/'
    if not os.path.exists(reduced_pacbio_fasta_directory):
        os.makedirs(reduced_pacbio_fasta_directory)
    reduced_oriented_pacbio_fasta_filename = reduced_pacbio_fasta_directory + sample_reduced + '.oriented.fasta'
    SeqIO.write(seq_dict_reduced, reduced_oriented_pacbio_fasta_filename, 'fasta')

    ## Save blast lineage as table in new utilities folder
    blast_lineage_directory = sample_reduced_directory + 'utilities/'
    if not os.path.exists(blast_lineage_directory):
        os.makedirs(blast_lineage_directory)
    blast_lineage_filename = blast_lineage_directory + 'blast_lineage_' + method + '_' + str(r) + '_reduction.csv'
    blast_lineage_reduced.to_csv(blast_lineage_filename, sep = '\t', index = False)
    return


###############################################################################################################
###############################################################################################################
# main function
###############################################################################################################
###############################################################################################################


def main():
    parser = argparse.ArgumentParser('Reduce species in PacBio sequenceing')
    # input 16S sequence
    parser.add_argument('input_fasta_filename', type = str, help = 'Input FASTA file containing 16S sequences')
    parser.add_argument('-t', '--target_rank', dest = 'target_rank', type = str, default = '.', help = 'The taxonomic level at which to design FISH probes')
    parser.add_argument('-s', '--similarity', dest = 'similarity', type = float, default = 90, help = 'Similarity threshold for grouping pacbio 16S sequences to produce consensus sequences')
    parser.add_argument('-m', '--method', dest = 'method', type = str, default = 'random', help = 'How should the community be reduced? Options: random, top, bottom')
    parser.add_argument('-r', '--reduction', dest = 'reduction', type = str, default = [0.5], help = 'Method random: What percent of original sample is retained (decimal); Method top: What percent of most common taxa to retain; Method bottom: What percent of least common taxa to retain')
    args = parser.parse_args()

    ## The goal is to reduce the number of species randomly. You also want to do
    ## some runs where you remove the least common species (looking at common
    ## species) as well as vice versa (looking at uncommon species).
    ## You will have to identify species alignment in the PacBio, then generate
    ## a new oriented blast database with reduced species. You also might have
    ## to create a new taxon consensus fasta.

    seq_dict = SeqIO.to_dict(SeqIO.parse(args.input_fasta_filename, 'fasta'))
    data_dir = re.sub('\w+/\w+/\w+.oriented.fasta','',args.input_fasta_filename)
    sample = re.search('\w+(?=.oriented.fasta)',args.input_fasta_filename)
    sample = sample.group()
    sample_dir = data_dir + '/' + sample
    reduction = [float(args.reduction)]
    method = args.method
    target_rank = args.target_rank
    similarity = args.similarity
    sim_dir = sample_dir + '/'+  target_rank + '/' + 's_' + str(similarity)
    taxon_consensus_output_directory = sim_dir + '/consensus'
    outdir = taxon_consensus_output_directory
    blast_lineage_filename = taxon_consensus_output_directory + '/blast_lineage.tab'
    blast_lineage = pd.read_csv(blast_lineage_filename, sep = '\t')
    ## Get all unique taxonomic ids
    all_taxids = blast_lineage[target_rank].unique()

    ## Reduce taxa depending on method
    if args.method == 'random':
        for r in reduction:
            ## New reduced library
            seq_dict_reduced, blast_lineage_reduced = reduce_taxon_fasta_random(all_taxids, r, blast_lineage, seq_dict, target_rank)
            save_reduced_things(seq_dict_reduced, blast_lineage_reduced, sample, method, r, data_dir)
    elif args.method == 'top':
        for r in reduction:
            seq_dict_reduced = reduce_taxon_fasta_top(all_taxids, r, blast_lineage, seq_dict)
            save_reduced_things(seq_dict_reduced, blast_lineage_reduced, sample, method, r, data_dir)
    return

if __name__ == '__main__':
    main()
