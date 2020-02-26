import pandas as pd
from Bio import SeqIO
import re

################################################################################
# Functions
################################################################################


def retrieve_cluster(cluster_lookup, taxon, taxon_fasta_filename,
                     consensus_fasta_filename, consensus_uc_filename,
                     fasta_num_seq, consensus_num_seq):
    # If there's only one sequence for the taxon, assign the taxon as the strain
    if fasta_num_seq == 1:
        fasta_record = SeqIO.read(taxon_fasta_filename, 'fasta')
        molecule_cluster = pd.DataFrame([fasta_record.id, str(taxon)]).transpose()
        molecule_cluster.columns = ['molecule_id', 'strain']
        cluster_lookup = pd.concat([cluster_lookup, molecule_cluster], ignore_index=True)

    # If clustering results in only one sequence, assign the taxon as the strain
    elif consensus_num_seq == 1:
        molecule_cluster = pd.DataFrame([[record.id, taxon]
                                        for record in SeqIO.parse(taxon_fasta_filename,
                                                                  'fasta')])
        molecule_cluster.columns = ['molecule_id', 'strain']
        cluster_lookup = pd.concat([cluster_lookup, molecule_cluster], ignore_index=True)

    # If clustering results in substrains, create strain names and assign to each molecule id
    else:
        taxon_uc = pd.read_csv(consensus_uc_filename, sep='\t', header=None)
        taxon_uc.columns = ['rec_type', 'cluster_num', 'seq_length', 'pid', 'strand',
                            'notused', 'notused', 'comp_alignment', 'query_label', 'target_label']
        # Get the appropriate row from each cluster ("or" operator?)
        taxon_uc_clusters = taxon_uc[(taxon_uc['rec_type'] == 'S') | (taxon_uc['rec_type'] == 'H')]
        molecule_cluster = pd.DataFrame([[taxon_uc_clusters.loc[i, 'query_label'],
                                         str(taxon) + '_Cluster_' + str(taxon_uc_clusters.loc[i,
                                                                        'cluster_num'])]
                                         for i in range(0, taxon_uc_clusters.shape[0])])
        molecule_cluster.columns = ['molecule_id', 'strain']
        cluster_lookup = pd.concat([cluster_lookup, molecule_cluster], ignore_index=True)
    return(cluster_lookup)


################################################################################
# Variables
################################################################################
# Snakemake variables
blast_lineage_filename = snakemake.input[0]
consensus_fasta_list = snakemake.input[1:]
blast_lineage_strain_filename = snakemake.output[0]
taxon_abundance_filename = snakemake.output[1]
target_rank = snakemake.params[0]
otu = snakemake.params[1]

# # Test variables
# blast_lineage_filename = "blast_output/synthetic_2019-07-18/blast_lineage.csv"
# taxon_fasta_list = ["taxon_consensus/synthetic_2019-07-18/1280.grouped.fasta",
#                     "taxon_consensus/synthetic_2019-07-18/1747.grouped.fasta",
#                     "taxon_consensus/synthetic_2019-07-18/562.grouped.fasta",
#                     "taxon_consensus/synthetic_2019-07-18/1282.grouped.fasta"]
# blast_lineage_strain_filename = "blast_output/synthetic_2019-07-18/blast_lineage_strain.csv"
# taxon_abundance_filename = "readouts/synthetic_2019-07-18/taxon_abundance.csv"
# target_rank = 'species'
# otu = 'F'

################################################################################
# Script
################################################################################
# Read in blast lineage dataframe
blast_lineage = pd.read_csv(blast_lineage_filename, sep='\t')

# Create output DataFrame
cluster_lookup = pd.DataFrame(columns=['molecule_id', 'strain'])

# Iterate through input filenames
for consensus_fasta_filename in consensus_fasta_list:
    # Get the consensus filename and summary filename for the taxon
    taxon_fasta_filename = re.sub(r'\w+(?=.fasta)', 'grouped', consensus_fasta_filename)
    consensus_uc_filename = re.sub('.consensus.fasta', '.uc', consensus_fasta_filename)

    # Determine the taxon
    taxon = re.search(r'\w+(?=.consensus.fasta)', consensus_fasta_filename)
    taxon = taxon.group(0)

    # Count the number of sequences in the taxon fasta and consensus file
    fasta_num_seq = sum(1 for record in SeqIO.parse(taxon_fasta_filename, 'fasta'))
    consensus_num_seq = sum(1 for record in SeqIO.parse(consensus_fasta_filename, 'fasta'))

    # Retrieve cluster information and assign a strain to each molecule id
    cluster_lookup = retrieve_cluster(cluster_lookup, taxon, taxon_fasta_filename,
                                      consensus_fasta_filename, consensus_uc_filename,
                                      fasta_num_seq, consensus_num_seq)

# Merge strain assignments with blast lineage DataFrame and save as csv
blast_lineage_strain = blast_lineage.merge(cluster_lookup, on='molecule_id')
blast_lineage_strain.to_csv(blast_lineage_strain_filename, sep='\t', index=False)

# Get abundance for each taxon and save to a csv file
taxon_abundance = blast_lineage_strain[target_rank].value_counts().reset_index()
taxon_abundance.columns = ['taxid', 'counts']
taxon_abundance['portion_of_whole'] = round(taxon_abundance.counts / taxon_abundance.counts.sum(),
                                            5)
taxon_abundance.to_csv(taxon_abundance_filename, sep='\t', index=False)
