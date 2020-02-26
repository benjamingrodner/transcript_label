import pandas as pd
from Bio import SeqIO


################################################################################
# Functions
################################################################################

def write_taxon_fasta(df, taxon, seq_dict, outdir):
    # specify file names
    if taxon == 'strain':
        taxon = 'species'
    taxon_name = str(df[taxon].unique()[0])
    taxon_fasta_filename = outdir + '/' + taxon_name + '.grouped.fasta'
    # get sequences belonging to a particular taxon
    taxon_seqs = [seq_dict[mid] for mid in df['molecule_id']]
    # write taxon specific sequences to file
    SeqIO.write(taxon_seqs, taxon_fasta_filename, 'fasta')


################################################################################
# Variables
################################################################################
# Snakemake variables
oriented_fasta_filename = snakemake.input[0]
blast_lineage_filename = snakemake.input[1]
taxon_consensus_directory = snakemake.output[0]
target_rank = snakemake.params[0]

# # Test variables
# oriented_fasta_filename = "dbases/synthetic_2019-07-18/synthetic_2019-07-18.oriented.fasta"
# blast_lineage_filename = "blast_output/synthetic_2019-07-18_blast_lineage.csv"
# taxon_consensus_directory = "taxon_consensus/"
# target_rank = "species"

################################################################################
# Script
################################################################################

# Get blast lineage
blast_lineage = pd.read_csv(blast_lineage_filename, sep='\t')

# Get sample as a dictionary
seq_dict = SeqIO.to_dict(SeqIO.parse(oriented_fasta_filename, 'fasta'))

# Write sequences grouped by target rank
if target_rank == 'strain':
    blast_lineage.groupby(['species']).apply(write_taxon_fasta,
                                             taxon=target_rank, seq_dict=seq_dict,
                                             outdir=taxon_consensus_directory)
else:
    blast_lineage.groupby([target_rank]).apply(write_taxon_fasta,
                                               taxon=target_rank, seq_dict=seq_dict,
                                               outdir=taxon_consensus_directory)
