import pandas as pd
from ete3 import NCBITaxa


################################################################################
# Functions
################################################################################

def get_lineage_at_desired_ranks(taxid, desired_ranks):
    'Retrieve lineage information at desired taxonomic ranks'
    # initiate an instance of the ncbi taxonomy database
    ncbi = NCBITaxa()

    # retrieve lineage information for each full length 16S molecule
    lineage = ncbi.get_lineage(taxid)
    lineage2ranks = ncbi.get_rank(lineage)
    ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
    ranki = [ranks2lineage.get(x) for x in desired_ranks]
    ranks = [x if x is not None else 0 for x in ranki]
    return(ranks)


################################################################################
# Variables
################################################################################

# Snakemake variables
fasta_filename = snakemake.input[0]
blast_output_filename = snakemake.input[1]
blast_lineage_filename = snakemake.output[0]
taxid_map_filename = snakemake.output[1]

# # Test variables
# fasta_filename = "test_data/synthetic_2019-07-18.fasta"
# blast_output_filename = "blast_output/synthetic_2019-07-18.blast.out"
# blast_lineage_filename = "blast_output/synthetic_2019-07-18_blast_lineage.csv"
# taxid_map_filename = "dbases/synthetic_2019-07-18_taxid_map.txt"


################################################################################
# Script
################################################################################

# read in blast result
blast_result = pd.read_csv(blast_output_filename, header=None, sep="\t")
blast_result.columns = ['molecule_id', 'reference_id', 'pid', 'qcovhsp', 'length',
                        'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send',
                        'evalue', 'bitscore', 'staxids']

# initiate an instance of ncbi taxonomy database
ncbi = NCBITaxa()

# retrieve lineage information for each full length 16S molecule
desired_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
ranks = pd.DataFrame(columns=['staxids'] + desired_ranks)
blast_result_staxids = blast_result['staxids'].unique()
ranks['staxids'] = blast_result_staxids
for i in range(0, blast_result_staxids.shape[0]):
    taxid = blast_result_staxids[i]
    if not str(taxid).isdigit():
        taxid = taxid.split(';')[0]
    ranks.ix[i, 1:len(desired_ranks) + 1] = get_lineage_at_desired_ranks(taxid, desired_ranks)

# merge lineage information with PacBio 16S blast results
blast_lineage = blast_result.merge(ranks, on='staxids', how='left')
taxid_map = blast_lineage[['molecule_id', 'species']]

# Write blast lineage and taxid_map to the output files as a csv and txt respectively
blast_lineage.to_csv(blast_lineage_filename, sep='\t', index=False)
taxid_map.to_csv(taxid_map_filename, sep=' ', index=None, header=None, mode='a')
