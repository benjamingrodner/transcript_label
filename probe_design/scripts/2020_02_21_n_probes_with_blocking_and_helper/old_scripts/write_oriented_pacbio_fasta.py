# write_oriented_pacbio_fasta
# hiprfish probe design script
# Ben Grodner and Hao Shi
# Edited 2019_07_12

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np

# Snakemake variables
fasta_filename = snakemake.input[0]
blast_output_filename = snakemake.input[1]
oriented_fasta_filename = snakemake.output[0]
# # Test variables
# fasta_filename = 'data/ecoli16s.fasta'
# blast_output_filename = 'blast_output/ecoli16s.blast.out'
# oriented_fasta_filename = 'data/ecoli16s.oriented.fasta'

# Create a dictionary from the input fasta file
seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_filename, 'fasta'))
# Read the blast results into a pandas object
blast_result = pd.read_csv(blast_output_filename, sep='\t', header=None)
# Column names were specified in the blastn command
blast_result.columns = ['molecule_id', 'reference_id', 'pid', 'qcovhsp',
                        'length', 'mismatch', 'gapopen', 'qstart', 'qend',
                        'sstart', 'send', 'evalue', 'bitscore', 'staxids']
# Create a column labeling already oriented sequences
blast_result['oriented'] = (blast_result['sstart'] < blast_result['send'])
# Select oriented sequences, reverse complement other sequences, and add the
# sequence to a list
oriented_seqs = []
midz = blast_result['molecule_id']
for mid in midz.unique():
    bool_oriented = blast_result[midz == mid]['oriented'].values
    frac = sum(bool_oriented) / len(bool_oriented)
    if np.all(bool_oriented) or frac >= 0.5:
        oriented_seqs.extend([seq_dict[mid]])
    elif not np.all(bool_oriented) or frac < 0.5:
        oriented_seqs.extend(
            [SeqRecord(seq_dict[mid].reverse_complement().seq,
                       id=mid, description='')])
# Write the oriented list to a file
SeqIO.write(oriented_seqs, oriented_fasta_filename, 'fasta')
