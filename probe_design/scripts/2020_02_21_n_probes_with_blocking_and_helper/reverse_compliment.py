# Generate the reverse compliment form a genomic fasta file, then add sense and antisense strands to a new fasta file
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import re


def generate_reverse_complement_fasta(fasta_filename):
    fasta = SeqIO.parse(fasta_filename, "fasta")
    reverse = []
    for record in fasta:
        # print(re.sub('\|','_',record.description))
        # print(record.id)
        rc = SeqRecord(record.seq.reverse_complement())
        rc.id = re.sub('\|','_rc_',record.id)
        rc.description = re.sub('\|','_rc_',record.description)
        reverse.append(rc)
    # combined_seqs
    dir, name = os.path.split(fasta_filename)
    name = re.sub('.fasta|.fna','',name)
    # name = re.sub('.fna','',name)
    reverse_filename = '{}/{}_reverse_complement.fasta'.format(dir, name)
    # combined_filename
    with open(reverse_filename, 'w') as output_handle:
        SeqIO.write(reverse, output_handle, 'fasta')

coding_fasta_filename = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/databases/e_coli_mg1655_coding/e_coli_coding_sequences.fasta'
os.path.exists(coding_fasta_filename)
generate_reverse_complement_fasta(coding_fasta_filename)

reverse_filename = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/databases/e_coli_mg1655_coding/e_coli_coding_sequences_reverse_complement.fasta'
dict = SeqIO.to_dict(SeqIO.parse(reverse_filename, "fasta"))
for record in dict:
    print(record)
# def main():
#     parser = argparse.ArgumentParser('Generate the reverse compliment form a genomic fasta file, then add sense and antisense strands to a new fasta file')
#     # Target sequences
#     parser.add_argument('gene_sequences_filename', type = str, help = 'Input genomic FASTA file')
#     # parser.add_argument('snapgene_filename', type = str, help = 'Input FASTA file containing 16S sequences')
#     # parser.add_argument('custom_blast_db', type = str, help = 'Blast database')
#     # Csv containing information from hao Zhou
#     parser.add_argument('custom_blast_db', type = str, help = 'Blast database')
#     # parser.add_argument('blast_db_csv_filename', type = str, help = 'Blast database')
#     parser.add_argument('-l', '--length_threshold', dest = 'length_threshold', type = int, default = 50, help = 'Length threshold')
#     parser.add_argument('-b', '--include_start', dest = 'include_start', type = int, default = 1, help = 'Starting position of included region for probe design, measured from the beginning of 16S rRNA molecules')
#     parser.add_argument('-e', '--include_end', dest = 'include_end', type = int, default = 1, help = 'Ending position of included region for probe design, measured from the end of 16S rRNA molecules')
#     args = parser.parse_args()

def combine_sense_antisense_fasta(genomic_fasta_filename):
    genomic_fasta = SeqIO.parse(genomic_fasta_filename, "fasta")
    combined_seqs = []
    for record in genomic_fasta:
        rc = SeqRecord(record.seq.reverse_complement(), id = str(record.id + '_reverse_complement'), description = record.description)
        combined_seqs.append(record)
        combined_seqs.append(rc)
    # combined_seqs
    dir, name = os.path.split(genomic_fasta_filename)
    database_dir, database_dir_named = os.path.split(dir)
    database_dir_named = database_dir + '/' + database_dir_named + '_sense_antisense_combined'
    name = re.sub('.fasta|.fna','',name)
    # name = re.sub('.fna','',name)
    combined_filename = '{}/{}_sense_antisense_combined.fasta'.format(database_dir_named, name)
    if ~os.path.exists(database_dir_named):
        os.mkdir(database_dir_named)
    # combined_filename
    with open(combined_filename, 'w') as output_handle:
        SeqIO.write(combined_seqs, output_handle, 'fasta')
    # SeqIO.write(combined_seqs, combined_filename, 'fasta')

genomic_fasta_filename = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/databases/e_coli_mg1655_genomic/GCF_000005845.2_ASM584v2_genomic.fna'

combine_sense_antisense_fasta(genomic_fasta_filename)




def count_seqs(fasta_filename):
    dict = SeqIO.to_dict(SeqIO.parse(fasta_filename, "fasta"))
    return(len(dict))

def make_dict(fasta_filename):
    dict = SeqIO.to_dict(SeqIO.parse(fasta_filename, "fasta"))
    return(dict)

# cds_from_genomic_filename = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/databases/GCF_000005845.2_ASM584v2_cds_from_genomic.fna'
count_seqs(cds_from_genomic_filename)
count_seqs(coding_fasta_filename)

# Get dict ids from a fasta file, input a list of genes of interest
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


def sum_seq_lengths(fasta_filename):
    fasta = SeqIO.parse(fasta_filename, "fasta")
    count = 0
    for record in fasta:
        count += len(record.seq)
    return(count)

sum_seq_lengths(cds_from_genomic_filename)
sum_seq_lengths(coding_fasta_filename)
sum_seq_lengths(genomic_fasta_filename)

# Replace characters in fasta id such as "|" make sure to escape, i.e. "\|"
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

cds_from_genomic_filename = '/workdir/bmg224/hiprfish/mobile_elements/probe_design/databases/e_coli_mg1655_coding/GCF_000005845.2_ASM584v2_cds_from_genomic.fna'
c_glutamicum_database_filename = '/workdir/bmg224/hiprfish/transcript_label/probe_design/database/c_glutamicum_coding/c_glutamicum_coding_sequences.fasta'
e_coli_database_filename = '/workdir/bmg224/hiprfish/transcript_label/probe_design/database/e_coli_mg1655_coding/e_coli_coding_sequences.fasta'
ec_gene_ids = get_specific_gene_id(e_coli_database_filename, ['pck','Cgl1951'])
cg_gene_ids = get_specific_gene_id(c_glutamicum_database_filename, ['Cgl1951','pck'])
