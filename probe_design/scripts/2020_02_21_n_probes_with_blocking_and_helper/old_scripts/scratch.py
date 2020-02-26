

def consensus_fasta_input(wildcards):
    checkpoint_output = checkpoints.write_taxon_fasta.get(**wildcards).output[0]
    return expand("taxon_consensus/{sample}/{taxon}.consensus.fasta",
                  sample=wildcards.sample,
                  taxon=glob_wildcards(os.path.join(checkpoint_output,
                                                    "{taxon}.grouped.fasta")).taxon)
