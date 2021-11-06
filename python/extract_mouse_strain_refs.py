# Extract reference files for nominated mouse strains

import argparse

from Bio import SeqIO
import simple_bio_seq as simple

imgt_file = 'imgt/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP'
strains = ['musculus', 'castaneus', 'domesticus']
gene_types = ['IGHV', 'IGHD', 'IGHJ']


strain_seqs = {}

for strain in strains:
    strain_seqs[strain] = {}

recs = SeqIO.parse(imgt_file, 'fasta')

for rec in recs:
    for gene_type in gene_types:
        imgt_species = rec.description.split('|')[2]
        if '-REGION' in rec.description and gene_type in rec.description and 'Mus musculus' in imgt_species:
            name = rec.description.split('|')[1]
            if imgt_species == 'Mus musculus':
                if 'musculus' in strains:
                    strain_seqs['musculus'][name] = str(rec.seq).upper()
            else:
                for strain in strains:
                    if strain in imgt_species:
                        strain_seqs[strain][name] = str(rec.seq).upper()
            
for strain in strains:
    simple.write_fasta(strain_seqs[strain], 'imgt/%s_imgt_ungapped.fasta' % strain)

