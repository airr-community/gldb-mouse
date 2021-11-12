# Extract imgt reference files for nominated mouse strains

# Put everything in to musculus, and other strains into their own files if listed

import argparse

from Bio import SeqIO
import simple_bio_seq as simple

imgt_file = 'imgt/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP'
strains = ['musculus']
gene_types = ['IGHV', 'IGHD', 'IGHJ']


strain_seqs = {}

for strain in strains:
    strain_seqs[strain] = {}
    strain_seqs[strain]['gapped'] = {}
    strain_seqs[strain]['ungapped'] = {}

recs = SeqIO.parse(imgt_file, 'fasta')

for rec in recs:
    for gene_type in gene_types:
        imgt_species = rec.description.split('|')[2]
        if '-REGION' in rec.description and gene_type in rec.description and 'Mus musculus' in imgt_species:
            name = rec.description.split('|')[1]
            if 'Mus musculus' in imgt_species:
                if 'musculus' in strains:
                    if 'V' in gene_type:
                        strain_seqs['musculus']['gapped'][name + '|' + imgt_species.replace(' ', '_')] = str(rec.seq).upper()

                    strain_seqs['musculus']['ungapped'][name + '|' + imgt_species.replace(' ', '_')] = str(rec.seq).upper().replace('.', '')

                for strain in strains:
                    if strain != 'musculus' and strain in imgt_species:
                        if 'V' in gene_type:
                            strain_seqs[strain]['gapped'][name + '|' + imgt_species] = str(rec.seq).upper()

                        strain_seqs[strain]['ungapped'][name + '|' + imgt_species] = str(rec.seq).upper().replace('.', '')
            
for strain in strains:
    simple.write_fasta(strain_seqs[strain]['gapped'], 'imgt/%s_imgt_V_gapped.fasta' % strain)
    simple.write_fasta(strain_seqs[strain]['ungapped'], 'imgt/%s_imgt_ungapped.fasta' % strain)

