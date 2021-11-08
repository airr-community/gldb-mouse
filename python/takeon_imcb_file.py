# take on multi-strain data from imcb csv file, creating subdirectories and splitting files

import csv
import os.path
import os
import simple_bio_seq as simple

imcb_file = 'imcb12288-sup-0006-tables2.csv'
strain_seqs = {}

# read strains from file

with open(imcb_file, 'r') as fi:
    reader = csv.DictReader(fi)
    for row in reader:
        strain, lab = row['Strain_ID'].split('/')
        strain = '%s_%s' % (strain, lab)
        name = '%s_%s' % (strain, row['Sequence_ID'].split('_')[1])
        if strain not in strain_seqs:
            strain_seqs[strain] = {}
        strain_seqs[strain][name] = row['Sequence']


# establish directories for each strain

for strain in strain_seqs.keys():
    if not os.path.isdir(strain):
        os.mkdir(strain)

# create fasta for each strain

for strain, seqs in strain_seqs.items():
    simple.write_fasta(seqs, os.path.join(strain, 'imcb_%s.fasta' % strain))