# consolidate info from different sources, imgt

import argparse
import simple_bio_seq as simple
import csv
import number_ighv

parser = argparse.ArgumentParser(description='Create a summary file that brings together info from the iglabel database with the sources of its sequences and IMGT names')
parser.add_argument('database', help='sequence database (csv)')
parser.add_argument('imgt_ref', help='imgt reference file (ungapped)')
parser.add_argument('imgt_v_gapped', help='imgt reference file (gapped V genes)')
parser.add_argument('subgroup', help='species subgroup')
parser.add_argument('-source', help='subdirectory/filename of source to add (multiple sources can be specified)', action='append')
parser.add_argument('output_file', help='output file (csv)')
args = parser.parse_args()

sources = {}

for source in args.source:
    genes = simple.read_fasta(source)
    source_name = source.split('/')[0]
    if source_name not in sources:
        sources[source_name] = {}
    if 'watson' in source_name:
        for name, seq in genes.items():
            if 'IG' not in name:
                breakpoint()
            type = 'IG' + name.split('IG')[1][:2]
            if seq not in sources[source_name]:
                sources[source_name][seq] = []
            sources[source_name][seq].append({'type': type, 'gene_name': name})
    if 'collins' in source_name:
        for name, seq in genes.items():
            type = 'IG' + name.split('IG')[1][:2]
            if '|' in name:
                name = name.split('|')[0]
            if seq not in sources[source_name]:
                sources[source_name][seq] = []
            sources[source_name][seq].append({'type': type, 'gene_name': name})


imgt_ref = simple.read_fasta(args.imgt_ref)
imgt_v_gapped = simple.read_fasta(args.imgt_v_gapped)
sources['imgt'] = {}

for name, seq in imgt_ref.items():
    if seq not in sources['imgt']:
        sources['imgt'][seq] = []
    type = 'IG' + name.split('IG')[1][:2]
    sources['imgt'][seq].append({'type': type, 'gene_name': name})

with open(args.database, 'r') as fi, open(args.output_file, 'w', newline='') as fo:
    reader = csv.DictReader(fi)
    headers = ['gene_label', 'type', 'functionality', 'inference_type', 'species_subgroup', 'subgroup_type']
    headers.extend(sources.keys())
    headers.extend(['alt_names', 'notes', 'affirmation', 'sequence', 'sequence_gapped', 'j_codon_frame', 'j_cdr3_end'])
    writer = csv.DictWriter(fo, fieldnames=headers)
    writer.writeheader()

    untracked_sequences = 0

    for row in reader:
        for seq in row['sequences'].split(','):
            source_found = False
            rec = {'gene_label': row['label'], 'type': None, 'functionality': 'F', 'inference_type': 'Rearranged',
                    'species_subgroup': args.subgroup, 'subgroup_type': 'strain', 'alt_names': '', 'notes': '', 'affirmation': '1', 'sequence': seq, 'sequence_gapped': seq}
            alt_names = []
            for source in sources.keys():
                rec[source] = ''
                if source == 'imgt':        # take a substring match
                    for imgt_seq in sources[source].keys():
                        if seq in imgt_seq:
                            rec[source] = ','.join([rec['gene_name'] for rec in sources[source][imgt_seq]])
                            rec['type'] = sources[source][imgt_seq][0]['type']
                            rec['gene_label'] = rec['type'] + '-' + row['label']
                else:
                    if seq in sources[source]:
                        rec[source] = ','.join([rec['gene_name'] for rec in sources[source][seq]])
                        rec['type'] = sources[source][seq][0]['type']
                        rec['gene_label'] = rec['type'] + '-' + row['label']
                        alt_names.extend([source + ':' + rec['gene_name'] for rec in sources[source][seq]])
                        source_found = True

            if source_found:
                if 'IGHV' in rec['type']:
                    res, aa, notes = number_ighv.gap_sequence(rec['sequence'], imgt_v_gapped, imgt_ref)
                    rec['sequence_gapped'] = res
                    rec['notes'] += notes
                rec['alt_names'] = ','.join(alt_names)
                writer.writerow(rec)
            else:
                untracked_sequences += 1

    if untracked_sequences:
        print(f'{untracked_sequences} records found in the database were not listed in the source and have not been added to the summary.')