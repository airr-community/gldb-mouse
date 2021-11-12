# consolidate info from different sources, imgt

import argparse
import simple_bio_seq as simple
import csv
import number_ighv

parser = argparse.ArgumentParser(description='Find valid genes in a contig given blast matches')
parser.add_argument('database', help='output file (csv)')
parser.add_argument('imgt_ref', help='imgt reference file (ungapped)')
parser.add_argument('imgt_v_gapped', help='imgt reference file (gapped V genes)')
parser.add_argument('subgroup', help='species subgroup')
parser.add_argument('-source', help='subdirectory/filename of source to add)', action='append')
parser.add_argument('output_file', help='output file (csv)')
args = parser.parse_args()

sources = {}

for source in args.source:
    genes = simple.read_fasta(source)
    source_name = source.split('/')[0]
    sources[source_name] = {}
    if 'watson' in source_name:
        for name, seq in genes.items():
            type = 'IG' + name.split('IG')[1][:2]
            if seq not in sources[source_name]:
                sources[source_name][seq] = []
            sources[source_name][seq].append({'type': type, 'gene_name': name})

imgt_ref = simple.read_fasta(args.imgt_ref)
imgt_v_gapped = simple.read_fasta(args.imgt_v_gapped)
sources['imgt'] = {}

for name, seq in imgt_ref.items():
    if seq not in sources['imgt']:
        sources['imgt'][seq] = []
    sources['imgt'][seq].append({'type': type, 'gene_name': name})

with open(args.database, 'r') as fi, open(args.output_file, 'w', newline='') as fo:
    reader = csv.DictReader(fi)
    headers = ['gene_label', 'type', 'functional', 'inference_type', 'species_subgroup']
    headers.extend(sources.keys())
    headers.extend(['alt_names', 'notes', 'sequence', 'sequence_gapped'])
    writer = csv.DictWriter(fo, fieldnames=headers)
    writer.writeheader()

    for row in reader:
        for seq in row['sequences'].split(','):
            rec = {'gene_label': row['label'], 'type': None, 'functional': 'Y', 'inference_type': 'Rearranged', 
                    'species_subgroup': args.subgroup, 'alt_names': '', 'notes': '', 'sequence': seq, 'sequence_gapped': ''}
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
            if 'IGHV' in rec['type']:
                res, aa, notes = number_ighv.gap_sequence(rec['sequence'], imgt_v_gapped, imgt_ref)
                rec['sequence_gapped'] = res
                rec['notes'] += notes
            rec['alt_names'] = ','.join(alt_names)
            

            writer.writerow(rec)