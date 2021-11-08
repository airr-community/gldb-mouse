# consolidate info from different sources, imgt

import argparse
import simple_bio_seq as simple
import csv

parser = argparse.ArgumentParser(description='Find valid genes in a contig given blast matches')
parser.add_argument('database', help='output file (csv)')
parser.add_argument('imgt_ref', help='output file (csv)')
parser.add_argument('-source', help='subdirectory/filename of source to add)', action='append')
parser.add_argument('output_file', help='output file (csv)')
args = parser.parse_args()

sources = {}

def add_name_and_type(rec, label, type):
    if rec['type'] is None:
        rec['type'] = type
        rec['gene_label'] = 'IGH' + type[0] + '-' + label

    return


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
sources['imgt'] = {}

for name, seq in imgt_ref.items():
    if seq not in sources['imgt']:
        sources['imgt'][seq] = []
    sources['imgt'][seq].append({'type': type, 'gene_name': name})

with open(args.database, 'r') as fi, open(args.output_file, 'w', newline='') as fo:
    reader = csv.DictReader(fi)
    headers = ['gene_label', 'type']
    headers.extend(sources.keys())
    headers.append('sequence')
    writer = csv.DictWriter(fo, fieldnames=headers)
    writer.writeheader()

    for row in reader:
        for seq in row['sequences'].split(','):
            rec = {'gene_label': row['label'], 'type': None, 'sequence': seq}
            for source in sources.keys():
                rec[source] = ''
                if source == 'imgt':        # take a substring match
                    for imgt_seq in sources[source].keys():
                        if seq in imgt_seq:
                            rec[source] = ','.join([rec['gene_name'] for rec in sources[source][imgt_seq]])
                            rec['type'] = sources[source][imgt_seq][0]['type']
                else:
                    if seq in sources[source]:
                        rec[source] = ','.join([rec['gene_name'] for rec in sources[source][seq]])
                        rec['type'] = sources[source][seq][0]['type']
            writer.writerow(rec)