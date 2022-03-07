# consolidate info from jackson et al for b6

from receptor_utils import simple_bio_seq as simple

db_seqs = simple.read_csv('B6/db/mouse_b6.csv')
ungapped_seqs = simple.read_fasta('B6/jackson_et_al/b6_ighv_ungapped.fasta')
gapped_seqs = simple.read_fasta('B6/jackson_et_al/b6_ighv_gapped.fasta')
seq_name_lookup = {v: k for k, v in ungapped_seqs.items()}


results = []

non_functional = {
    'IGHV1-13*01': 'P',
    'IGHV1-20*01': 'P',
    'IGHV1-21-1*01': 'P',
    'IGHV1-83*01': 'P',
    'IGHV4-2*01': 'P',
    'IGHV5-21*01': 'P',
    'IGHV8-7*01': 'P',
    'IGHV1-23*01': 'ORF',
    'IGHV1-62-3*01': 'ORF',
    'b6IGHV040': 'ORF',
}

alt_imgt = {
    'IGHV2-9-1*01': 'IGHV2-9*02',
}

def process_alt(name):
    alt_name = ''
    imgt_name = ''
    notes = ''

    if name[:3] == 'IGH':
        imgt_name = name

    else:
        alt_name = f'JACKSON: {name}'
        notes = f'Reported in Jackson et al. (2022) as {name}'

    return {'alt_name': alt_name, 'imgt_name': imgt_name, 'notes': notes}


for rec in db_seqs:
    alt_names = []
    seqs = rec['sequences'].split(',')
    for seq in seqs:
        alt_names.append(seq_name_lookup[seq])

    result = {
        'gene_label': 'IGHV-' + rec['label'],
        'type': 'IGHV',
        'functionality': 'F' if seq_name_lookup[rec['longest_seq']] not in non_functional else non_functional[seq_name_lookup[rec['longest_seq']]],
        'inference_type': 'Rearranged',
        'species_subgroup': 'C57BL/6',
        'subgroup_type': 'strain',
        'imgt': [],
        'alt_names': [],
        'notes': [],
        'affirmation': '1',
        'sequence': rec['longest_seq'],
        'sequence_gapped': gapped_seqs[seq_name_lookup[rec['longest_seq']]],
    }

    for alt_name in alt_names:
        alt_fields = process_alt(alt_name)
        result['alt_names'].append(alt_fields['alt_name'])
        result['imgt'].append(alt_fields['imgt_name'])
        result['notes'].append(alt_fields['notes'])

    if alt_fields['imgt_name'] in alt_imgt:
        result['imgt'].append(alt_imgt[alt_fields['imgt_name']])

    result['alt_names'] = ','.join(result['alt_names'])
    result['imgt'] = ','.join(result['imgt'])
    result['notes'] = '\r'.join(result['notes'])
    results.append(result)


simple.write_csv('B6/summary.csv', results)
