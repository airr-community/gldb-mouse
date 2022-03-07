# consolidate info from jackson et al for balb-c

from receptor_utils import simple_bio_seq as simple

db_seqs = simple.read_csv('BALBC/db/mouse_balbc.csv')
ungapped_seqs = simple.read_fasta('BALBC/jackson_et_al/balbc_ighv_ungapped.fasta')
gapped_seqs = simple.read_fasta('BALBC/jackson_et_al/balbc_ighv_gapped.fasta')
seq_name_lookup = {v: k for k, v in ungapped_seqs.items()}


results = []

non_functional = {
    'balbIGHV041': 'ORF',
    'IGHV1S10*01': 'P',
    'IGHV2-2-1*01': 'P',
    'IGHV14-2*02': 'P',
}

def process_alt(name):
    alt_name = ''
    imgt_name = ''
    notes = ''

    if name[:3] == 'IGH':
        imgt_name = name

    elif name[0] == 'J':
        alt_name = f'Haines: {name}'
        notes = f'Reported by Haines et al. (2001) as {name}'

    elif 'mus' in name:
        alt_name = f'VBASE2: {name}'
        notes = f'Reported in VBASE2 as {name}. Not previously reported in Haines et al. (2001)'

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
        'species_subgroup': 'BALB/c',
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

    result['alt_names'] = ','.join(result['alt_names'])
    result['imgt'] = ','.join(result['imgt'])
    result['notes'] = '\r'.join(result['notes'])
    results.append(result)


simple.write_csv('BALBC/summary.csv', results)
