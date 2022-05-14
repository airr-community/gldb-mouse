# consolidate light j chain info from kos et al

from receptor_utils import simple_bio_seq as simple





for chain in ['IGK', 'IGL']:
    results = []
    db_seqs = simple.read_csv(f'J_genes/db/mouse_J_{chain}.csv')
    ungapped_seqs = simple.read_fasta(f'J_genes/kos_et_al_2022_light_chain/{chain}J.fasta')
    seq_name_lookup = {v: k for k, v in ungapped_seqs.items()}

    for rec in db_seqs:
        alt_names = []
        seqs = rec['sequences'].split(',')
        for seq in seqs:
            alt_names.append(seq_name_lookup[seq])

        if '_' not in alt_names[0]:
            imgt_name = alt_names[0]
        else:
            imgt_name = ''


        result = {
            'gene_label': f"{chain}J-{rec['label']}",
            'type': f'{chain}J',
            'functionality': 'F',
            'inference_type': 'Rearranged',
            'species_subgroup': '',
            'subgroup_type': '',
            'imgt': imgt_name,
            'alt_names': 'Kos et al: ' + ','.join(alt_names),
            'notes': '',
            'affirmation': '1',
            'sequence': rec['longest_seq'],
            'sequence_gapped': '.',
        }

        results.append(result)


    simple.write_csv(f'J_genes/{chain}_summary.csv', results)
