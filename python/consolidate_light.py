# consolidate info from jackson et al for balb-c
from collections import namedtuple

from receptor_utils import simple_bio_seq as simple
import os

ldir = 'kos_et_al_2022_light_chain'


results = []

annotation_notes = {
    'IGKV4-52*01_S5019_C3H':'Conserved Trp not found',
    'IGKV9-129*01_C3H':'Conserved Trp not found',
    'IGKV13-84*01_S4500_PWD':'Conserved Trp not found',
    'IGKV4-92*01_S5510_PWD':'Second cysteine not found',
    'IGKV4-52*01_S5019_SJL':'Conserved Trp not found',
    'IGKV4-52*01_SJL':'Conserved Trp not found',
    'IGKV4-52*01_NZB':'Conserved Trp not found',
    'IGKV4-52*01_S5019_NZB':'Conserved Trp not found',
    'IGKV9-123*01_S4862_NZB':'Second cysteine not found',
    'IGKV4-52*01_S2975_B6':'Conserved Trp not found',
    'IGKV9-129*01_B6':'Conserved Trp not found',
    'IGKV5-45*01_S2588_MSM':'Conserved Trp not found',
    'IGKV13-84*01_S4500_MSM':'Conserved Trp not found',
    'IGKV4-52*01_S5019_129':'Conserved Trp not found',
    'IGKV9-129*01_AJ':'Conserved Trp not found',
    'IGKV4-52*01_AJ':'Conserved Trp not found',
    'IGKV4-52*01_S5019_AJ':'Conserved Trp not found',
    'IGKV9-129*01_DBA2':'Conserved Trp not found',
    'IGKV4-52*01_DBA2':'Conserved Trp not found',
    'IGKV4-52*01_S5019_DBA2':'Conserved Trp not found',
    'IGKV4-52*01_S2975_DBA1':'Conserved Trp not found',
    'IGKV9-129*01_DBA1':'Conserved Trp not found',
    'IGKV4-52*01_S5019_DBA1':'Conserved Trp not found',
    'IGKV9-129*01_AKR':'Conserved Trp not found',
    'IGKV4-52*01_BALB':'Conserved Trp not found',
    'IGKV9-129*01_BALB':'Conserved Trp not found',
    'IGKV10-96*06_CAST':'Conserved Trp not found',
    'IGLV5*01_CAST':'Stop codon in V-REGION',
    'IGLV6*01_CAST':'Stop codon in V-REGION',
    'IGLV4*01_CAST':'Stop codon in V-REGION',
    'IGLV8*01_CAST':'Stop codon in V-REGION',
    'IGLV7*01_CAST':'Stop codon in V-REGION',
    'IGLV5*01_MSM':'Stop codon in V-REGION',
    'IGLV6*03_MSM':'Stop codon in V-REGION',
}


def process_alt(name):
    Ret = namedtuple('Ret', 'func, notes')
    func = 'F'
    notes = ''

    if name in annotation_notes:
        if 'Stop' not in annotation_notes[name]:
            func = 'ORF'
        notes = annotation_notes[name]

    return Ret(func, notes)


def seq_name_lookup(seq, inf_recs, strain, chain):
    Ret = namedtuple('Ret', 'paper_name, imgt_name')
    imgt_name = ''
    paper_name = ''

    for rec in inf_recs:
        if seq == rec['Sequence'].upper() and rec['Strain_ID'] == strain and chain + 'V' == rec['Segment_Type']:
            if rec['Percent_Identity_to_IMGT'] == '100':
                imgt_name = rec['Closest_IMGT_Allele']
            paper_name = rec['Seq_Name']

    return Ret(paper_name, imgt_name)


def process_db(db_seqs, inf_recs, chain, strain, gapped_seqs, seq_lookup):
    results = []

    for rec in db_seqs:
        alt_names = []
        imgt_names = []
        rec_notes = []
        func = ''
        seqs = rec['sequences'].split(',')
        for seq in seqs:
            names = seq_name_lookup(seq, inf_recs, strain, chain)
            alt_names.append(names.paper_name)
            imgt_names.append(names.imgt_name)
            seq_notes = process_alt(names.paper_name)
            rec_notes.append(seq_notes.notes)
            func = seq_notes.func

        result = {
            'gene_label': chain + 'V-' + rec['label'],
            'type': chain + 'V',
            'functionality': func,
            'inference_type': 'Rearranged',
            'species_subgroup': strain,
            'subgroup_type': 'strain',
            'imgt': ','.join(imgt_names),
            'alt_names': 'Kos et al: ' + ','.join(alt_names),
            'notes': ','.join(rec_notes),
            'affirmation': '1',
            'sequence': rec['longest_seq'],
            'sequence_gapped': gapped_seqs[seq_lookup[rec['longest_seq']]],
        }

        results.append(result)

    return results


def main():
    for chain in ['IGK', 'IGL']:
        results = []
        inf_recs = simple.read_csv(os.path.join(ldir, f'{chain}V_Inferences.csv'))
        strain_map = {rec['Strain_ID']: rec['Strain_ID'].replace('/', '_').replace(' ', '') for rec in inf_recs}

        for orig_strain, strain in strain_map.items():
            dir = strain
            if 'C57BL' in strain:
                dir = 'B6'
            if 'BALB' in strain:
                dir = 'BALBC'
            if os.path.isdir(dir):
                print(f'{strain} : {dir}')
                db_seqs = simple.read_csv(os.path.join(dir, 'db', f'mouse_{dir}_{chain}.csv'))
                gapped_seqs = simple.read_fasta(os.path.join(dir, ldir, f'{chain}_V_gapped.fasta'))
                ungapped_seqs = simple.read_fasta(os.path.join(dir, ldir, f'{chain}_V.fasta'))
                seq_lookup = {v: k for k, v in ungapped_seqs.items()}
                results.extend(process_db(db_seqs, inf_recs, chain, orig_strain, gapped_seqs, seq_lookup))

        simple.write_csv(os.path.join(ldir, f'{chain}V_summary.csv'), results)


main()
