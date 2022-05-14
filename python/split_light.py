# Assign/create directories for each light chain and create fastas

import receptor_utils.simple_bio_seq as simple
import os


ldir = 'kos_et_al_2022_light_chain'
commands = []

# read v sequences, create directories if necessary, save seqs as fasta
for chain in ['IGK', 'IGL']:
    recs = simple.read_csv(os.path.join(ldir, f'{chain}V_Inferences.csv'))
    strain_map = {rec['Strain_ID']: rec['Strain_ID'].replace('/', '_').replace(' ', '') for rec in recs}

    for orig_strain, strain in strain_map.items():
        dir = strain
        if 'C57BL' in strain:
            dir = 'B6'
        if 'BALB' in strain:
            dir = 'BALBC'
        if os.path.isdir(dir):
            print(f'{strain} -> {dir}')
        else:
            os.mkdir(dir)
            print(f'{strain}: creating -> {dir}')

        if not os.path.isdir(os.path.join(dir, ldir)):
            os.mkdir(os.path.join(dir, ldir))

        if not os.path.isdir(os.path.join(dir, 'db')):
            os.mkdir(os.path.join(dir, 'db'))

        seqs = {}
        for rec in recs:
            if rec['Strain_ID'] == orig_strain:
                seqs[rec['Seq_Name']] = rec['Sequence'].upper()

        simple.write_fasta(seqs, os.path.join(dir, ldir, f'{chain}_V.fasta'))
        commands.append(f'cd {dir}')
        commands.append(f'rm db/mouse_{dir}_{chain}.csv')
        commands.append(f'python ../../iglabel/iglabel.py create db/mouse_{dir}_{chain}.csv')
        commands.append(f'python ../../iglabel/iglabel.py query db/mouse_{dir}_{chain}.csv {ldir}/{chain}_V.fasta results.csv actions.csv')
        commands.append(f'python ../../iglabel/iglabel.py add db/mouse_{dir}_{chain}.csv actions.csv "Kos et al., 2022"')
        commands.append(f'gap_inferred {ldir}/{chain}_V.fasta ../imgt/Mus_musculus_{chain}V_gapped.fasta {ldir}/{chain}_V_gapped.fasta')
        commands.append('cd ..')

commands = '\n'.join(commands)
with open('kos_light_chain_commands.bat', 'w') as f:
    f.write(commands)


