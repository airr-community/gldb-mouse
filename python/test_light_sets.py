# download sets from OGRDB and test against Junstin's file

import os
import requests
import receptor_utils.simple_bio_seq as simple

ldir = 'kos_et_al_2022_light_chain'
root_fetch = 'https://ogrdb.airr-community.org/download_germline_set/26/ungapped'

with open(os.path.join(ldir, 'test', 'germline_sets.html'), 'r') as fi:
    markup = fi.read()

addrs = markup.split('<a href="')
addrs = [x.split('</a>')[0] for x in addrs]

germline_sets = {}

master = {}
master['IGKV'] = simple.read_csv(os.path.join(ldir, 'IGKV_inferences.csv'))
master['IGLV'] = simple.read_csv(os.path.join(ldir, 'IGLV_inferences.csv'))

ogrdb = {}
ogrdb['IGKV'] = {}
ogrdb['IGLV'] = {}

for addr in addrs:
    if 'germline_set' in addr and 'span' not in addr:
        addr = addr.split('">')
        germline_sets[addr[1]] = addr[0].replace('/germline_set/', '')

for set_name, set_id in germline_sets.items():
    file_name = set_name.replace(' ', '_').replace('/', '_')
    set_name, chain = set_name.split(' ')

    if chain not in ['IGKV', 'IGLV']:
        continue

    r = requests.get(f'https://ogrdb.airr-community.org/download_germline_set/{set_id}/ungapped')
    open(os.path.join(ldir, 'test', f'{file_name}.fasta'), 'wb').write(r.content)
    ogrdb[chain][set_name] = simple.read_fasta(os.path.join(ldir, 'test', f'{file_name}.fasta'))
    print(set_name)

for chain in ['IGKV', 'IGLV']:
    for set_name, set_data in ogrdb[chain].items():
        for s_name, s_seq in set_data.items():
            found = False
            for rec in master[chain]:
                if rec['Strain_ID'].replace(' ', '') == set_name and rec['Sequence'].upper() == s_seq:
                    found = True
            if not found:
                print(f'{chain} {set_name} {s_name}: sequence found in OGRDB but not in Justin\'s file')

    for rec in master[chain]:
        if rec['Strain_ID'] not in ogrdb[chain]:
            # print(f'{chain} {rec["Strain_ID"]}: set found in Justin\'s file but not in OGRDB')
            pass
        else:
            for s_name, s_seq in ogrdb[chain][rec['Strain_ID']].items():
                found = False
                if rec['Sequence'].upper() == s_seq:
                    found = True
                    break
            if not found:
                print(f'{chain} {rec["Strain_ID"]} {rec["Seq_Name"]}: sequence found in Justin\'s file but not in OGRDB')

