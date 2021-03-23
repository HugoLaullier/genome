from Bio import Entrez
from lxml import etree
import pandas as pd
import os
import shutil

# Entrez.email = 'adibou@example.com'  # Always tell NCBI who you are
# handle = Entrez.efetch(db='bioproject', id='PRJNA589917')
# record_xml = handle.read()
# print(record_xml)
# root = etree.fromstring(record_xml)
# tree = etree.ElementTree(root)
# print(tree)
# id_bioproj = tree.xpath('/RecordSet/DocumentSummary/Project/ProjectID/ArchiveID')[0].get('id')
# handle = Entrez.esearch(db='nuccore', term=id_bioproj+'[BioProject]')
# search_results = Entrez.read(handle)
# print(search_results)
# for result in search_results['IdList']:
#     entry = Entrez.efetch(db='nuccore', id=result, rettype='fasta')
#     this_seq_in_fasta = entry.read()
#     print(this_seq_in_fasta)

# exit()

def get_pb_access(organism, f_kingdom):
    first_row = True
    bp_access = 'none'
    for row in f_kingdom:
        if first_row:
            first_row=False
            continue
        parsed_row = row.split('\t')
        # print(parsed_row)
        if parsed_row[0].replace(' ','_').replace('/','_') == organism:
            bp_access = parsed_row[2]
    return bp_access

def write_cds(organism, kingdom, path, f_euk, f_vic, f_pla, f_pro):
    print(organism)
    bp_access = get_pb_access(organism, f_euk)
    print('F_EUK')
    if bp_access == 'none':
        bp_access = get_pb_access(organism, f_vic)
        print('F_VIC')
    if bp_access == 'none':
        bp_access = get_pb_access(organism, f_pla)
        print('F_PLA')
    if bp_access == 'none':
        bp_access = get_pb_access(organism, f_pro)
        print('F_PRO')

    if bp_access == 'none':
        print('ERROR')
        exit()
    
    print(bp_access)

    with open(path + '/' + organism + '_CDS.txt', 'w+') as cds:
        cds.write('Test')

if os.path.exists('../Results'):
    shutil.rmtree('../Results')

with open('../GENOME_REPORTS/overview.txt') as f, open('../GENOME_REPORTS/eukaryotes.txt') as euk, open('../GENOME_REPORTS/viruses.txt') as vir, open('../GENOME_REPORTS/plasmids.txt') as pla, open('../GENOME_REPORTS/prokaryotes.txt') as pro:
    first_row = True
    for row in f:
        if first_row:
            first_row=False
            continue
        parsed_row = row.split('\t')
        organism = parsed_row[0].replace(' ','_').replace('/','_')
        kingdom = parsed_row[1].replace(' ','_').replace('/','_')
        group = parsed_row[2].replace(' ','_').replace('/','_')
        subgroup = parsed_row[3].replace(' ','_').replace('/','_')
        path = '../Results/' + kingdom +'/' + group +'/' + subgroup +'/' + organism
        if not os.path.exists(path):
            os.makedirs(path)
            
        write_cds(organism, kingdom, path, euk, vir, pla, pro)
        with open(path + '/' + organism + '_intron.txt', 'w+') as intron:
            intron.write('Test')
        