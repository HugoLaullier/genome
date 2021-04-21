from Bio import Entrez
from lxml import etree
import pandas as pd
import os
import shutil
from Bio import SeqIO

# GREP overview.txt

Entrez.email = 'adibou@example.com'  # Always tell NCBI who you are
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
    #print(organism)
    bp_access = get_pb_access(organism, f_euk)
    #print('F_EUK')
    if bp_access == 'none':
        bp_access = get_pb_access(organism, f_vic)
    #    print('F_VIC')
    if bp_access == 'none':
        bp_access = get_pb_access(organism, f_pla)
    #    print('F_PLA')
    if bp_access == 'none':
        bp_access = get_pb_access(organism, f_pro)
    #    print('F_PRO')

    if bp_access == 'none':
        pass
    #    print('ERROR')
        #exit()
    
    #print(bp_access)

    with open(path + '/' + organism + '_CDS.txt', 'w+') as cds:
        cds.write('Test')

if os.path.exists('../Results'):
    shutil.rmtree('../Results')


organisms_overview = []
with open('../GENOME_REPORTS/overview.txt') as f, open('../GENOME_REPORTS/eukaryotes.txt') as euk, open('../GENOME_REPORTS/viruses.txt') as vir, open('../GENOME_REPORTS/plasmids.txt') as pla, open('../GENOME_REPORTS/prokaryotes.txt') as pro:
    first_row = True
    count_rows = 1
    for row in f:
        print(count_rows, " / 59674")
        count_rows += 1
        if first_row:
            first_row=False
            continue
        parsed_row = row.split('\t')
        organism = parsed_row[0].replace(' ','_').replace('/','_')
        kingdom = parsed_row[1].replace(' ','_').replace('/','_')
        group = parsed_row[2].replace(' ','_').replace('/','_')
        subgroup = parsed_row[3].replace(' ','_').replace('/','_')
        path = '../Results/' + kingdom +'/' + group +'/' + subgroup +'/' + organism
        organisms_overview.append(parsed_row[0])
        if not os.path.exists(path):
            os.makedirs(path)
        
        write_cds(organism, kingdom, path, euk, vir, pla, pro)
        with open(path + '/' + organism + '_intron.txt', 'w+') as intron:
            intron.write('Test')

# GREP IDS
ids_files = os.listdir('../GENOME_REPORTS/IDS/')
organisms_ids = []
individual_dict = {}
i = 0
for ids in ids_files:
    i += 1
    print(str(i) + ' ' * (1 if i >= 10 else 2) + '/ ' + str(len(ids_files)) + ' : ' + ids)
    with open('../GENOME_REPORTS/IDS/' + ids) as f:
        for row in f:
            parsed_row = row.replace('\n', '').split('\t')
            if (parsed_row[1][0:2] != 'NC'):
                continue
            if (not (parsed_row[5] in organisms_overview)):
                continue
            #print(parsed_row)
            organisms_ids.append(parsed_row[5])
            if parsed_row[5] in individual_dict.keys():
                individual_dict[parsed_row[5]].append(parsed_row[1])
            else :
                individual_dict[parsed_row[5]] = [parsed_row[1]]
print()
print(len(individual_dict))
# print(ids_files)
# print(individual_dict)

print("----------------------------")

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

print(len(organisms_overview))
print(len(individual_dict))
print("Abelson murine leukemia virus" in organisms_ids)
print("Abelson murine leukemia virus" in organisms_overview)
# print(intersection(organisms_overview, organisms_ids))
print(len(intersection(organisms_overview, organisms_ids)))
# path = '../Results/'
# Entrez.email = 'goret@example.com'
# i = 0
# for ind in individual_dict.keys():
#     i += 1
#     print(str(i) + ' / ' + str(len(individual_dict)), str(i/len(individual_dict)*100)+'%')
#     handle = Entrez.efetch(db="nucleotide", id=individual_dict[ind], retmode="xml")
#     print(individual_dict[ind])
#     records = Entrez.read(handle)
#     handle.close()
#     final_path = path + records[0]["GBSeq_taxonomy"].replace('; ','/')+'/'+records[0]["GBSeq_organism"]
#     print(records[0].keys())
#     print('NOTE : ' + str(records[0]["GBSeq_project"]))
#     # if not os.path.exists(final_path):
#     #     os.makedirs(final_path)
    


#print(handle.read())
# print(handle.read())
# print()
# print(record.id)
# print()
# print(record.name)
# print()
# print(record.description)
# print()
# print(len(record.features))
# print()
# print(repr(record.seq))
# print()
# print(record)