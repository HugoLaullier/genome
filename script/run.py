from Bio import Entrez
from lxml import etree
import pandas as pd
import os
import shutil
from Bio import SeqIO
import pickle
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

save_pickle = False
debug = True
Entrez.email = 'porcinet@example.com'

def reset_tree():
    """
    reset the tree stored locally
    return organism_df
    """
    # delete previous tree
    if os.path.exists('../Results'):
        shutil.rmtree('../Results')
    
    # parse overview.txt
    organism_names = []
    organism_paths = []
    with open('../GENOME_REPORTS/overview.txt') as f:
        first_row = True
        count_rows = 1
        for row in f:
            if debug:
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
            organism_names.append(parsed_row[0])
            organism_paths.append('../Results/' + kingdom +'/' + group +'/' + subgroup +'/')
            if not os.path.exists(path):
                os.makedirs(path)

    # parse ids files
    ids_files = os.listdir('../GENOME_REPORTS/IDS/')

    organism_names_ids = []
    organism_paths_ids = []
    organism_NC_ids = []

    i = 0
    for ids in ids_files:
        i += 1
        if debug:
            print(str(i) + ' ' * (1 if i >= 10 else 2) + '/ ' + str(len(ids_files)) + ' : ' + ids)
        with open('../GENOME_REPORTS/IDS/' + ids) as f:
            for row in f:
                parsed_row = row.replace('\n', '').split('\t')
                if (parsed_row[1][0:2] != 'NC'):
                    continue
                try:
                    index = organism_names.index(parsed_row[5])
                except ValueError:
                    continue
                try:
                    organism_NC_ids[organism_names_ids.index(organism_names[index])].append(parsed_row[1])
                except ValueError:
                    organism_names_ids.append(organism_names[index])
                    organism_paths_ids.append(organism_paths[index])
                    organism_NC_ids.append([parsed_row[1]])
    
    organism_df = pd.DataFrame({
                "name":organism_names_ids,
                "path":organism_paths_ids,
                "NC":organism_NC_ids})
    # create pickle file saving the dataframe
    if not os.path.exists("../pickle"):
        os.makedirs("../pickle")
    with open("../pickle/organism_df", 'wb') as f:
        pickle.dump(organism_df, f)
    return organism_df

def load_df_from_pickle():
    """
    load pickle dataframe and return it
    """
    with open("../pickle/organism_df", 'rb') as f:
        organism_df = pickle.load(f)
    return organism_df

def load_data_from_NC(index, name, path, NC_list):
    """
    download data of an organism from genbank using the API
    """
    NC_i = 1
    for NC in NC_list:
        name = name.replace(" ", "_")
        NC_filename = str(name) + "_CDS_NC_" + str(NC_i) + ".txt"
        NC_i += 1
        if debug:
            print("NC id  =", NC)
            print("----------------------------")
        handle_fasta = Entrez.efetch(db="nucleotide", id=NC, rettype="fasta", retmode="text")
        record_fasta = SeqIO.read(handle_fasta, "fasta")
        if debug:
            print(record_fasta)
            print("----------------------------")
        handle_fasta.close()
        handle_text = Entrez.efetch(db="nucleotide", id=NC, retmode="xml")
        record = Entrez.read(handle_text)
        handle_text.close()
        with open(path + name + "/" + NC_filename, 'w+') as out:
            for i in range(len(record[0]["GBSeq_feature-table"])):
                feature_location = record[0]["GBSeq_feature-table"][i]["GBFeature_location"]
                if debug:
                    print(i+1, "/", len(record[0]["GBSeq_feature-table"]))
                    print(feature_location)
                # TODO Tests sur les regions (partie 2.3)
                out.write("CDS " + feature_location + "\n")
                if feature_location.find("complement")!= -1 and feature_location.find("join") != -1:
                    feature_location = feature_location[16:-1]
                    x = feature_location.split(",")
                    fn = []
                    for xi in x:
                        xi = xi.split("..")
                        fn.append(FeatureLocation(int(xi[0]), int(xi[1])))
                    f = CompoundLocation(fn)
                    if debug:
                        print("COMPLEMENT JOIN")
                        print(f.extract(record_fasta.seq).complement())
                    out.write(str(f.extract(record_fasta.seq).complement()))
                
                elif feature_location.find("complement")!= -1:
                    feature_location = feature_location[11:-1]
                    x = feature_location.split("..")
                    f = SeqFeature(FeatureLocation(int(x[0]), int(x[1])), type="domain")
                    if debug:
                        print("COMPLEMENT")
                        print(f.extract(record_fasta.seq).complement())
                    out.write(str(f.extract(record_fasta.seq).complement()))
                
                elif feature_location.find("join") != -1:
                    feature_location = feature_location[5:-1]
                    x = feature_location.split(",")
                    fn = []
                    for xi in x:
                        xi = xi.split("..")
                        fn.append(FeatureLocation(int(xi[0]), int(xi[1])))
                    f = CompoundLocation(fn)
                    if debug:
                        print("JOIN")
                        print(f.extract(record_fasta.seq))
                    out.write(str(f.extract(record_fasta.seq)))
                
                else:
                    x = feature_location.split("..")
                    f = SeqFeature(FeatureLocation(int(x[0]), int(x[1])), type="domain")
                    if debug:
                        print("EXTRACT")
                        print(f.extract(record_fasta.seq))
                    out.write(str(f.extract(record_fasta.seq)))
            out.write("\n")

if __name__ == "__main__":
    if (save_pickle): # reset local tree and pickle
        organism_df = reset_tree()
    else: # load data from pickle
        organism_df = load_df_from_pickle()
    
    if debug:
        print("----------------------------")
        print(organism_df.head(5))
        print(organism_df.tail(5))
        print("----------------------------")
    
    for (index, name, path, NC_list) in organism_df.itertuples():
        load_data_from_NC(index, name, path, NC_list)
        if debug: # only load the first organism
            exit()
