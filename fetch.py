#!/usr/bin/env python3
from Bio import Entrez
from lxml import etree
import pandas as pd
import os
import shutil
from Bio import SeqIO
import pickle
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
import random
import string
import signal
from ftplib import FTP
import re
import time
import os.path
import datetime
from threading import Thread
import functools

save_pickle = False
DEBUG = False


def timeout(timeout):
    def deco(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            res = [Exception('function [%s] timeout [%s seconds] exceeded!' % (func.__name__, timeout))]
            def newFunc():
                try:
                    res[0] = func(*args, **kwargs)
                except Exception as e:
                    res[0] = e
            t = Thread(target=newFunc)
            t.daemon = True
            try:
                t.start()
                t.join(timeout)
            except Exception as je:
                print ('error starting thread')
                raise je
            ret = res[0]
            if isinstance(ret, BaseException):
                raise ret
            return ret
        return wrapper
    return deco

def reset_tree(root, progress = None, window = None):
    """
    reset the tree stored locally
    return organism_df
    """
    if (progress != None and window != None):
        progress['value'] = 0
        window.update_idletasks()
    # delete previous tree
    if os.path.exists('../Results'):
        shutil.rmtree('../Results')

    # update with new report from website
    t = os.path.getmtime("../GENOME_REPORTS/overview.txt")
    pred_time = datetime.datetime.fromtimestamp(t)
    now = datetime.datetime.now()

    if not(pred_time.month == now.month and pred_time.day == now.day and pred_time.year == now.year) :
        try: shutil.rmtree('../GENOME_REPORTS')
        except: pass
        os.mkdir("../GENOME_REPORTS")
        os.chdir('../GENOME_REPORTS')
        os.mkdir('IDS')

        with FTP('ftp.ncbi.nlm.nih.gov') as ftp:
            ftp.login()  # connecter au FTP

            ftp.cwd('genomes/GENOME_REPORTS')
            ftp.retrbinary('RETR overview.txt', open("overview.txt",'wb').write)

            ftp.cwd('IDS')  # changer de répertoire courant
            for filename in ftp.nlst():
                ftp.retrbinary('RETR '+ filename, open("IDS/" + filename, 'wb').write)


        # parse overview.txt
        organism_names = []
        organism_paths = []

        os.chdir('../script')
        with open('../GENOME_REPORTS/overview.txt') as f:
            first_row = True
            count_rows = 1
            for row in f:
                if DEBUG:
                    print(count_rows, " / 59674")

                if (progress != None and window != None and count_rows % 500 == 0):
                    progress['value'] = (count_rows/59674)*50
                    window.update_idletasks()
                count_rows += 1
                if first_row:
                    first_row=False
                    continue
                parsed_row = row.split('\t')

                try :
                    organism = parsed_row[0].replace(' ','_').replace('/','_')
                    kingdom = parsed_row[1].replace(' ','_').replace('/','_')
                    group = parsed_row[2].replace(' ','_').replace('/','_')
                    subgroup = parsed_row[3].replace(' ','_').replace('/','_')
                    path = '../Results/' + kingdom +'/' + group +'/' + subgroup +'/' + organism
                    organism_names.append(parsed_row[0])
                    organism_paths.append('../Results/' + kingdom +'/' + group +'/' + subgroup +'/')
                except IndexError : pass
        # parse ids files
        ids_files = os.listdir('../GENOME_REPORTS/IDS/')

        organism_names_ids = []
        organism_paths_ids = []
        organism_NC_ids = []
        i = 0
        for ids in ids_files:
            i += 1
            if (progress != None and window != None):
                    progress['value'] = 50 + (i/10)*50
                    window.update_idletasks()
            if DEBUG:
                print(str(i) + ' ' * (1 if i >= 10 else 2) + '/ ' + str(len(ids_files)) + ' : ' + ids)
            with open('../GENOME_REPORTS/IDS/' + ids) as f:
                for row in f:
                    parsed_row = row.replace('\n', '').split('\t')
                    if (parsed_row[1][0:2] != 'NC'):
                        continue
                    try:
                        index = organism_names.index(parsed_row[5])
                    except ValueError:
                        #if ids == "Archaea.ids" :
                        #    print(parsed_row)
                        for i, org in enumerate(organism_names) :
                            if org in parsed_row[5] :
                                index = i
                                break

                    try:
                        organism_NC_ids[organism_names_ids.index(organism_names[index])].append(parsed_row[1])
                    except ValueError:
                        organism_names_ids.append(organism_names[index])
                        organism_paths_ids.append(organism_paths[index])
                        organism_NC_ids.append([parsed_row[1]])
                        name = organism_names[index].replace(" ", "_")
                        name = name.replace("[", "_")
                        name = name.replace("]", "_")
                        name = name.replace(":", "_")
                        path = organism_paths[index] + name + "/"
                        if not os.path.exists(path):
                            os.makedirs(path)
        organism_df = pd.DataFrame({
                    "name":organism_names_ids,
                    "path":organism_paths_ids,
                    "NC":organism_NC_ids})
        # create pickle file saving the dataframe
        if not os.path.exists("../pickle"):
            os.makedirs("../pickle")
        with open("../pickle/organism_df", 'wb') as f:
            pickle.dump(organism_df, f)
    root.destroy()


def load_df_from_pickle():
    """
    load pickle dataframe and return it
    """
    try:
        with open("../pickle/organism_df", 'rb') as f:
            organism_df = pickle.load(f)
    except IOError:
        print("Pickle file not accessible")
        return reset_tree()
    if not os.path.exists('../Results'):
        for i in range(len(organism_df)):
            name = organism_df["name"][i].replace(" ", "_")
            name = name.replace("[", "_")
            name = name.replace("]", "_")
            name = name.replace(":", "_")
            path = organism_df["path"][i] + name + "/"
            if not os.path.exists(path):
                os.makedirs(path)
    return organism_df

def join(outfile, header_str, selected_region, feature_location, record_fasta, is_complement):
    write_buffer = header_str + feature_location
    if is_complement:
        feature_location = feature_location[16:-2]
    else:
        feature_location = feature_location[5:-1]
    
    x = feature_location.split(",")

    is_valid = True

    indexes = []
    for xi in x:
        xi = xi.split("..")
        try:
            indexes.append([int(xi[0]),int(xi[1])])
        except ValueError:
            is_valid = False

    # check si les index sont dans le bon ordre
    for i in range(len(indexes) - 1):
        if indexes[i] > indexes[i + 1]:
            return False
    indexes.sort(key=lambda r:r[0])


    if selected_region == "intron":
        fn = []
        for i in range(len(indexes) - 1):
            if(check_inf_sup(indexes[i][1],indexes[i+1][0]) == False):
                is_valid = False
            fn.append(FeatureLocation(indexes[i][1]-1, indexes[i+1][0]))
        if not is_valid:
            return
    else:
        fn = []
        for xi in indexes:
            if(check_inf_sup(xi[0],xi[1]) == False):
                is_valid = False
            else :
                fn.append(FeatureLocation(xi[0]-1, xi[1]))
        if not is_valid:
            return
    
    if len(fn) > 1:
        f = CompoundLocation(fn)
    else:
        f = SeqFeature(fn[0], type="domain")
    
    if selected_region != "intron":
        if is_complement:
            if DEBUG:
                print("COMPLEMENT JOIN")
                print(f.extract(record_fasta.seq).reverse_complement())
            write_buffer += "\n" +str(f.extract(record_fasta.seq).reverse_complement())
        else:
            if DEBUG:
                print("JOIN")
                print(f.extract(record_fasta.seq))
            write_buffer += "\n" +str(f.extract(record_fasta.seq))
    
    if len(fn) >= 1:
        for i in range(len(fn)):
            type_seq = ": exon "
            if selected_region == "intron":
                type_seq = ": intron "
            if is_complement:
                write_buffer += "\n" + header_str + feature_location + type_seq+ str(i+1)+ "\n" + str(SeqFeature(fn[i], type="domain").extract(record_fasta.seq).reverse_complement())
            else:
                write_buffer += "\n" + header_str + feature_location + type_seq+ str(i+1)+ "\n" + str(SeqFeature(fn[i], type="domain").extract(record_fasta.seq))
    outfile.write(write_buffer + "\n")
    return

def extract(outfile, header_str, selected_region, feature_location, record_fasta, is_complement):
    if selected_region == "intron":
        return
    write_buffer = header_str + feature_location + "\n"
    if is_complement:
        feature_location = feature_location[11:-1]
    x = feature_location.split("..")
    try:
        (int(x[0]),int(x[1]))
    except ValueError:
        return
    else:
        if(check_inf_sup(int(x[0]),int(x[1])) == False):
            return
        f = SeqFeature(FeatureLocation(int(x[0])-1, int(x[1])), type="domain")
        if is_complement:
            if DEBUG:
                print("COMPLEMENT")
                print(f.extract(record_fasta.seq).reverse_complement())
            write_buffer += str(f.extract(record_fasta.seq).reverse_complement())
        else:
            if DEBUG:
                print("EXTRACT")
                print(f.extract(record_fasta.seq))
            write_buffer += str(f.extract(record_fasta.seq))
    outfile.write(write_buffer + "\n")
    return

def handler():
    print("Timeout, NC skipped")
    raise Exception("end of time")

def load_data_from_NC(index, name, path, NC_list, selected_region):
    """
    download data of an organism from genbank using the API
    """
    Entrez.email = ''.join(random.choice(string.ascii_lowercase) for i in range(20))+'@'+''.join(random.choice(string.ascii_lowercase) for i in range(20))+ '.com'
    nb_region_found = 0
    print()
    print("downloading [" + name + "]")
    # remove duplicate NC
    NC_list = list(dict.fromkeys(NC_list))
    print(NC_list)
    for NC in NC_list:
        print("\t", NC)
        if DEBUG == True :
            print(str(NC) + " / " + str(len(NC_list)))
        name = name.replace(" ", "_")
        name = name.replace("[", "_")
        name = name.replace("]", "_")
        name = name.replace(":", "_")
        if DEBUG:
            print("NC id  =", NC)
            print("----------------------------")
        handle_fasta = Entrez.efetch(db="nucleotide", id=NC, rettype="fasta", retmode="text")
        '''
        signal.signal(signal.SIGALRM, handler)
        signal.alarm(10)
        try:
            record_fasta = SeqIO.read(handle_fasta, "fasta")
        except Exception:
            continue
        signal.alarm(0)
        '''
        record_fasta = timeout(timeout=10)(SeqIO.read)(handle_fasta, "fasta")

        if DEBUG:
            print(record_fasta)
            print("----------------------------")
            
        handle_fasta.close()
        handle_text = Entrez.efetch(db="nucleotide", id=NC, retmode="xml")
        record = Entrez.read(handle_text)
        handle_text.close()
        list_file = []
        for i in range(len(record[0]["GBSeq_feature-table"])):
            if DEBUG == True :
                print("\tfeature : " + str(i + 1) + " / " + str(len(record[0]["GBSeq_feature-table"])))
            feature_location = record[0]["GBSeq_feature-table"][i]["GBFeature_location"]
            feature_key = record[0]["GBSeq_feature-table"][i]["GBFeature_key"]
            if feature_key != selected_region and not (selected_region == "intron" and feature_key == "CDS"):
                continue
            nb_region_found += 1
            NC_filename = selected_region + "_" + str(name) + "_" + str(NC) + ".txt"
            if len(list_file) != 0 :
                if NC_filename not in list_file:
                    os.remove(path + name + "/" + NC_filename)
                    list_file.append(NC_filename)
            else :
                if os.path.isfile(path + name + "/" + NC_filename):
                    os.remove(path + name + "/" + NC_filename)
                list_file.append(NC_filename)
            with open(path + name + "/" + NC_filename, 'a+') as out:
                if DEBUG:
                    print(i+1, "/", len(record[0]["GBSeq_feature-table"]))
                    print(feature_location)
                
                header_str = selected_region + ' ' + name + ' ' + NC + ': '

                if feature_location.find("complement") != -1 and feature_location.find("join") != -1:
                    join(out, header_str, selected_region, feature_location, record_fasta, True)

                elif feature_location.find("complement") != -1:
                    extract(out, header_str, selected_region, feature_location, record_fasta, True)

                elif feature_location.find("join") != -1:
                    join(out, header_str, selected_region, feature_location, record_fasta, False)

                else:
                    extract(out, header_str, selected_region, feature_location, record_fasta, False)
            if os.stat(path + name + "/" + NC_filename).st_size == 0:
                        os.remove(path + name + "/" + NC_filename)
                        nb_region_found -= 1

    if nb_region_found == 0:
        print("Selected functional region not found for organism : [" + name + "]")
        return 0
    print(name + " downloaded")
    return nb_region_found

def check_inf_sup(inf,sup):
    if(inf <= sup):
        return True
    else:
        return False

# TODOLIST
# print chaque exons pour les joins ??? cf sujet
# update arbre plus rapidement
# pas de bleu (parce que!)
# Viruses missing files and directory ...