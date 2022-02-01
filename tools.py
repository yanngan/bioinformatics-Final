
from Bio import GenBank ,SeqIO
import matplotlib.pyplot as plt
import numpy as np


# ---------- const ------------------
dictionary_Hydrophobic_amino_acids = {"A", "P", "L", "I", "V", "M", "P", "T"}




def import_data_from_DB_file( file_path):
    '''
    import data from local DG file
    :param file_path: the path leading to the file
    :return: a list of GenBank records
    '''
    list_records = []
    with open(file_path) as handle:
        for record in GenBank.parse(handle):
            list_records.append(record)
    return list_records

def import_data_from_UniProt_server(fasta_path):
    '''
    import data from UniProt local file (FASTA format)
    :param fasta_path: the path leading to the file
    :return: a list of GenBank records
    '''
    list_records = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        list_records.append(record)
    return list_records

def get_dictionary_of_protain_UniProt_file(fasta_path):
    '''
    :param fasta_path: the path to the fasta file
    :return: dictionary: key = GOA identifier, value = protain seq
    '''
    uniprot_protein_list = import_data_from_UniProt_server(fasta_path)
    GOA_dictionary_to_return = {}
    Name_dictionary_to_return = {}
    for protain in uniprot_protein_list:
        GOA_dictionary_to_return[protain.name.split("|")[1]] = protain.seq
        for val in protain.description.split():
            if "GN=" in val:
                Name_dictionary_to_return[val[3:]] = protain.seq
    return GOA_dictionary_to_return,Name_dictionary_to_return

def get_dictionary_of_protain_from_GB_file(file_path, type_to_return = "amino_acid"):
    '''
    :param: file_path: the path to the DB file
    :param: type_to_return: the value to put in the dictionary value
    :return:  dictionary: key = GOA identifier, value = protain seq
    '''
    count_entry = 0
    input_handle = open(file_path, "r")
    GOA_dictionary_to_return = {}
    Name_dictionary_to_return = {}
    all_protein = {}
    for seq_record in SeqIO.parse(input_handle, "genbank"):
        #print("Dealing with GenBank record %s" % seq_record.id)
        for seq_feature in seq_record.features :
            if seq_feature.type == "CDS":
                assert len(seq_feature.qualifiers['translation']) == 1
                count_entry += 1
                if type_to_return == "amino_acid":
                    if 'db_xref' in seq_feature.qualifiers and any("GOA" in s for s in seq_feature.qualifiers['db_xref']):
                        GOA_dictionary_to_return[seq_feature.qualifiers['db_xref'][2].split(":")[1]] = seq_feature.qualifiers['translation'][0]
                    if seq_feature.qualifiers.get('gene') != None:
                        Name_dictionary_to_return[seq_feature.qualifiers.get('gene')[0]] = seq_feature.qualifiers['translation'][0]
                    all_protein[len(all_protein)]=[seq_feature.qualifiers['translation'][0],seq_feature.location,seq_feature]
                else:
                    temp_seq = seq_record.seq[seq_feature.location.start:seq_feature.location.end]
                    dna = temp_seq if seq_feature.strand == 1 else temp_seq.reverse_complement()
                    if 'db_xref' in seq_feature.qualifiers and any("GOA" in s for s in seq_feature.qualifiers['db_xref']):
                        GOA_dictionary_to_return[seq_feature.qualifiers['db_xref'][2].split(":")[1]] = dna
                    if seq_feature.qualifiers.get('gene') != None:
                        Name_dictionary_to_return[seq_feature.qualifiers.get('gene')[0]] = dna
                    all_protein[len(all_protein)]=[seq_feature.qualifiers['translation'][0],seq_feature.location,seq_feature]



    #print("len of full CDS GB records = " + str(count_entry))
    return GOA_dictionary_to_return,Name_dictionary_to_return,all_protein

def get_transmembrane_from_UniProt_file(fasta_path):
    uniprot_protein_list = import_data_from_UniProt_server(fasta_path)
    GOA_dictionary_to_return = {}
    Name_dictionary_to_return = {}
    for protain in uniprot_protein_list:
        if protain.description.find("transmembrane") != -1:
            GOA_dictionary_to_return[protain.name.split("|")[1]] = protain.seq
            for val in protain.description.split():
                if "GN=" in val:
                    Name_dictionary_to_return[val[3:]] = protain.seq
    return GOA_dictionary_to_return,Name_dictionary_to_return


def get_GC_and_len_statistic(dic,GC_count,global_len):
    #   avg
    dic_len_avg = 0
    dic_count_GC_avg = 0
    #   median
    dic_len_median = 0
    dic_count_GC_median = 0
    #   min
    dic_len_min = len(list(dic.values())[0])
    dic_count_GC_min = GC_count
    #   max
    dic_len_max = len(list(dic.values())[0])
    dic_count_GC_max = 0

    sum_for_len_avg = 0
    lens = []
    list_GC_count = []
    for protein in dic.values():
        temp_GC_count = protein.count('G')
        temp_GC_count += protein.count('C')
        if len(protein) > dic_len_max:
            dic_len_max = len(protein)
        if len(protein) < dic_len_min:
            dic_len_min = len(protein)
        if temp_GC_count > dic_count_GC_max:
            dic_count_GC_max = temp_GC_count
        if temp_GC_count < dic_count_GC_min:
            dic_count_GC_min = temp_GC_count
        sum_for_len_avg += len(protein)
        lens.append(len(protein))
        list_GC_count.append(temp_GC_count)
    dic_len_avg = sum_for_len_avg / len(dic)
    lens.sort()
    dic_len_median = lens[int(len(dic) / 2)]

    dic_count_GC_avg = (GC_count / global_len) * 100
    list_GC_count.sort()
    dic_count_GC_median = lens[int(len(list_GC_count) / 2)]

    return {
        "len": {
            "avg": dic_len_avg,
            "median": dic_len_median,
            "min": dic_len_min,
            "max": dic_len_max
        },
        "gc": {
            "avg": dic_count_GC_avg,
            "median": dic_count_GC_median,
            "min": dic_count_GC_min,
            "max": dic_count_GC_max
        },
        "lens": lens,
        "list_GC_count": list_GC_count,
    }


def printSection(part,sub_par):
    print("--------------------------------- part "+part+"---------------------------------")
    print("--------------------------------sub-par "+sub_par+"---------------------------------")

def create_histogram(val_list, title="histogram",show = True,bins ='auto'):
    n, bins, patches = plt.hist(x=val_list, bins=bins, color='#0504aa')
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title(title)
    if show:
        plt.show()
