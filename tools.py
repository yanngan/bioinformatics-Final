
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

def get_dictionary_of_protain_from_GB_file(file_path):
    '''
    :param: file_path: the path to the DB file
    :return:  dictionary: key = GOA identifier, value = protain seq
    '''
    count_entry = 0
    input_handle = open(file_path, "r")
    GOA_dictionary_to_return = {}
    Name_dictionary_to_return = {}
    for seq_record in SeqIO.parse(input_handle, "genbank"):
        print("Dealing with GenBank record %s" % seq_record.id)
        for seq_feature in seq_record.features :
            if seq_feature.type == "CDS":
                assert len(seq_feature.qualifiers['translation']) == 1
                count_entry += 1
                if any("GOA" in s for s in seq_feature.qualifiers['db_xref']):
                    GOA_dictionary_to_return[seq_feature.qualifiers['db_xref'][2].split(":")[1]] = seq_feature.qualifiers['translation'][0]
                if seq_feature.qualifiers.get('gene') != None:
                    Name_dictionary_to_return[seq_feature.qualifiers.get('gene')[0]] = seq_feature.qualifiers['translation'][0]
    print("len of full CDS GB records = " + str(count_entry))
    return GOA_dictionary_to_return,Name_dictionary_to_return

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


def create_histogram(val_list, title="histogram"):
    n, bins, patches = plt.hist(x=val_list, bins='auto', color='#0504aa')
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title(title)
    plt.show()


