
from Bio import GenBank ,SeqIO



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
    dictionary_to_return = {}
    for protain in uniprot_protein_list:
        dictionary_to_return[protain.name.split("|")[1]] = protain.seq
        if protain.name.split("|")[1] == "O32050":
            print("hello")
    return dictionary_to_return

def get_dictionary_of_protain_from_GB_file(file_path):
    '''
    :param: file_path: the path to the DB file
    :return:  dictionary: key = GOA identifier, value = protain seq
    '''
    debug = 0
    input_handle = open(file_path, "r")
    dictionary_to_return = {}
    for seq_record in SeqIO.parse(input_handle, "genbank"):
        print("Dealing with GenBank record %s" % seq_record.id)
        for seq_feature in seq_record.features :
            if seq_feature.type == "CDS":
                assert len(seq_feature.qualifiers['translation']) == 1
                debug += 1
                # output_handle = open("faa_filename.txt", "a")
                # output_handle.write(">%s from %s\n%s\n" % (
                #     str(seq_feature),#.qualifiers['locus_tag'][0],
                #     seq_record.name,
                #     seq_feature.qualifiers['translation'][0]))
                if any("GOA" in s for s in seq_feature.qualifiers['db_xref']):
                    dictionary_to_return[seq_feature.qualifiers['db_xref'][2].split(":")[1]] = seq_feature.qualifiers['translation'][0]
    print("len of full CDS DB records = " + str(debug))
    print("len of CDS DB records having in 'db_xref' the identifier 'GOA' = " + str(len(dictionary_to_return)))
    return dictionary_to_return

def get_transmembrane_from_UniProt_file(fasta_path):
    uniprot_protein_list = import_data_from_UniProt_server(fasta_path)
    dictionary_to_return = {}
    for protain in uniprot_protein_list:
        if protain.description.find("transmembrane") != -1:
            dictionary_to_return[protain.name.split("|")[1]] = protain.seq
    return dictionary_to_return
