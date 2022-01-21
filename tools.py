
from Bio import GenBank ,SeqIO



def import_data_from_DG_file( file_path):
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