import part_a
import tools

def cross_reference_GB_UniProt(UniProt_file_path, GB_file_path):
    '''
    Making Cross-reference between UniProt and GB
    the Cross-reference is based on 'GOA' id
    print a report of the Cross-reference
    :param UniProt_file_path: the path to the UniProt file
    :param GB_file_path: the path to the GB file
    :return: nothing
    '''
    uniprot_protein_dictionary = tools.get_dictionary_of_protain_UniProt_file(UniProt_file_path)
    gb_protein_dictionary = tools.get_dictionary_of_protain_from_GB_file(GB_file_path)

    count_uniprot_in_gb = 0
    for key, val in uniprot_protein_dictionary.items():
        if key in gb_protein_dictionary:
            count_uniprot_in_gb += 1

    count_gb_in_uniprot = 0
    for key, val in gb_protein_dictionary.items():
        if key in uniprot_protein_dictionary:
            count_gb_in_uniprot += 1

    print("---------count uniprot in gb--------")
    print("count_uniprot_in_gb = " + str(count_uniprot_in_gb))
    print("uniprot len = " + str(len(uniprot_protein_dictionary)))

    print("---------count gb in uniprot--------")
    print("count_gb_in_uniprot = " + str(count_gb_in_uniprot))
    print("gb len = " + str(len(gb_protein_dictionary)))


def main_b():
    print("in b")
    # ------ b.a ------
    cross_reference_GB_UniProt("uniprot-organism Bacillus+subtilis+(strain+168)+[224308] .fasta", "BS168.gb")
    # ------ b.b ------
    print(str(tools.get_transmembrane_from_UniProt_file("uniprot-organism Bacillus+subtilis+(strain+168)+[224308] .fasta")))




