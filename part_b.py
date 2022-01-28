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
    GOA_uniprot_protein_dictionary, Name_uniprot_protein_dictionary = tools.get_dictionary_of_protain_UniProt_file(UniProt_file_path)
    GOA_gb_protein_dictionary, Name_gb_protein_dictionary= tools.get_dictionary_of_protain_from_GB_file(GB_file_path)

    count_uniprot_in_gb_by_GOA = 0
    for key, val in GOA_uniprot_protein_dictionary.items():
        if key in GOA_gb_protein_dictionary:
            count_uniprot_in_gb_by_GOA += 1

    count_uniprot_in_uniprot_by_Name = 0
    for key, val in Name_uniprot_protein_dictionary.items():
        if key in Name_gb_protein_dictionary:
            count_uniprot_in_uniprot_by_Name += 1

    count_gb_in_uniprot_by_GOA = 0
    for key, val in GOA_gb_protein_dictionary.items():
        if key in GOA_uniprot_protein_dictionary:
            count_gb_in_uniprot_by_GOA += 1

    count_gb_in_uniprot_by_Name = 0
    for key, val in Name_gb_protein_dictionary.items():
        if key in Name_uniprot_protein_dictionary:
            count_gb_in_uniprot_by_Name += 1

    print("---------count uniprot in gb--------")
    print("count_uniprot_in_gb by GOA = " + str(count_uniprot_in_gb_by_GOA))
    print("count_uniprot_in_gb by Name = " + str(count_uniprot_in_uniprot_by_Name))
    print("uniprot GOA len = " + str(len(GOA_uniprot_protein_dictionary)))
    print("uniprot Name len = " + str(len(Name_uniprot_protein_dictionary)))

    print("---------count gb in uniprot--------")
    print("count_gb_in_uniprot by GOA = " + str(count_gb_in_uniprot_by_GOA))
    print("count_gb_in_uniprot by Name = " + str(count_gb_in_uniprot_by_Name))
    print("gb GOA len = " + str(len(GOA_gb_protein_dictionary)))
    print("gb Name len = " + str(len(Name_gb_protein_dictionary)))






def research_transmembrane(UniProt_file_path):
    _,name_dictionary_to_return = tools.get_transmembrane_from_UniProt_file(UniProt_file_path)
    list_len = []
    for key,val in name_dictionary_to_return.items():
        list_len.append(len(val))
    all_values = name_dictionary_to_return.values()
    print("---------------------------")
    print("max length of protein containing transmembrane part = " + str(max(list_len)))
    print("min length of protein containing transmembrane part = " + str(min(list_len)))
    print("AVG length of protein containing transmembrane part = " + str(sum(list_len) / len(list_len)))

    tools.create_histogram(list_len, "Histogram - length of protein containing Transmembrane part")

    map_distribution = {}
    for key,val in name_dictionary_to_return.items():
        count_hydrophobic = 0
        for char in val:
            if char in tools.dictionary_Hydrophobic_amino_acids:
                count_hydrophobic += 1
        map_distribution[key] = count_hydrophobic/len(val)
    print("Distribution of hydrophobic amino acids for each one of the transmembrane protein : (GOA,Distribution)")
    print(str(map_distribution))
    print("total avg of hydrophobic amino acids in those transmembrane protein = " + str(sum(map_distribution.values())/len(name_dictionary_to_return)))



def research_GC(UniProt_file_path,gb_file_path):
    _, dictionary_of_transmembrane = tools.get_transmembrane_from_UniProt_file(UniProt_file_path)
    _, A_Name_gb_protein_dictionary = tools.get_dictionary_of_protain_from_GB_file(gb_file_path)

    B = {}
    global_len = 0
    GC_count = 0
    for key, val in dictionary_of_transmembrane.items():
        if key in A_Name_gb_protein_dictionary:
            B[key] = val
            global_len += len(val)
            GC_count += val.count('G')
            GC_count += val.count('C')

    print("GC Distribution in B group is  = " + str((GC_count/global_len)*100) +"%")










def main_b():
    print("in b")

    uniprot_file_path = "uniprot-organism Bacillus+subtilis+(strain+168)+[224308] .fasta"
    gb_file_path = "BS168.gb"

    # ------ b.a ------
    cross_reference_GB_UniProt(uniprot_file_path, gb_file_path )
    # ------ b.b ------
    research_transmembrane(uniprot_file_path)
    # ------ b.b ------
    research_GC(uniprot_file_path,gb_file_path)





