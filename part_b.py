import part_a
import tools
import numpy as np
import matplotlib.pyplot as plt

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
    GOA_gb_protein_dictionary, Name_gb_protein_dictionary,_= tools.get_dictionary_of_protain_from_GB_file(GB_file_path)

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
    print("Distribution of hydrophobic amino acids for each one of the transmembrane protein : (Name,Distribution)")
    print(str(map_distribution))
    print("total avg of hydrophobic amino acids in those transmembrane protein = " + str(sum(map_distribution.values())/len(name_dictionary_to_return)))



def research_GC(UniProt_file_path,gb_file_path):
    _, dictionary_of_transmembrane = tools.get_transmembrane_from_UniProt_file(UniProt_file_path)
    _, A_Name_gb_protein_dictionary,_ = tools.get_dictionary_of_protain_from_GB_file(gb_file_path)

    B = {}
    global_len_B = 0
    GC_count_B = 0
    list_len_pro_in_B = []
    for key, val in dictionary_of_transmembrane.items():
        if key in A_Name_gb_protein_dictionary:
            B[key] = val
            global_len_B += len(val)
            GC_count_B += val.count('G')
            GC_count_B += val.count('C')
    print("GC Distribution in B group is  = " + str((GC_count_B / global_len_B) * 100) + "%")

    A_not_in_B = {}
    global_len_A_not_in_B = 0
    GC_count_A = 0
    GC_count_A_not_in_B = 0
    for key, val in A_Name_gb_protein_dictionary.items() :
        if key not in dictionary_of_transmembrane:
            A_not_in_B[key] = val
            global_len_A_not_in_B += len(val)
            GC_count_A_not_in_B += val.count('G')
            GC_count_A_not_in_B += val.count('C')
        GC_count_A += val.count('G')
        GC_count_A += val.count('C')




    res_B = tools.get_GC_and_len_statistic(B,GC_count_B,global_len_B)
    res_A = tools.get_GC_and_len_statistic(A_Name_gb_protein_dictionary,GC_count_A,global_len_B)
    res_A_not_in_B = tools.get_GC_and_len_statistic(A_not_in_B,GC_count_A_not_in_B,global_len_B)
    print("---len---")
    print("len min = "+str(res_B['len']["min"]))
    print("len max = "+str(res_B['len']["max"]))
    print("len avg = "+str(res_B['len']["avg"]))
    print("len median = "+str(res_B['len']["median"]))

    print("---gc---")
    print("gc min = "+str(res_B['gc']["min"]))
    print("gc max = "+str(res_B['gc']["max"]))
    print("gc avg = "+str(res_B['gc']["avg"]))
    print("gc median = "+str(res_B['gc']["median"]))

    # tools.create_histogram(res_B['lens'], 'lens')
    # tools.create_histogram(res_B['list_GC_count'], 'lens')


    max_len = max(res_B['len']["max"],res_A['len']["max"],res_A_not_in_B['len']["max"])
    max_len_GC = max(res_B['gc']["max"],res_A['gc']["max"],res_A_not_in_B['gc']["max"])
    to_histograms_len = {
        "lens A": res_A['lens'],
        "lens B": res_B['lens'],
        "lens A not in B": res_A_not_in_B['lens']
    }

    to_histograms_count_GC = {
        "count BC A": res_A['list_GC_count'],
        "count BC B": res_B['list_GC_count'],
        "count BC A not in B": res_A_not_in_B['list_GC_count']
     }

    bins_max_len = np.linspace(0, max_len,50)
    bins_max_len_GC = np.linspace(0, max_len_GC,50)

    index = 1
    for key, value in to_histograms_len.items():
        plt.subplot(2, 2, index)
        tools.create_histogram(value, key, False,bins_max_len)
        index += 1

    plt.subplot(2, 2, index)
    plt.hist([to_histograms_len['lens B'],to_histograms_len['lens A not in B']], bins_max_len, alpha=0.5, label=['lens B','lens A not in B'])
    plt.title("lens B & lens A not in B")
    plt.legend(loc='upper right')
    plt.show()


    index = 1
    for key, value in to_histograms_count_GC.items():
        plt.subplot(2, 2, index)
        tools.create_histogram(value, key, False,bins_max_len_GC)
        index += 1

    plt.subplot(2, 2, index)
    plt.hist([to_histograms_count_GC['count BC B'],to_histograms_count_GC['count BC A not in B']], bins_max_len_GC, alpha=0.5, label=['count BC B','count BC A not in B'])
    plt.title("count GC B & count A not in B")
    plt.legend(loc='upper right')
    plt.show()











def main_b():
    print("in b")

    uniprot_file_path = "uniprot-organism Bacillus+subtilis+(strain+168)+[224308] .fasta"
    gb_file_path = "BS168.gb"

    # ------ b.a ------
    tools.printSection("B","1")
    cross_reference_GB_UniProt(uniprot_file_path, gb_file_path )
    # ------ b.b ------
    tools.printSection("B", "2")
    research_transmembrane(uniprot_file_path)
    # ------ b.b ------
    tools.printSection("B", "3")
    research_GC(uniprot_file_path,gb_file_path)





