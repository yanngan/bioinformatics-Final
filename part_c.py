from Bio import Entrez, SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
genetic_code=""

def NewCode(currentCode, changeAT):
    nextlatter = nextLatter[currentCode[changeAT]]
    currentCode = currentCode[0:changeAT] + nextlatter + currentCode[changeAT + 1:3]
    return currentCode
    pass


nextLatter={
    'A':'C', 'C':'G', 'G':'T', 'T':'A'
}

gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}



import tools
def open_GB( GB_file_path):

    gb_protein_dictionary = tools.import_data_from_DB_file(GB_file_path)
    return  gb_protein_dictionary[0].sequence

# def main_c():
#     print("in c")
#     #a    ----->
#     open_GB("corona_2020.gb")
#     print(count_mutation_by_type(1,"synonymous"))
def countSismogramAt(kodon):
    count=0
    #first latter
    if gencode[kodon]==gencode[kodon[0:2]+nextLatter[kodon[2]]]:
        count += 1
    if gencode[kodon]==gencode[kodon[0:2]+nextLatter[nextLatter[kodon[2]]]]:
        count += 1
    if gencode[kodon]==gencode[kodon[0:2]+nextLatter[nextLatter[nextLatter[kodon[2]]]]]:
        count+=1

    #sec latter
    if gencode[kodon]==gencode[kodon[0:1]+nextLatter[kodon[1]]+kodon[2]]:
        count += 1
    if gencode[kodon]==gencode[kodon[0:1]+nextLatter[nextLatter[kodon[1]]]+kodon[2]]:
        count += 1
    if gencode[kodon]==gencode[kodon[0:1]+nextLatter[nextLatter[nextLatter[kodon[1]]]]+kodon[2]]:
        count+=1

    #third latter
    if gencode[kodon]==gencode[nextLatter[kodon[0]]+kodon[1:3]]:
        count += 1
    if gencode[kodon]==gencode[nextLatter[nextLatter[kodon[0]]]+kodon[1:3]]:
        count += 1
    if gencode[kodon]==gencode[nextLatter[nextLatter[nextLatter[kodon[0]]]]+kodon[1:3]]:
        count+=1

    return count

def countSismogram(DATA):
    dic={}
    print(len(DATA))
    for i in range(int(len(DATA)/3)):
        dic[i]=countSismogramAt(DATA[3*i:3*i+3])
    return dic

print("in c")
#a    ----->
data=open_GB("corona_2020.gb")
print(data)
dic=countSismogram(data)
print(dic)

