from Bio import Entrez, SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds



# def NewCode(currentCode, changeAT):
#     nextlatter = nextLatter[currentCode[changeAT]]
#     currentCode = currentCode[0:changeAT] + nextlatter + currentCode[changeAT + 1:3]
#     return currentCode
#     pass


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



def addGaps(seq1, seq2):
    alignment = pairwise2.align.globalxx(seq1, seq2)[0]
    return alignment.seqA,alignment.seqB,alignment.score
    pass


def OnlyInFirst(gen,seq1,seq2):
    newGens={}
    for i in range(len(seq1)):
        if seq2[i]=='-':
            newGens[i]=seq1[i]
    return newGens
    pass
def compare(first_RNA,sec_RNA):
    pass


def getNumberOfEqual(howMany,seqA,seqB):
    dic={}
    for i in range(len(seqA)):
        if howMany==0:
            i=len(seqA)
            break
        if seqA[i]==seqB[i]:
            dic[i]=seqA[i]
            howMany-=1
    return dic
def calculateDNDS(seqA,seqB):
    dN, dS=cal_dn_ds(seqA,seqB)
    dN_dS_ratio= float(dN,dS)
    return dN,dS,dN_dS_ratio
    pass
def findGenFrom(start,seq):
    index=findAUG(seq[start:])
    pass

# def main_b():
#     print("in c")
#     #a    ----->
#     data=open_GB("corona_2020.gb")
#     dic=countSismogram(data)
#     print(dic)
#     #b    ----->
#
data_corona=open_GB("corona_2020.gb")
print(data_corona)
dic=countSismogram(data_corona)
print(dic)

data_corona_recent=open_GB("corona_2022.gb")
TESTA="ACG"
TESTB="AAG"

corona_2020,corona_2022,score=addGaps(TESTA,TESTB)
print(corona_2020,corona_2022,score)
corona_2020_old_gens=OnlyInFirst(corona_2020,corona_2022)
corona_2022_new_gens=OnlyInFirst(corona_2022,corona_2020)
print(corona_2020_old_gens)
print(corona_2022_new_gens)
equal_gens=getNumberOfEqual(5,corona_2020,corona_2022)
gens={}
From_2020=0
From_2022=0
for i in range(len(corona_2020)):
    genA,From_2020=findGenFrom(From_2020,corona_2020)
    genB,From_2022=findGenFrom(From_2022,corona_2020)
    gens[i]=[genA,genB]


print(equal_gens)
dnds_dic={}
for i in range(gens):
    dN,dS,dN_dS_ratio=calculateDNDS(gens[i][0],gens[i][1])
    dnds_dic[i]=[dN,dS,dN_dS_ratio]

matches,additional_gens=compare(corona_2020,corona_2022)