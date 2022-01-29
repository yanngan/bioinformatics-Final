from Bio import Entrez, SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
import tools
from tabulate import tabulate

#--------------------------------------global values-------------------------------------#
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


#.........................................................................................#
#.........................................................................................#
#-----------------------------------------------------------------------------------------#



#--------------------------------------load data------------------------------------------#
def open_GB( GB_file_path):
    return tools.import_data_from_DB_file(GB_file_path)

def get_sequence(GB_file_path):
    return  GB_file_path[0].sequence

def get_protain(GB_file_path):
    GOA_dictionary_to_return,Name_dictionary_to_return,all_protein=tools.get_dictionary_of_protain_from_GB_file(GB_file_path)
    return all_protein
#------------------------------------------A----------------------------------------------#

def countSynonymsAt(kodon):
    # count a certain unit how many synonyms he has
    count=0
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

def countSynonyms(DATA):
    dic={}
    print(len(DATA))
    for i in range(int(len(DATA)/3)):
        # check each unit how many Synonyms he has
        dic[i]=countSynonymsAt(DATA[3*i:3*i+3])
    return dic

#.........................................................................................#
#.........................................................................................#
#------------------------------------------B.a----------------------------------------------#
def count_equal_genes(genesA,genesB):
    count=0
    for i in range(len(genesA)):
        if genesA[i][0]==genesB[i][0]:
            count=+1
    return count


#.........................................................................................#
#------------------------------------------B.b----------------------------------------------#
def calculateDNDS(seqA,seqB):
    dN, dS=cal_dn_ds(CodonSeq(seqA),CodonSeq(seqB))
    dN_dS_ratio=-1
    if dS!=0:
        dN_dS_ratio= float(dN/dS)
    return dN,dS,dN_dS_ratio
    pass

def getKODON(gene,data):
    end=gene[1].end
    start=gene[1].start
    return data[start:end]
    pass

def getGeneNames(gene):
    names=[]
    for i in range(len(gene)):
        gene_name=gene[i][2].qualifiers['gene'][0]
        if gene_name not in names:
            names.append(gene_name)
    return names
def compareOrganizm(genA,genB):
    equalGenes=0
    unqualGenes=0
    gen_names_A=getGeneNames(genA)
    gen_names_B = getGeneNames(genB)
    for i in range(min(len(gen_names_A),len(gen_names_B))):
        if gen_names_A[i]==gen_names_B[i]:
            equalGenes+=1
        else:
            unqualGenes+=1


    return equalGenes,unqualGenes

    pass

#.........................................................................................#
#.........................................................................................#
#.........................................................................................#
#-----------------------------------------------------------------------------------------#

#load files
data_corona=open_GB("corona_2020.gb")
data_corona_recent=open_GB("corona_2022.gb")




#     A      ----->
dic=countSynonyms(get_sequence(data_corona))
print(dic)




#     B.a    ----->
print()
genes_2020=get_protain("corona_2020.gb")
genes_2022=get_protain("corona_2022.gb")

count_of_equal_genes=count_equal_genes(genes_2020,genes_2022)
print("exacaly same genes:",count_of_equal_genes)
equalSum,unequalSum=compareOrganizm(genes_2020,genes_2022)
print("same genes overall:",equalSum)
if unequalSum==0:
    if equalSum!=count_of_equal_genes:
        print("both have the same genes, some protains are diffrent")
    else:
        print("both have the same exacaly gene's")
else:
    print("the amount of unequal genes:",unequalSum)



#     B.b    ----->
print()
dndsList=[]
toPrint=[[]*5]
for i in range(3,8):
    dndsList.append(calculateDNDS(getKODON(genes_2020[i],get_sequence(data_corona))[:-3],getKODON(genes_2022[i],get_sequence(data_corona_recent))[:-3]))
    dN, dS, dN_dS_ratio = dndsList[i-3]
    gene=genes_2020[i][2].qualifiers['gene']
    locus_tag=genes_2020[i][2].qualifiers['locus_tag']
    protein_id=genes_2020[i][2].qualifiers['protein_id']
    product=genes_2020[i][2].qualifiers['product']
    Type=genes_2020[i][2].type
    selectionType=""

    if dN_dS_ratio>1.25 or dN_dS_ratio==-1:
        selectionType="Positive selection"
    elif dN_dS_ratio<0.75:
        selectionType="Negative selection"
    else:
        selectionType="Neutral selection"
    if dN_dS_ratio==-1:
        dN_dS_ratio="unvalid"
    toPrint.append([gene[0],locus_tag[0],protein_id[0],product,Type,dN,dS,dN_dS_ratio,selectionType])


print(tabulate(toPrint, headers=['gene', 'locus_tag','protein_id','product','type','dN','dS','dN_dS_ratio','prefences'], tablefmt='orgtbl'))



