from Bio import Entrez, SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
import tools
from tabulate import tabulate
from Bio.Align import substitution_matrices
from Bio import Align
import pprint


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
def oneOfTheTwo(a,b,equal):
    if a==equal or b==equal:
        return True
    return False
def getN(kodon):# return how many mutations can be without the start or end kodons
    n=9
    if oneOfTheTwo('_','M',gencode[kodon[0:2]+nextLatter[kodon[2]]]):
        n-=1
    if oneOfTheTwo('_','M',gencode[kodon[0:2]+nextLatter[nextLatter[kodon[2]]]]):
        n-=1
    if oneOfTheTwo('_','M',gencode[kodon[0:2]+nextLatter[nextLatter[nextLatter[kodon[2]]]]]):
        n-=1

    #sec latter
    if oneOfTheTwo('_','M',gencode[kodon[0:1]+nextLatter[kodon[1]]+kodon[2]]):
        n-=1
    if oneOfTheTwo('_','M',gencode[kodon[0:1]+nextLatter[nextLatter[kodon[1]]]+kodon[2]]):
        n-=1
    if oneOfTheTwo('_','M',gencode[kodon[0:1]+nextLatter[nextLatter[nextLatter[kodon[1]]]]+kodon[2]]):
        n-=1

    #third latter
    if oneOfTheTwo('_','M',gencode[nextLatter[kodon[0]]+kodon[1:3]]):
        n-=1
    if oneOfTheTwo('_','M',gencode[nextLatter[nextLatter[kodon[0]]]+kodon[1:3]]):
        n-=1
    if oneOfTheTwo('_','M',gencode[nextLatter[nextLatter[nextLatter[kodon[0]]]]+kodon[1:3]]):
        n-=1

    return n
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


    return count*3/getN(kodon)

def countSynonyms():# return fs*3/n => synonyms count
    dic={}
    kodon="AAA"
    for i in range(4):
        for j in range(4):
            for k in range(4):
                dic[kodon]=countSynonymsAt(kodon)
                kodon=kodon[0:2]+nextLatter[kodon[2]]
            kodon=kodon[0] +nextLatter[kodon[1]]+kodon[2]
        kodon= nextLatter[kodon[0]]+kodon[1:3]
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
def correctKodon(seq,gene):
    fixed_seq=""
    seq_counter=0
    for i in range(len(gene)):
        if gene[i]=='-':
            fixed_seq+='---'
        else:
            fixed_seq+=seq[seq_counter:seq_counter+3]
            seq_counter+=3
    return fixed_seq
    pass
def calculateDNDS(seqA,seqB,geneA,geneB):# return dn , ds , and dnds ratio
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

    alignments = aligner.align(geneA, geneB)
    string=str(alignments[0])


    index=string.find('\n')
    in_order_a=correctKodon(seqA,string[0:index])
    in_order_b=correctKodon(seqB,string[index*2+2:-1])
    dN, dS=cal_dn_ds(CodonSeq(in_order_a[0:min(len(in_order_a),len(in_order_b))]),CodonSeq(in_order_b[0:min(len(in_order_a),len(in_order_b))]))
    dN_dS_ratio=-1
    if dS!=0:
        dN_dS_ratio= float(dN/dS)
    return dN,dS,dN_dS_ratio
    pass

def getKODON(gene,data):# return the kodons that made the gene
    end=gene[1].end
    start=gene[1].start
    return data[start:end]
    pass

def getGeneNames(gene):# return names of the genes (remove dups as well)
    names=[]
    for i in range(len(gene)):
        gene_name=gene[i][2].qualifiers['gene'][0]
        if gene_name not in names:
            names.append(gene_name)
    return names
def compareOrganizm(genA,genB):# count how many proteins are equal
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
def main_c():
    print()

    #load files
    data_corona=open_GB("corona_2020.gb")
    data_corona_recent=open_GB("corona_2022.gb")



    #     A      ----->
    dic=countSynonyms()
    count =0
    for sub_dict  in dic.keys():
        spaces=""
        num=str(round(float(dic[sub_dict]), 3))
        for i in range(5-len(num)):
            spaces=spaces+" "
        print(sub_dict+":"+num,end=spaces+"       ")
        count+=1
        if(count==4):
            count=0
            print()







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
            print("both have the same genes, some proteins are different")
        else:
            print("both have the same exactly gene's")
    else:
        print("the amount of unequal genes:",unequalSum)



    #     B.b    ----->
    print()
    dndsList=[]
    toPrint=[[]*5]
    F=2
    T=7
    for i in range(F,T):
        dndsList.append(calculateDNDS(getKODON(genes_2020[i],get_sequence(data_corona)),getKODON(genes_2022[i],get_sequence(data_corona_recent)),genes_2020[i][0],genes_2022[i][0]))
        dN, dS, dN_dS_ratio = dndsList[i-F]
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
        toPrint.append([gene[0],locus_tag[0],protein_id[0],product[0],Type,dN,dS,dN_dS_ratio,selectionType])


    print(tabulate(toPrint, headers=['gene', 'locus_tag','protein_id','product','type','dN','dS','dN_dS_ratio','prefences'], tablefmt='orgtbl'))# print the information table


main_c()
