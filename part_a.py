from Bio import SeqIO
from Bio.Seq import Seq
import os
import pandas as pd
from statistics import mean
from matplotlib import pyplot as plt
import re
import tools as t


def checkForNonNucleotide(flag,gene):
    '''
        :param flag: bool flag. indicate if there is non nucleotide in whole DNA
        :param gene: the gene in which we test for non nucleotide
        :return: "yes" if we found non necleotide char and "no" if not
        '''
    if flag:
        pattren = '[^ACGT]'
        if re.findall(pattren, str(gene)) != []:
            return "yes"
    return "no"


def digStartEnd(features,i):
    '''
           :param features: all DNA features
           :param i: specifies the  feature index
           :return: the locations of the gene in DNA.
           '''
    start = features[i].location.start  # calculate position in DNA seq
    end = features[i].location.end
    strand = features[i].location.strand
    return start,end, strand


def checksStart_EndCodon(gene):
    '''
               :param gene: a singal gene from DNA
               :return: string that indicate if start codon or stop codon are missing
               '''
    gene=str(gene)
    start_codon= 'ATG'
    end_codon=['TGA','TAA','TAG']
    e = ""
    if gene[:3]!= start_codon:
        e="missing start codon"
    if gene[len(gene)-3:len(gene)] not in end_codon:
        if e == "":
            return "missing stop codon"
        else:
            return "missing stop and start codon"
    return e


def sanityChecks(features,DNA_seq):

    '''
    function that check for inconsistent value in genbank file and write it to
    'gene_exceptions.csv' file
     :param features: all DNA features
     :param DNA_seq: whole DNA seq
     :return: None
               '''
    pattren = '[^ACGT]'
    nonNucleotide = False
    if re.findall(pattren,str(DNA_seq))==[]:
        nonNucleotide = True
    proteinSanity=[]
    i = 0
    for index in range(1, len(features) - 1):
        if features[index].type == 'gene':
            if features[index + 1].type == 'CDS':
                start, end, strand = digStartEnd(features, index)
                coding_gene = Seq(str(DNA_seq)[start:end])
                if strand ==1:
                    protein = coding_gene.translate(features[index + 1].qualifiers['transl_table'][0],to_stop=True)
                else:
                    reverseGene = Seq(coding_gene.reverse_complement())
                    protein = reverseGene.translate(features[index + 1].qualifiers['transl_table'][0],to_stop=True)


                if not (features[index + 1].qualifiers['translation'][0] == protein):

                    length = len(features[index + 1].qualifiers['translation'][0]) == len(protein)
                    if length == False:
                        proteinSanity.append(
                            [digGeneName(features[index]), 'CDS-gene', start, end, strand, 'yes', 'no',
                             checksStart_EndCodon(coding_gene),checkForNonNucleotide(nonNucleotide,str(coding_gene))])
                    else:
                        proteinSanity.append(
                            [digGeneName(features[index]), 'CDS-gene', start, end, strand, 'no', 'yes',
                             checksStart_EndCodon(coding_gene),checkForNonNucleotide(nonNucleotide,str(coding_gene))])
        i = i+1
    df=createDataFrame(proteinSanity,['name', 'type','start', 'end','strand','inconsistent protein length'
        ,'wrong nucleotide','missing start/end codon','nonNucleotide in gene'])
    print(df)
    df.to_csv("gene_exceptions.csv")


def getAvgCodingGCAndMoreInfo(features,DNA_seq):
    '''
                   :param features: all DNA features
                   :param DNA_seq: whole DNA seq
                   :return: the %GC of coding gene, list of all %GC of every coding gene,
                   data about coding area that will be a dataframe
                   '''
    GC_list = []  # list of all  coding gene %GC
    coding_data = []  # will be dataframe with data about evrey single coding gene
    sum = 0  # length of all coding gene
    gc_sum = 0  # GC counter of all coding gene area
    for index in range(1, len(features) - 1):
        if features[index].type == 'gene' and features[index + 1].type == 'CDS':
            start = features[index].location.start  # calculate position in DNA seq
            end = features[index].location.end
            coding_gene = str(DNA_seq)[start:end-1]

            sum = sum + len(coding_gene)

            temp_gc = str(coding_gene).count('C') + str(coding_gene).count('G')  # count GC of single coding gene
            gc_sum = gc_sum + temp_gc

            coding_data.append([digGeneName(features[index]), str(start), str(end), features[index].location.strand,
                                gcPrcnt(coding_gene)])

            GC_list.append(gcPrcnt(coding_gene))

    return gc_sum/sum, GC_list, coding_data


def gcPrcnt(DNA_seq):
    '''
        :param DNA_seq: whole DNA seq
        :return: the DNA %GC
                       '''
    gc_count = str(DNA_seq).count('C') + str(DNA_seq).count('G')
    return gc_count / len(DNA_seq)


def digGeneName(gene):
    '''
               :param gene: specific gene from DNA
               :return: name of the gene

               '''
    try:
      return  gene.qualifiers['gene'][0]
    except:
        return ""


def flatten(t):
    '''
         :param t: deep list
         :return: flat list
                   '''
    return [item for sublist in t for item in sublist]


def plotHistograms(gene, title, type):
    gene.sort()
    if type=="DNA":
        l = range(0, gene[len(gene) - 1], 100)
        plt.ylabel('length incidence')
        plt.xlabel('gene length')
    else:
        l=50
        plt.ylabel('%GC incidence')
        plt.xlabel('%GC of coding gene')

    plt.hist(gene,bins=l,edgecolor='black')

    plt.title(title)
    plt.tight_layout()
    plt.show()


def getGroupStatisticDataFrame(features,gene_data):
    '''
           :param features: all features from GB file
           :param gene_data: list of list with all  gene lengths
           :return: dataframe with statistic about coding and not coding gene
                          '''
    data_index = 0
    non_coding_DNA = []
    coding_DNA = []
    for index in range(1, len(features) - 1):
        if features[index].type == 'gene':
            if features[index + 1].type == 'CDS':
                coding_DNA.append(gene_data[data_index][0])
                data_index += 1
            else:
                non_coding_DNA.append(gene_data[data_index][0])
                data_index += 1

    non_coding_DNA.sort()
    coding_DNA.sort()
    non_coding_avg = mean(non_coding_DNA)
    coding_avg = mean(coding_DNA)
    length_statistic = [['non coding', non_coding_avg, non_coding_DNA[0], non_coding_DNA[len(non_coding_DNA) - 1]],
                        ['coding', coding_avg, coding_DNA[0], coding_DNA[len(coding_DNA) - 1]]]

    return (createDataFrame(length_statistic, ['group', 'avg', 'min', 'max']),
            non_coding_DNA,coding_DNA)


def readGBFile(path):
    '''
               :param path: path to GB file
               :return: SeqIO object that can be handle by python code
                              '''

    assert (os.path.exists(path))  # Sanity check
    return SeqIO.read(path, "genbank")


def createDataFrame(l,col):
    '''
                   :param l: list of list with data
                   :param col: the  dataframes columns

                   :return:dataframe object
                                  '''
    return pd.DataFrame(l, columns=col)


def countDnaArea(feature_gene):
    '''
          :param feature_gene: all feature of GB file
          :return: dictionary with all types of DNA area as key,
          and count number of each area as value
                                 '''
    dic={}
    for index in range(1, len(feature_gene) - 1):
        if feature_gene[index].type not in dic:
            dic[feature_gene[index].type] = 1
        else:
            dic[feature_gene[index].type] += 1
    return dic


def getAllGeneLength(ftrs):
    '''
              :param ftrs: all features of DNA from GB file
              :return: list of list with all  gene lengths
                             '''
    gene_data = []
    gene_id = -1
    for index in range(1, len(ftrs) - 1):
        if ftrs[index].type == 'gene':
            gene_id += 1
            length = (ftrs[index].location.end - ftrs[index].location.start)  # get the gene length
            gene_data.append([length])
    return gene_data




def main_a():
    print("in a")

    # read gb file and get the features
    file = "BS168.gb"
    record_gb = readGBFile("BS168.gb")
    features = record_gb.features

    t.printSection("a", "1")
    # create dict in form -> area type:number
    dic = countDnaArea(features)
    print(dic)
    t.printSection("a", "2")
    gene_data = getAllGeneLength(features)
    length_df = createDataFrame(gene_data, ['length'])
    print(length_df)


    statistic_length_df, non_coding_DNA, coding_DNA = getGroupStatisticDataFrame(features,gene_data)
    print(statistic_length_df)
    # plot histograms of gene length

    plotHistograms(coding_DNA, 'coding DNA hist',"DNA")
    plotHistograms(non_coding_DNA, 'non coding DNA hist',"DNA")
    plotHistograms(flatten(gene_data), 'all gene hist',"DNA")

    t.printSection("a", "3")

    #  calculate %GC of whole DNA

    DNA_seq = record_gb.seq.upper()
    print("the %GC of whole DNA seq is ", str(gcPrcnt(DNA_seq)*100)+"%")


    # calculate the average %GC of coding DNA and more data

    GC_avg, GC_list, coding_data = getAvgCodingGCAndMoreInfo(features,DNA_seq)

    # show the  five richest and lowest genes in GC

    coding_data = createDataFrame(coding_data, ['name', 'start', 'end', 'strand', 'avg'])
    coding_data.to_csv("part_a.csv")
    d = coding_data.sort_values(by=['avg'])
    print("the  five richest and lowest genes in GC in coding area:")
    print(d)

    print("the %GC over all coding genes is: ",str(GC_avg*100)+"%")

    # plot  GC histogram of coding gene

    plotHistograms(GC_list, 'coding DNA GC hist',"GC")

    t.printSection("a", "4")

    sanityChecks(features, DNA_seq)
