from Bio import SeqIO
import os
import pandas as pd
from statistics import mean
from matplotlib import pyplot as plt


def getAvgCodingGCAndMoreInfo(features,DNA_seq):
    GC_list = []  # list of all  coding gene %GC
    coding_data = []  # will be dataframe with data about evrey single coding gene
    sum = 0  # length of all coding gene
    gc_sum = 0  # GC counter of all coding gene area
    for index in range(1, len(features) - 1):
        if features[index].type == 'gene' and features[index + 1].type == 'CDS':
            start = features[index].location.start  # calculate position in DNA seq
            end = features[index].location.end
            coding_gene = str(DNA_seq)[start:end]

            sum = sum + len(coding_gene)

            temp_gc = str(coding_gene).count('C') + str(coding_gene).count('G')  # count GC of single coding gene
            gc_sum = gc_sum + temp_gc

            coding_data.append([digGeneName(features[index]), str(start), str(end), features[index].location.strand,
                                gcPrcnt(coding_gene)])

            GC_list.append(gcPrcnt(coding_gene))

    return gc_sum/sum, GC_list, coding_data


def gcPrcnt(DNA_seq):
    gc_count = str(DNA_seq).count('C') + str(DNA_seq).count('G')
    return gc_count / len(DNA_seq)


def digGeneName(gene):
    try:
      return  gene.qualifiers['gene'][0]
    except:
        print("missing gene name")
        return ""


def flatten(t):
    return [item for sublist in t for item in sublist]


def plotHistograms1(gene, title):
    gene.sort()

    plt.hist(gene,bins=50,edgecolor='black')
    plt.ylabel('%GC incidence')
    plt.xlabel('%GC of coding gene')
    plt.title(title)
    plt.tight_layout()
    plt.show()


def plotHistograms(gene, title):
    gene.sort()
    l=range(0,gene[len(gene)-1],100)
    plt.hist(gene,bins=l,edgecolor='black')
    plt.ylabel('length incidence')
    plt.xlabel('gene length')
    plt.title(title)
    plt.tight_layout()
    plt.show()


def getGroupStatisticDataFrame(features,gene_data):
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
    assert (os.path.exists(path))  # Sanity check
    return SeqIO.read(path, "genbank")


def createDataFrame(l,col):
    return pd.DataFrame(l, columns=col)


def countDnaArea(feature_gene):
    dic={}
    for index in range(1, len(feature_gene) - 1):
        if feature_gene[index].type not in dic:
            dic[feature_gene[index].type] = 1
        else:
            dic[feature_gene[index].type] += 1
    return dic


def getAllGeneLength(ftrs):
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

    # create dict in form -> area type:number
    dic = countDnaArea(features)

    gene_data = getAllGeneLength(features)
    length_df = createDataFrame(gene_data, ['length'])

    statistic_length_df, non_coding_DNA, coding_DNA = getGroupStatisticDataFrame(features,gene_data)

    # plot histograms of gene length

    plotHistograms(coding_DNA, 'coding DNA hist')
    plotHistograms(non_coding_DNA, 'non coding DNA hist')
    plotHistograms(flatten(gene_data), 'all gene hist')

    #  calculate %GC of whole DNA

    DNA_seq = record_gb.seq.upper()
    print("the %GC of whole DNA seq is ", gcPrcnt(DNA_seq))

    # calculate the average %GC of coding DNA and more data

    GC_avg, GC_list, coding_data = getAvgCodingGCAndMoreInfo(features,DNA_seq)

    # show the  five richest and lowest genes in GC

    coding_data = createDataFrame(coding_data, ['name', 'start', 'end', 'strand', 'avg'])
    d = coding_data.sort_values(by=['avg'])
    print(d)

    print(GC_avg)

    # plot  GC histogram of coding gene

    plotHistograms1(GC_list, 'coding DNA GC hist')

main_a()