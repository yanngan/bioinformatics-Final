import part_a
import tools


def main_b():
    part_a.main_a()
    print("in b")
    tools.import_data_from_DG_file("BS168.gb")
    tools.import_data_from_UniProt_server("uniprot-organism Bacillus+subtilis+(strain+168)+[224308] .fasta")
