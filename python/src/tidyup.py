from glob import glob
from Bio import SeqIO
import shutil

name_symbol_list = [
    ("pogoniulus pusillus", "PPU"),
    ("gallus gallus", "GGA"),
    ("lacerta agilis", "LAG"),
    ("corvus hawaiiensis", "CHA"),
    ("rhynochetos jubatus", "RJU"),
    ("crotalus viridis", "CVI"),
    ("naja naja", "NNA"),
    ("anolis carolinensis", "ACA"),
    ("falco naumanni", "FNA"),
    ("clarias gariepinus", "CGA"),
    ("homo sapiens", "HSA")
]
chr_name_symbol_list = [
    f"chromosome {a.lower()}"
    for a in ['1', '2', '3', '4', '5', '6', '7', '8', '9',
    '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '90', '91', '92', '93', '94', '95', '96', '97', '98', '99',
                      'W', 'X', 'Y', 'Z', 'MIC_1', 'MIC_2', 'MIC_3', 'MIC_4', 'MIC_5', 'MIC_6', 'MIC_7', 'MIC_8', 'MIC_9',
    'MIC_10', 'MIC_11', 'MIC_12', 'MIC_13', 'MT']
]
chr_name_symbol_list = chr_name_symbol_list+[
    f"linkage group {a.lower()}"
    for a in ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']
]
chr_name_symbol_list.reverse()
chr_name_symbol_list = chr_name_symbol_list+[" mitochondrion"]

if __name__ == '__main__':
    # tidy up data workspace
    file_list = glob("../data/tray/*")
    for file_name in file_list:
        rec = SeqIO.read(file_name, "fasta")

        # The symbol
        symbol = [a[1] for a in name_symbol_list if a[0] in rec.description.lower()]
        if len(symbol) == 0:
            raise "Error"
        symbol = symbol[0]

        chr_symbol = [a for a in chr_name_symbol_list if a in rec.description.lower()]
        if len(chr_symbol) == 0:
            raise "Error"
        chr_symbol = chr_symbol[0]
        chr_symbol = chr_symbol[chr_symbol.index(" ")+1:]
        chr_symbol = chr_symbol.replace("mitochondrion", "mt")
        chr_symbol = chr_symbol.replace("_", "")
        destination_name = symbol+chr_symbol+".fasta"
        shutil.move(file_name, f"../data/fasta/{symbol}/{destination_name}")
        print(f"Move {file_name} to ../data/fasta/{destination_name}")