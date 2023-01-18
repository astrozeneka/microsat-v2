import microsatellite
from Bio import SeqIO
from glob import glob
import os

chromosome_orders = ['01', '02', '03', '04', '05', '06', '07', '08', '09',
    '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '90', '91', '92', '93', '94', '95', '96', '97', '98', '99',
                     'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'W', 'X', 'Y', 'Z', 'MIC01', 'MIC02', 'MIC03', 'MIC04', 'MIC05', 'MIC06', 'MIC07', 'MIC08', 'MIC09',
    'MIC10', 'MIC11', 'MIC12', 'MIC13', 'MT']

def search_microsat(keyword, min_repeats, use_cache):
    file_name = f"../data/fasta/{keyword[:3]}/{keyword}.fasta"
    print(f"Search microsatellite in {keyword}")
    record = SeqIO.read(file_name, "fasta")
    tsv_file = f"../data/tmp/{keyword}.tsv"
    if use_cache and os.path.isfile(tsv_file):
        return
    res = microsatellite.searchssr(str(record.seq), min_repeats)
    with open(tsv_file, "w") as f:
        f.write("\n".join(["\t".join([str(a) for a in b]) for b in res]))
    return res
