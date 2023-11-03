GENOMES = ["LAG"]
USE_CACHE = True
MIN_REPEATS = (12,7,5,4,4,4)

import microsatellite
from Bio import SeqIO
from glob import glob
import os


import argparse
parser = argparse.ArgumentParser(description='MicrosatNavigator: Analyze microsatellite')
parser.add_argument('--fasta', help='Fisher association data')
parser.add_argument('--output', help='Output file')
args = parser.parse_args()

def search_microsat(keyword, min_repeats, use_cache):
    file_name = args.fasta
    print(f"Search microsatellites in {keyword}")
    record = SeqIO.read(file_name, "fasta")
    tsv_file = f"../data/tmp/{keyword}.tsv"
    if use_cache and os.path.isfile(tsv_file):
        return
    res = microsatellite.searchssr(str(record.seq), min_repeats)
    tsv_file = args.output
    #tsv_file = args.fasta.replace(".fasta", ".tsv").replace(".fna", ".tsv").replace(".fa", ".tsv")
    with open(tsv_file, "w") as f:
        f.write("\n".join(["\t".join([str(a) for a in b]) for b in res]))
    return res

def genome_chromosome_list(genome):
    output = glob(f"../data/fasta/{genome}/*.fasta")
    print()

if __name__ == '__main__':

    # If not using cache, clear it
    if not USE_CACHE:
        tmp_list = glob("../data/tmp/*")
        for fn in tmp_list:
            os.unlink(fn)

    GENOME = GENOMES[0]
    chr_list = genome_chromosome_list(GENOME)

    #rec = SeqIO.read("../data/LAG01.fasta", "fasta")
    search_microsat("LAG1", MIN_REPEATS, USE_CACHE)
    print()