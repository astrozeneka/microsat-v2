
USE_CACHE = True
MIN_REPEATS = (12,7,5,4,4,4)

import microsatellite
from Bio import SeqIO
from glob import glob
import os
from Bio.Seq import Seq


import argparse
parser = argparse.ArgumentParser(description='MicrosatNavigator: Analyze microsatellite')
parser.add_argument('--fasta', help='FASTA input file')
parser.add_argument('--output', help='Output file')
parser.add_argument('--fastq', help='FASTQ input file')
parser.add_argument('--standardize', default=True)
args = parser.parse_args()

_STANDARDS={}
def get_motif_standard(motif):
    if motif == "":
        return ""
    if motif in _STANDARDS.keys():
        return _STANDARDS[motif]
    options = []
    for i in range(len(motif)):
        ex = motif[i:] + motif[:i]
        options.append(ex)
        options.append(str(Seq(ex).complement()))
    def cost(motif):
        motif = motif.upper()
        output = 0
        for i in range(len(motif)):
            output += 2**(len(motif)-i-1) * "ACTGNMKWRY".index(motif[i])
        return (motif, output)
    cost_list = [cost(m) for m in options]
    standard = min(cost_list)[0]
    _STANDARDS[motif] = standard
    return standard


def search_microsat_fastq(keyword, min_repeats):
    print(f"Search microsatellites in {keyword}")
    output = []
    for record in SeqIO.parse(args.fastq, 'fastq'):
        res = microsatellite.searchssr(str(record.seq), min_repeats)
        res = [[record.name] + list(a) for a in res]
        if args.standardize:
            for i, row in enumerate(res):
                res[i][1] = get_motif_standard(row[1])
        output = output+res
    with open(args.output, "w") as f:
        f.write("\n".join(["\t".join([str(a) for a in b]) for b in output]))
    print("Done")

def search_microsat(keyword, min_repeats, use_cache):
    with open(args.output, "w") as f:
        f.write("")
    print(f"Search microsatellites in {keyword}")
    counter = 0
    for record in SeqIO.parse(args.fasta, 'fasta'):
        res = microsatellite.searchssr(str(record.seq), min_repeats)
        res = [[record.name] + list(a) for a in res]
        if args.standardize:
            for i, row in enumerate(res):
                res[i][1] = get_motif_standard(row[1])
        tsv = "\n".join(['\t'.join([str(a) for a in b]) for b in res])
        tsv += "\n" if len(tsv) > 0 else ""
        with open(args.output, 'a') as f:
            f.write(tsv)
        counter += 1
        if counter % 1000000 == 0:
            print(f"{int(counter/1000000)}M reads ...")
    print("Analysis done")

def genome_chromosome_list(genome):
    output = glob(f"../data/fasta/{genome}/*.fasta")
    print()

if __name__ == '__main__':

    # If not using cache, clear it
    if not USE_CACHE:
        tmp_list = glob("../data/tmp/*")
        for fn in tmp_list:
            os.unlink(fn)

    # chr_list = genome_chromosome_list()

    #rec = SeqIO.read("../data/LAG01.fasta", "fasta")
    if args.fastq:
        search_microsat_fastq(os.path.basename(args.fastq), MIN_REPEATS)
    elif args.fasta:
        search_microsat(os.path.basename(args.fasta), MIN_REPEATS, USE_CACHE)
    else:
        print("No input file specified")
    print()