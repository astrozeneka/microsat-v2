GENOMES = ["LAG"]
USE_CACHE = True
MIN_REPEATS = (12,7,5,4,4,4)

import microsatellite
from Bio import SeqIO
from glob import glob
import os

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



if __name__ == '__main__':

    # If not using cache, clear it
    if not USE_CACHE:
        tmp_list = glob("../data/tmp/*")
        for fn in tmp_list:
            os.unlink(fn)

    rec = SeqIO.read("../data/LAG01.fasta", "fasta")
    search_microsat("LAG1", MIN_REPEATS, USE_CACHE)
    # res = microsatellite.searchssr(str(rec.seq), min_repeats)