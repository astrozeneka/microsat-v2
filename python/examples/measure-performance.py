
from src.utils import search_microsat
from time import time

MIN_REPEATS = (12,7,5,4,4,4)

if __name__ == '__main__':
    t = time()
    res = search_microsat("GGA1", MIN_REPEATS, False)
    print(f"Run GGA1 in {int(time() - t)} s")