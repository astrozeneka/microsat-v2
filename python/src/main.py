import microsatellite

if __name__ == '__main__':
    sequence = "CCTGGACCT"
    res = microsatellite.searchssr(sequence, (4, 5, 7))
    print()