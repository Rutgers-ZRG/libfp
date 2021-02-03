#!/usr/bin/python env python3

import numpy as np
import fplib 

def readvasp(vp):
    buff = []
    with open(vp) as f:
        for line in f:
            buff.append(line.split())

    lat = np.array(buff[2:5], float)
    try:
        typt = np.array(buff[5], int)
    except:
        del(buff[5])
        typt = np.array(buff[5], int)
    nat = sum(typt)
    pos = np.array(buff[7:7 + nat], float)
    types = []
    for i in range(len(typt)):
        types += [i+1]*typt[i]
    types = np.array(types, int)
    rxyz = np.dot(pos, lat)
    return lat, rxyz, types

def test():
    znucl = [14, 8]
    lat1, rxyz1, types1 = readvasp('struct1.vasp')
    lat2, rxyz2, types2 = readvasp('struct2.vasp')

    cell1 = (lat1, rxyz1, types1, znucl)
    cell2 = (lat2, rxyz2, types2, znucl)
    fp1 = fplib.get_lfp(cell1, cutoff=6, log=True)
    fp2 = fplib.get_lfp(cell2, cutoff=6, log=True)

    dist = fplib.get_fp_dist(fp1, fp2, types1)

    print ('fingerprint distance between struct1 and struct2: ', dist)

if __name__ == "__main__":
    test()
