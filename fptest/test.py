#!/usr/bin/python

import numpy as np
import fppy


def readvasp(poscar):
    buff = []
    with open(poscar) as f:
        for line in f:
            buff.append(line.split())

    lat = np.array(buff[2:5], float)
    try:
        typt = np.array(buff[5], int)
    except:
        del(buff[5])
        typt = np.array(buff[5], int)
    nat = sum(typt)
    pos = np.array(buff[7:7+nat], float)
    rxyz = np.dot(pos, lat)
    types = []
    for i in range(len(typt)):
        types += [i + 1] * typt[i]
    types = np.array(types, int)

    return (lat, rxyz, types)


def main():
    lmax = 0
    natx = 300
    znucl = [20, 5, 6]
    znucl = np.array(znucl)
    cutoff = 4

    ntyp = len(znucl)
    (lat1, rxyz1, types) = readvasp('pos1')
    (lat2, rxyz2, types) = readvasp('pos2')
    (sfp1, lfp1) = fppy.fp_periodic(lat1, rxyz1, types, znucl, lmax, natx, cutoff)
    (sfp2, lfp2) = fppy.fp_periodic(lat2, rxyz2, types, znucl, lmax, natx, cutoff)
    fp1 = np.array(sfp1)
    fp2 = np.array(sfp2)

    (dist, m) = fppy.fp_dist(ntyp, types, fp1, fp2)
    print 'dist', dist
    print m


if __name__ == '__main__':
    main()



