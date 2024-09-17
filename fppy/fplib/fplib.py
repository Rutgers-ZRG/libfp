#/*******************************************************************************
# * Copyright (C) 2021 Zhu Research Group @ Rutgers-Newark
# * All rights reserved.
# *
# * This file is part of fplib.
# *
# * Permission is hereby granted, free of charge, to any person obtaining a copy
# * of this software and associated documentation files (the "Software"), to deal
# * in the Software without restriction, including without limitation the rights
# * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# * copies of the Software, and to permit persons to whom the Software is
# * furnished to do so, subject to the following conditions:
# *
# * The above copyright notice and this permission notice shall be included in
# * all copies or substantial portions of the Software.
# *
# * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# * THE SOFTWARE.
# * ****************************************************************************/

from . import _fplib as fp
import numpy as np

def get_version():
    return tuple(fp.version())


def get_lfp(cell,
            cutoff=4.0,
            log=True,
            orbital='s',
            natx=300):
    '''
    cell : tuple (lattice, rxyz, types, znucl)
    '''
    (lat, rxyz, types, znucl) = _expand_cell(cell)
    if log is True:
        wlog = 1
    else:
        wlog = 0
    if orbital == 's':
        lmax = 0
    elif orbital == 'sp':
        lmax = 1
    else:
        print ('Warning: wrong type of orbital')
        lmax = 0

    (sfp, lfp) = fp.fp_periodic(0, 0, wlog, lat, rxyz, types, znucl, lmax, natx, cutoff)
    return np.array(lfp)


def get_sfp(cell,
            cutoff=4.0,
            log=True,
            orbital='s',
            natx=300):
    '''
    cell : tuple (lattice, rxyz, types, znucl)
    '''
    lat, rxyz, types, znucl = _expand_cell(cell)
    if log is True:
        wlog = 1
    else:
        wlog = 0
    if orbital == 's':
        lmax = 0
    elif orbital == 'sp':
        lmax = 1
    else:
        print ('Warning: wrong type of orbital')
        lmax = 0

    (sfp, lfp) = fp.fp_periodic(1, 0, wlog, lat, rxyz, types, znucl, lmax, natx, cutoff)
    return np.array(sfp)

def get_dfp(cell,
            cutoff=4.0,
            log=True,
            orbital='s',
            natx=300):
    '''
    cell : tuple (lattice, rxyz, types, znucl)
    '''
    (lat, rxyz, types, znucl) = _expand_cell(cell)
    if log is True:
        wlog = 1
    else:
        wlog = 0
    if orbital == 's':
        lmax = 0
    elif orbital == 'sp':
        lmax = 1
    else:
        print ('Warning: wrong type of orbital')
        lmax = 0

    (sfp, lfp, dfp) = fp.fp_periodic(0, 1, wlog, lat, rxyz, types, znucl, lmax, natx, cutoff)
    return np.array(lfp), np.array(dfp)


def get_fp_dist(lfp1, lfp2, types, assignment=False):
    ntyp = len(set(types))
    dist, m = fp.fp_dist(ntyp, types, lfp1, lfp2)
    if assignment:
        return dist, m
    else:
        return dist

def get_nfp(rxyz, types, znucl):
    nfp = fp.fp_nonperiodic(rxyz, types, znucl)
    return np.array(nfp)

def get_nfp_dist(fp1, fp2):
    d = fp1 - fp2
    return np.sqrt(np.dot(d, d))

def _expand_cell(cell):
    lat = np.array(cell[0], float)
    rxyz = np.array(cell[1], float)
    types = np.array(cell[2], int)
    znucl = np.array(cell[3], int)
    if len(rxyz) != len(types) or len(set(types)) != len(znucl):
        raise ValueError('Something wrong with rxyz / types / znucl.')
    return lat, rxyz, types, znucl

