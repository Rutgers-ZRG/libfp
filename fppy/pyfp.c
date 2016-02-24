/*******************************************************************************
 * Copyright (C) 2015 Li Zhu 
 * All rights reserved. 
 * 
 * pyfp.c
 * This file is part of fplib.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 * ****************************************************************************/

#include <Python.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
#include "../src/fplib.h"

static PyObject * get_fppy_periodic(PyObject *self, PyObject *args);
static PyObject * get_fppy_dist_periodic(PyObject *self, PyObject *args);

static PyMethodDef functions[] = {
    {"fp_periodic", get_fppy_periodic, METH_VARARGS, "get fingerprint" },
    {"fp_dist", get_fppy_dist_periodic, METH_VARARGS, "get fingerprint dist"},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initfppy(void)
{
    Py_InitModule3("fppy", functions, "C-extension for fplib\n\n...\n");
    return ;
}

static PyObject * get_fppy_periodic(PyObject *self, PyObject *args)
{
    int i, j, lmax, l, lseg, natx;
    double cutoff;
    double **sfp, **lfp;
    PyArrayObject* lattice;
    PyArrayObject* position;
    PyArrayObject* atom_type;
    PyArrayObject* znu;
    PyObject* array, *vec, *shortfp, *longfp;

    if (!PyArg_ParseTuple(args, "OOOOiid", &lattice, &position, &atom_type, &znu, 
                &lmax, &natx, &cutoff))
        return NULL;

    double (*lat)[3] = (double(*)[3])lattice->data;
    double (*rxyz)[3] = (double(*)[3])position->data;
    int nat = position->dimensions[0];
    int ntyp = znu->dimensions[0];
    long* znuclong = (long*)znu->data;
    long* typeslong = (long*)atom_type->data;

    int types[nat];
    for ( i = 0; i < nat; i++)
        types[i] = (int)typeslong[i];

    int znucl[ntyp];
    for ( i = 0; i < ntyp; i++)
        znucl[i] = (int)znuclong[i];

    if (lmax == 0) {
        lseg = 1; l = 1;
    } else if (lmax == 1) {
        lseg = 4; l = 2;
    } else {
        fprintf(stderr, "ORBITAL ERROR."); return NULL;
    }


    sfp = (double **) malloc(sizeof(double)*nat);
    lfp = (double **) malloc(sizeof(double)*nat);

    for ( i = 0; i < nat; i++ ) {
        sfp[i] = (double *) malloc(sizeof(double) * l * (ntyp + 1));
        lfp[i] = (double *) malloc(sizeof(double) * (natx * lseg));
    }

    get_fp_periodic(lmax, nat, ntyp, types, lat, rxyz, znucl, natx,  cutoff, sfp, lfp);
    shortfp = PyList_New(0);
    longfp = PyList_New(0);
    for ( i = 0; i < nat; i++ ) {
        vec = PyList_New(0);
        for ( j = 0; j < l*(ntyp+1); j++){
            PyList_Append(vec, PyFloat_FromDouble( sfp[i][j] ) );
        }
        PyList_Append( shortfp, vec );
        vec = PyList_New(0);
        for (j = 0; j < natx*lseg; j++) {
            PyList_Append(vec, PyFloat_FromDouble( lfp[i][j] ) );
        }
        PyList_Append( longfp, vec);
    }
    
    array = PyList_New(0);
    PyList_Append( array, shortfp );
    PyList_Append( array, longfp );

    for (i = 0; i < nat; i++) {
        free(sfp[i]);
        free(lfp[i]);
    }
    free(sfp);
    free(lfp);

    return array;

}


static PyObject * get_fppy_dist_periodic(PyObject *self, PyObject *args)
{
    int i, j, ntyp;
    double fpdist;
    PyArrayObject* tfp1;
    PyArrayObject* tfp2;
    PyArrayObject* atom_type;

    if (!PyArg_ParseTuple(args, "iOOO", &ntyp, &atom_type, &tfp1, &tfp2))
        return NULL;

    // double (*
    int nat = tfp1->dimensions[0];
    int fp_len = tfp1->dimensions[1];
    long* typeslong = (long*)atom_type->data;
    int types[nat];
    for ( i = 0; i < nat; i++) {
        types[i] = (int)typeslong[i];
    }

    double (*fptt1)[fp_len] = (double(*)[fp_len])tfp1->data;
    double (*fptt2)[fp_len] = (double(*)[fp_len])tfp2->data;

    double **fp1;
    double **fp2;

    fp1 = (double **) malloc(sizeof(double)*nat);
    fp2 = (double **) malloc(sizeof(double)*nat);
    for ( i = 0; i < nat; i++ ) {
        fp1[i] = (double *) malloc(sizeof(double) * fp_len);
        fp2[i] = (double *) malloc(sizeof(double) * fp_len);
    }


    for (i = 0; i < nat; i++){
        for (j = 0; j < fp_len; j++){
            fp1[i][j] = fptt1[i][j];
            fp2[i][j] = fptt2[i][j];
        }
    }


    fpdist = get_fpdistance_periodic(nat, ntyp, types, fp_len, fp1, fp2);

    free(fp1);
    free(fp2);

    // printf("%d %d \n ", nat, nat2);
    return PyFloat_FromDouble(fpdist);

}







