/*******************************************************************************
 * Copyright (C) 2021 Zhu Research Group @ Rutgers-Newark
 * All rights reserved.
 *
 * pyfp.c
 * This file is part of libfp.
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
#include <assert.h>
#include <numpy/arrayobject.h>
#include <libfp.h>


static PyObject * py_get_version(PyObject *self, PyObject *args); 
static PyObject * py_get_nonperiodic(PyObject *self, PyObject *args);
static PyObject * py_get_periodic(PyObject *self, PyObject *args);
static PyObject * py_get_dist_periodic(PyObject *self, PyObject *args);

struct module_state
{
    PyObject *error;
};
static struct module_state _state;
static PyObject *
error_out(PyObject *m) {
  struct module_state *st = (struct module_state*)PyModule_GetState(m);
  PyErr_SetString(st->error, "something bad happened");
  return NULL;
}

static PyMethodDef _libfp_methods[] = {
    {"version", py_get_version, METH_VARARGS, "Fplib version"},
    {"fp_nonperiodic", py_get_nonperiodic, METH_VARARGS, "get fingerprint non" },
    {"fp_periodic", py_get_periodic, METH_VARARGS, "get fingerprint" },
    {"fp_dist", py_get_dist_periodic, METH_VARARGS, "get fingerprint dist"},
    {NULL, NULL, 0, NULL}
};

#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))

static int _libfp_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int _libfp_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_libfp",
  NULL,
  sizeof(struct module_state),
  _libfp_methods,
  NULL,
  _libfp_traverse,
  _libfp_clear,
  NULL
};

#define INITERROR return NULL

PyObject *
PyInit__libfp(void)
{
  struct module_state *st;
  PyObject *module = PyModule_Create(&moduledef);

  if (module == NULL)
    INITERROR;

  st = GETSTATE(module);

  st->error = PyErr_NewException("_libfp.Error", NULL, NULL);
  if (st->error == NULL) {
    Py_DECREF(module);
    INITERROR;
  }
  return module;
}

static PyObject * py_get_version(PyObject *self, PyObject *args)
{   
    PyObject *array;
    int i;
    int version[3];

    version[0] = get_major_version();
    version[1] = get_minor_version();
    version[2] = get_micro_version();
    array = PyList_New(3);
    for (i = 0; i < 3; i++) {
        PyList_SetItem(array, i, PyLong_FromLong((long)version[i]));
    }
    return array;
}


static PyObject * py_get_nonperiodic(PyObject *self, PyObject *args)
{
    int i, lseg, nid;
    double *fp;
    PyArrayObject* py_positions;
    PyArrayObject* py_atom_type;
    PyArrayObject* py_znu;
    PyObject* vec;
    int nat, ntyp;
    PyObject* pytmp;


    if (!PyArg_ParseTuple(args, "OOO", &py_positions, &py_atom_type, &py_znu))
        return NULL;

    long* pznuc = (long*)PyArray_DATA(py_znu);
    long* ptypes = (long*)PyArray_DATA(py_atom_type);
    double (*rxyz)[3] = (double(*)[3])PyArray_DATA(py_positions);
    nat = PyArray_DIMS(py_positions)[0];
    ntyp = PyArray_DIMS(py_znu)[0];

    int types[nat];
    for (i = 0; i < nat; i++)
        types[i] = (int)ptypes[i];
    int znuc[ntyp];
    for (i = 0; i < ntyp; i++)
        znuc[i] = (int)pznuc[i];

    lseg = 4;
    nid = nat * lseg;
    fp = (double *) malloc(sizeof(double)*nid);

    get_fp_nonperiodic(nid, nat, ntyp, types, rxyz, znuc, fp);
    vec = PyList_New(0);
    for (i = 0; i < nid; i++) {
        pytmp = PyFloat_FromDouble( fp[i] );
        PyList_Append(vec, pytmp);
        Py_DECREF(pytmp);
    }
    free(fp);

    return vec;

}

static PyObject * py_get_periodic(PyObject *self, PyObject *args)
{
    int i, j, k, lmax, l, lseg, natx, flag, log, ldfp;
    double cutoff;
    double **sfp, **lfp, ****dfp;
    PyArrayObject* py_lattice;
    PyArrayObject* py_positions;
    PyArrayObject* py_atom_type;
    PyArrayObject* py_znu;
    PyObject* array, *vec, *vec1, *vec2, *shortfp, *longfp, ***outdfp;
    PyObject* pytmp;


    if (!PyArg_ParseTuple(args, "iiiOOOOiid", &flag, &ldfp, &log, &py_lattice, &py_positions, &py_atom_type, &py_znu,
                &lmax, &natx,  &cutoff))
        return NULL;


    double (*lat)[3] = (double(*)[3])PyArray_DATA(py_lattice);
    double (*rxyz)[3] = (double(*)[3])PyArray_DATA(py_positions);
    int nat = PyArray_DIMS(py_positions)[0];
    int ntyp = PyArray_DIMS(py_znu)[0];


    long* pznuc = (long*)PyArray_DATA(py_znu);
    long* ptypes = (long*)PyArray_DATA(py_atom_type);

    int types[nat];
    for ( i = 0; i < nat; i++){
        
        types[i] = (int)ptypes[i];
    }
    int znuc[ntyp];
    for ( i = 0; i < ntyp; i++)
        znuc[i] = (int)pznuc[i];

    if (lmax == 0) {
        lseg = 1; l = 1;
    } else if (lmax == 1) {
        lseg = 4; l = 2;
    } else {
        fprintf(stderr, "ORBITAL ERROR."); return NULL;
    }

   
    sfp = (double **) malloc(sizeof(double*)*nat);
    if (!sfp) {
        // Handle allocation failure
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory for sfp");
        return NULL;
    }
    lfp = (double **) malloc(sizeof(double*)*nat);
    if (!lfp) {
        // Handle allocation failure
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory for lfp");
        return NULL;
    }

    for ( i = 0; i < nat; i++ ) {
        sfp[i] = (double *) malloc(sizeof(double) * l * (ntyp + 1));
        lfp[i] = (double *) malloc(sizeof(double) * (natx * lseg));
    }

    // Zero-initialization
    for (i = 0; i < nat; i++) {
        memset(sfp[i], 0, sizeof(double) * l * (ntyp + 1));
        memset(lfp[i], 0, sizeof(double) * natx * lseg);
    }
    // dfp = (double ****) malloc(sizeof(double ***) * nat);
    // for ( i = 0; i < nat; i++ ) {
    //     dfp[i] = (double ***) malloc(sizeof(double **) * nat);
    //     for ( j = 0; j < nat; j++ ) {
    //         dfp[i][j] = (double **) malloc(sizeof(double *) * 3);
    //         for ( k = 0; k < 3; k++ ) {
    //             dfp[i][j][k] = (double *) malloc(sizeof(double) * (natx * lseg));
    //         }
    //     }
    // }
    // Assuming n_sphere is defined and represents the number of orbitals
    dfp = (double ****) malloc(sizeof(double ***) * nat);
    if (dfp == NULL) {
        fprintf(stderr, "Failed to allocate dfp\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < nat; i++) {
        dfp[i] = (double ***) malloc(sizeof(double **) * nat);
        if (dfp[i] == NULL) {
            fprintf(stderr, "Failed to allocate dfp[%d]\n", i);
            exit(EXIT_FAILURE);
        }
        for (j = 0; j < nat; j++) {
            dfp[i][j] = (double **) malloc(sizeof(double *) * 3);
            if (dfp[i][j] == NULL) {
                fprintf(stderr, "Failed to allocate dfp[%d][%d]\n", i, j);
                exit(EXIT_FAILURE);
            }
            for (k = 0; k < 3; k++) {
                dfp[i][j][k] = (double *) calloc(natx, sizeof(double));
                if (dfp[i][j][k] == NULL) {
                    fprintf(stderr, "Failed to allocate dfp[%d][%d][%d]\n", i, j, k);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    
    get_fp_periodic(flag, ldfp, log, lmax, nat, ntyp, types, lat, rxyz, znuc, natx,  cutoff, sfp, lfp, dfp);
  

    shortfp = PyList_New(0);
    longfp = PyList_New(0);
    for ( i = 0; i < nat; i++ ) {
        if (flag > 0) {
            vec = PyList_New(0);
            for ( j = 0; j < l*(ntyp+1); j++){
                pytmp = PyFloat_FromDouble( sfp[i][j] );
                PyList_Append(vec,  pytmp);
                Py_DECREF(pytmp);
            }
            PyList_Append( shortfp, vec );
            Py_DECREF(vec);
        }
        vec = PyList_New(0);
        for (j = 0; j < natx*lseg; j++) {
            pytmp = PyFloat_FromDouble( lfp[i][j] );
            PyList_Append(vec, pytmp );
            Py_DECREF(pytmp);
        }
        PyList_Append( longfp, vec);
        Py_DECREF(vec);
    }

    if (ldfp > 0) {
        outdfp = PyList_New(0);
        for ( i = 0; i < nat; i++ ) {
            vec = PyList_New(0);
            for ( j = 0; j < nat; j++ ) {
                vec1 = PyList_New(0);
                for ( k = 0; k < 3; k++ ) {
                    vec2 = PyList_New(0);
                    for ( l = 0; l < natx*lseg; l++ ) {
                        pytmp = PyFloat_FromDouble( dfp[i][j][k][l] );
                        PyList_Append(vec2, pytmp );
                        Py_DECREF(pytmp);
                    }
                    PyList_Append(vec1, vec2);
                    Py_DECREF(vec2);
                }
                PyList_Append(vec, vec1);
                Py_DECREF(vec1);

            }
            PyList_Append( outdfp, vec );
            Py_DECREF(vec);
        }
    }

    array = PyList_New(0);
    PyList_Append( array, shortfp );
    PyList_Append( array, longfp );
    if (ldfp > 0) {
        PyList_Append( array, outdfp );
        Py_DECREF(outdfp);
    }
    Py_DECREF(shortfp);
    Py_DECREF(longfp);

    for (i = 0; i < nat; i++) {
        free(sfp[i]);
        free(lfp[i]);
        sfp[i] = NULL;
        lfp[i] = NULL;
    }
    free(sfp);
    free(lfp);
    sfp = NULL;
    lfp = NULL;

    // free dfp
    for (i = 0; i < nat; i++) {
        for (j = 0; j < nat; j++) {
            for (k = 0; k < 3; k++) {
                free(dfp[i][j][k]);
                dfp[i][j][k] = NULL;
            }
            free(dfp[i][j]);
            dfp[i][j] = NULL;
        }
        free(dfp[i]);
        dfp[i] = NULL;
    }

    return array;

}


static PyObject * py_get_dist_periodic(PyObject *self, PyObject *args)
{
    int i, j, ntyp;
    double fpdist;
    PyArrayObject* py_fp1;
    PyArrayObject* py_fp2;
    PyArrayObject* py_atom_type;
    PyObject* array;
    PyObject* farray;
    int nat, fp_len;
    PyObject* pytmp;

    if (!PyArg_ParseTuple(args, "iOOO", &ntyp, &py_atom_type, &py_fp1, &py_fp2))
        return NULL;


    long* ptypes = (long*)PyArray_DATA(py_atom_type);   
    nat = PyArray_DIMS(py_fp1)[0];
    fp_len = PyArray_DIMS(py_fp1)[1];

    int types[nat], f[nat];
    for ( i = 0; i < nat; i++) {
        types[i] = (int)ptypes[i];
    }

    double (*fptt1)[fp_len] = (double(*)[fp_len])PyArray_DATA(py_fp1);
    double (*fptt2)[fp_len] = (double(*)[fp_len])PyArray_DATA(py_fp2);

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


    fpdist = get_fpdistance_periodic(nat, ntyp, types, fp_len, fp1, fp2, f);

    for ( i = 0; i < nat; i++) {
        free(fp1[i]);
        free(fp2[i]);
        fp1[i] = NULL;
        fp2[i] = NULL;
    }

    free(fp1);
    free(fp2);
    fp1 = NULL;
    fp2 = NULL;

    farray = PyList_New(0);
    for (i = 0; i < nat; i++){
        pytmp = PyLong_FromLong( (long)f[i] );
        PyList_Append(farray, pytmp);
        Py_DECREF(pytmp);
    }

    array = PyList_New(0);
    PyList_Append( array, PyFloat_FromDouble(fpdist));
    PyList_Append( array, farray);

    Py_DECREF(farray);

    // printf("%d %d \n ", nat, nat2);
    return array;

}







