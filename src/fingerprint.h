/*******************************************************************************
 * Copyright (C) 2015 Li Zhu 
 * All rights reserved. 
 * 
 * fingerprint.h
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

#define NC  2 
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
/*
void dsyev( char* jobz, char* uplo, int* n, double* a, int* lda,
        double* w, double* work, int* lwork, int* info );

void dsygv(int* itype, char* jobz, char* uplo, int* n, double* a, 
        int* lda, double* b, int* ldb, double* w, double* work, int* lwork, int* info);
*/


void get_fp(int flag, int ldfp, int log, int nat, int ntyp, int ixyz, int nx, int lseg, int l, double lat[3][3],
        double rxyz[][3], int types[], double rcov[], double cutoff, double **lfp, double **sfp, double ****dfp);

void creat_om(int lseg, int n_sphere, double rxyz_sphere[][3], double rcov_sphere[], 
         double **om);

void get_dom(int n_sphere, int icat, double rxyz_sphere[][3], double alpha[],
        double amp[], double damp[], double **om, double ****dom);

int get_ixyz(double lat[3][3], double cutoff);


