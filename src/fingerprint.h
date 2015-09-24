/* Copyright (C) 2015 Li Zhu */
/* All rights reserved. */

/* fingerprint.h */
/* This file is part of fplib. */

/* This file is distributed under the terms of the            */
/* GNN Lesser General Public License Version 3.0.             */
/* See ../LICENSE or http://www.gnu.org/licenses/lgpl-3.0.txt */

#define NC  3 
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

extern void dsyev( char* jobz, char* uplo, int* n, double* a, int* lda,
        double* w, double* work, int* lwork, int* info );

extern void dsygv(int* itype, char* jobz, char* uplo, int* n, double* a, 
        int* lda, double* b, int* ldb, double* w, double* work, int* lwork, int* info);

void get_fp(int flag, int nat, int ntyp, int ixyz, int nx, int lseg, int l, double lat[3][3],
        double rxyz[][3], int types[], double rcov[], double cutoff, double **lfp, double **sfp);


void creat_om(int lseg, int n_sphere, double rxyz_sphere[][3], double rcov_sphere[], 
        double amp[], double **om);

int get_ixyz(double lat[3][3], double cutoff);


