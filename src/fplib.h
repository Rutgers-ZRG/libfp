/* Copyright (C) 2015 Li Zhu */
/* All rights reserved. */

/* fplib.h */
/* This file is part of fplib. */

/* This file is distributed under the terms of the            */
/* GNN Lesser General Public License Version 3.0.             */
/* See ../LICENSE or http://www.gnu.org/licenses/lgpl-3.0.txt */


void get_fp_periodic(int lmax, int nat, int ntyp, int types[], double lat[3][3],
        double rxyz[][3], char *symb[], int natx, double cutoff, double **sfp, double **lfp);

void get_fp_periodic_short(int lmax, int nat, int ntyp, int types[], double lat[3][3],
        double rxyz[][3], char *symb[], int natx, double cutoff, double **sfp);

void get_fp_periodic_long(int lmax, int nat, int ntyp, int types[], double lat[3][3],
        double rxyz[][3], char *symb[], int natx, double cutoff, double **lfp);


