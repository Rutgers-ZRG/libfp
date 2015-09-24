/* Copyright (C) 2015 Li Zhu */
/* All rights reserved. */

/* fplib.c */
/* This file is part of fplib. */

/* This file is distributed under the terms of the            */
/* GNN Lesser General Public License Version 3.0.             */
/* See ../LICENSE or http://www.gnu.org/licenses/lgpl-3.0.txt */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fplib.h"
#include "rcov.h"
#include "fingerprint.h"

void get_fp_periodic(int lmax, int nat, int ntyp, int types[], double lat[3][3],
        double rxyz[][3], char *symb[], int natx, double cutoff, double **sfp, double **lfp)
{
    int i, ixyz, flag = 0;
    double rcov[nat];
    int lseg, l;

    if (lmax == 0){
        lseg = 1;
        l = 1;
    }else if (lmax == 1){
        lseg = 4;
        l = 2;
    }else {
        fprintf(stderr, "Error: ORBITAL.");
        exit(1);
    }

    for (i = 0; i < nat; i++)
        rcov[i] = get_rcov( symb[types[i]-1] );

    ixyz = get_ixyz(lat, cutoff);

    /* flag < 0: long fp only;  == 0: short fp only; > 0: long and short fp */
    get_fp(flag, nat, ntyp, ixyz, natx, lseg, l, lat, rxyz, types, rcov,  cutoff, lfp, sfp);

}

void get_fp_periodic_short(int lmax, int nat, int ntyp, int types[], double lat[3][3],
        double rxyz[][3], char *symb[], int natx, double cutoff, double **sfp)
{
    int i, ixyz, flag = 1;
    double rcov[nat], **lfp;
    int lseg, l;

    if (lmax == 0){
        lseg = 1;
        l = 1;
    }else if (lmax == 1){
        lseg = 4;
        l = 2;
    }else {
        fprintf(stderr, "Error: ORBITAL.");
        exit(1);
    }

    for (i = 0; i < nat; i++)
        rcov[i] = get_rcov( symb[types[i]-1] );

    ixyz = get_ixyz(lat, cutoff);

    /* flag < 0: long fp only; > 0 short fp only; == 0: long and short fp */
    get_fp(flag, nat, ntyp, ixyz, natx, lseg, l, lat, rxyz, types, rcov,  cutoff, lfp, sfp);
}

void get_fp_periodic_long(int lmax, int nat, int ntyp, int types[], double lat[3][3],
        double rxyz[][3], char *symb[], int natx, double cutoff, double **lfp)
{
    int i, ixyz, flag = -1;
    double rcov[nat], **sfp;
    int lseg, l;

    if (lmax == 0){
        lseg = 1;
        l = 1;
    }else if (lmax == 1){
        lseg = 4;
        l = 2;
    }else {
        fprintf(stderr, "Error: ORBITAL.");
        exit(1);
    }

    for (i = 0; i < nat; i++)
        rcov[i] = get_rcov( symb[types[i]-1] );

    ixyz = get_ixyz(lat, cutoff);

    /* flag < 0: long fp only;  == 0: short fp only; > 0: long and short fp */
    get_fp(flag, nat, ntyp, ixyz, natx, lseg, l, lat, rxyz, types, rcov,  cutoff, lfp, sfp);
}

