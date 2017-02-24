/*******************************************************************************
 * Copyright (C) 2015 Li Zhu 
 * All rights reserved. 
 * 
 * fplib.c
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fplib.h"
#include "rcov.h"
#include "fingerprint.h"
#include "apc.h"

void get_fp_periodic(int lmax, int nat, int ntyp, int types[], double lat[3][3],
        double rxyz[][3], int znucl[], int natx, double cutoff, double **sfp, double **lfp)
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
        rcov[i] = get_rcov( znucl[ types[i] - 1 ] );

    ixyz = get_ixyz(lat, cutoff);

    /* flag < 0: long fp only;  == 0: short fp only; > 0: long and short fp */
    get_fp(flag, nat, ntyp, ixyz, natx, lseg, l, lat, rxyz, types, rcov,  cutoff, lfp, sfp);

}

void get_fp_periodic_short(int lmax, int nat, int ntyp, int types[], double lat[3][3],
        double rxyz[][3], int znucl[], int natx, double cutoff, double **sfp)
{
    int i, ixyz, flag = 1;
    double rcov[nat], **lfp = NULL;
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
        rcov[i] = get_rcov( znucl[ types[i] - 1 ] );

    ixyz = get_ixyz(lat, cutoff);

    /* flag < 0: long fp only; > 0 short fp only; == 0: long and short fp */
    get_fp(flag, nat, ntyp, ixyz, natx, lseg, l, lat, rxyz, types, rcov,  cutoff, lfp, sfp);
}

void get_fp_periodic_long(int lmax, int nat, int ntyp, int types[], double lat[3][3],
        double rxyz[][3], int znucl[], int natx, double cutoff, double **lfp)
{
    int i, ixyz, flag = -1;
    double rcov[nat], **sfp = NULL;
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
        rcov[i] = get_rcov( znucl[ types[i] - 1 ] );

    ixyz = get_ixyz(lat, cutoff);

    /* flag < 0: long fp only;  == 0: short fp only; > 0: long and short fp */
    get_fp(flag, nat, ntyp, ixyz, natx, lseg, l, lat, rxyz, types, rcov,  cutoff, lfp, sfp);
}

double get_fpdistance_periodic(int nat, int ntyp, int types[], int fp_len, 
        double **fp1, double **fp2, int f[])
{
    double fpd, cc, tt, costmp[nat][nat], *a;
    int iat, jat, ityp, i, j, k, ii, jj, inat, *ft;

    fpd = 0.0;
    inat = 0;

    for (ityp = 1; ityp <= ntyp; ityp++) {
        i = 0;
        for (iat = 0; iat < nat; iat++) {
            if (types[iat] == ityp) {
                i++;
                j = 0;
                for (jat = 0; jat < nat; jat++){
                    if (types[jat] == ityp) {
                        j++;
                        tt = 0.0;
                        for (k = 0; k < fp_len; k++)
                            tt += (fp1[iat][k] - fp2[jat][k]) * (fp1[iat][k] - fp2[jat][k]); 
                        costmp[i-1][j-1] = sqrt(tt);
                    }
                }
            }
        }
        ft = (int *) malloc(sizeof(int)*i);
        a = (double *) malloc(sizeof(double)*i*i);
        //printf("nt %d %d\n", i,j);
       
        for (ii = 0; ii < i; ii++)
            for (jj = 0; jj < i; jj++)
                a[ii*i + jj] = costmp[ii][jj];
        apc(i, a, &cc, ft);

        for (k = 0; k < i; k++){
            f[inat] = ft[k];
            inat++;
        }

        free(ft);
        free(a);
        ft = NULL;
        a = NULL;
        fpd += cc;

    }
    fpd /= sqrt(nat);

    return fpd;
}

