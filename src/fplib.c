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
#include "version.h"

extern void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda,
                        double* w, double* work, int* lwork, int* info );

/*-----------------------------------------*/
/* Version: fplib-[major].[minor].[micro] */
/*-----------------------------------------*/
int get_major_version(void)
{
  return FPLIB_MAJOR_VERSION;
}

int get_minor_version(void)
{
  return FPLIB_MINOR_VERSION;
}

int get_micro_version(void)
{
  return FPLIB_MICRO_VERSION;
}



void get_fp_nonperiodic(int nid, int nat, int ntyp, int types[], double rxyz[][3], int znucl[], double fp[])
{
    
    int i, j;
    int lwork, lda, info;
    //double om[nid][nid];
    double **om;
    double amp[nat];
    double a[nid*nid];
    double w[nid];
    double wkopt;
    double *work;
    double rcov[nat];


    for (i = 0; i < nat; i++)
        rcov[i] = get_rcov( znucl[ types[i] - 1 ] );

    for (i = 0; i < nat; i++)
        amp[i] = 1.0;

    nid = 4 * nat;

            if ( (om = (double **) malloc(sizeof(double)*nid)) == NULL) {
            fprintf(stderr, "Memory could not be allocated.");
            exit(1);
        }
        for (i = 0; i < nid; i++){
            if ( (om[i] = (double *) malloc(sizeof(double)*nid)) == NULL) {
                fprintf(stderr, "Memory could not be allocated.");
                exit(1);
            }
        }
    creat_om(4, nat, rxyz, rcov, om);

    for (i = 0; i < nid; i++)
        for (j = 0; j< nid; j++)
            a[i*nid + j] = om[j][i];

    lda = nid;
    lwork = -1;
    dsyev_("V", "U", &nid, a, &lda, w, &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    work = (double*) malloc(lwork*sizeof(double));
    dsyev_("V", "U", &nid, a, &lda, w, work, &lwork, &info);
    if (info > 0 ) {
        fprintf(stderr, "Error: DSYEV 1");
        exit(1); }
    if (w[0] < -1E-12) {
        printf("w[0] = %g\n", w[0]);
        fprintf(stderr, "Error: Negative w");
        exit(1); }

    for (i = 0; i < nid; i++)
        fp[i] = w[nid-1-i];

    free(work);
            for (i = 0; i < nid; i++)
            free(om[i]);
        free(om);


}



void get_fp_periodic(int flag, int ldfp, int log, int lmax, int nat, int ntyp, int types[], double lat[3][3],
        double rxyz[][3], int znucl[], int natx, double cutoff, double **sfp, double **lfp, double ****dfp)
{
    int i, ixyz;
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
    {
        rcov[i] = get_rcov( znucl[ types[i] - 1 ] );
    }

    ixyz = get_ixyz(lat, cutoff);

    /* flag = 0: long fp only;  = 1: long and short fp */
    get_fp(flag, ldfp, log, nat, ntyp, ixyz, natx, lseg, l, lat, rxyz, types, rcov,  cutoff, lfp, sfp, dfp);

}
/*
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
                        costmp[i-1][j-1] = sqrt(tt/fp_len);
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
    fpd /= nat;

    return fpd;
}
*/
double get_fpdistance_periodic(int nat, int ntyp, int types[], int fp_len,
        double **fp1, double **fp2, int f[])
{
    if (nat <= 0 || fp_len <= 0) {
        // Handle invalid input
        return -1.0;
    }

    double fpd = 0.0, cc = 0.0, tt = 0.0, **costmp = NULL, *a = NULL;
    int iat, jat, ityp, i, j, k, ii, jj, inat = 0, *ft = NULL;

    // Allocate costmp
    costmp = (double **)malloc(nat * sizeof(double *));
    if (!costmp) {
        // Handle allocation failure
        return -1.0;
    }
    for (i = 0; i < nat; i++) {
        costmp[i] = (double *)malloc(nat * sizeof(double));
        if (!costmp[i]) {
            // Handle allocation failure
            for (j = 0; j < i; j++) free(costmp[j]);
            free(costmp);
            return -1.0;
        }
    }
   
    for (ityp = 1; ityp <= ntyp; ityp++) {
        i = 0;

        // Count number of atoms of type 'ityp'
        for (iat = 0; iat < nat; iat++) {
            if (types[iat] == ityp) {
                i++;
            }
        }

        if (i == 0) {
            continue; // No atoms of this type
        }

        // Allocate ft and a
        ft = (int *)malloc(sizeof(int) * i);
        if (!ft) {
            // Handle allocation failure
            // Free costmp
            for (j = 0; j < nat; j++) free(costmp[j]);
            free(costmp);
            return -1.0;
        }
        a = (double *)malloc(sizeof(double) * i * i);
        if (!a) {
            // Handle allocation failure
            free(ft);
            for (j = 0; j < nat; j++) free(costmp[j]);
            free(costmp);
            return -1.0;
        }
        
        int idx_i = 0;
        for (iat = 0; iat < nat; iat++) {
            if (types[iat] == ityp) {
                int idx_j = 0;
                for (jat = 0; jat < nat; jat++) {
                    if (types[jat] == ityp) {
                        tt = 0.0;
                        for (k = 0; k < fp_len; k++) {
                            double diff = fp1[iat][k] - fp2[jat][k];
                            tt += diff * diff;
                        }
                        costmp[idx_i][idx_j] = sqrt(tt / fp_len);
                        idx_j++;
                    }
                }
                idx_i++;
            }
        }
        
        // Copy costmp into a
        for (ii = 0; ii < i; ii++) {
            for (jj = 0; jj < i; jj++) {
                a[ii * i + jj] = costmp[ii][jj];
            }
        }

        // Call apc function
        
        int apc_status = apc(i, a, &cc, ft);
        
        if (apc_status != 0) {
            // Handle error
            free(ft);
            free(a);
            for (j = 0; j < nat; j++) free(costmp[j]);
            free(costmp);
            return -1.0;
        }

        // Update f and inat
        for (k = 0; k < i; k++) {
            if (inat >= nat) {
                // Prevent buffer overflow
                break;
            }
            f[inat] = ft[k];
            inat++;
        }

        free(ft);
        free(a);
        ft = NULL;
        a = NULL;
        fpd += cc;
    }

    if (nat == 0) {
        // Prevent division by zero
        for (i = 0; i < nat; i++) free(costmp[i]);
        free(costmp);
        return -1.0;
    }

    fpd /= nat;

    // Free costmp
    for (i = 0; i < nat; i++) free(costmp[i]);
    free(costmp);

    return fpd;
}