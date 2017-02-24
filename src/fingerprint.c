/*******************************************************************************
 * Copyright (C) 2015 Li Zhu
 * All rights reserved.
 *
 * fingerprint.c
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
#include "fingerprint.h"

void get_fp(int flag, int nat, int ntyp, int ixyz, int nx, int lseg, int l, double lat[3][3],
        double rxyz[][3], int types[], double rcov[], double cutoff, double **lfp, double **sfp)
{
    int iat, jat, ix, iy, iz, il, i, j;
    int n_sphere, ityp_sphere, nid, nids = l*(ntyp+1);
    int n_sphere_min = 1000000;
    int n_sphere_max = 0;
    int ind[nx*lseg];
    double xi, yi, zi, xj, yj, zj, d2;
    double cutoff2 = cutoff*cutoff;
    double rxyz_sphere[nx][3];
    double rcov_sphere[nx];
    double amp[nx];
    double fc, wc;

    double **om, *pvec;
    double omx[nids][nids], omy[nids][nids];

    int lda, info, lwork;
    double wkopt;
    double *work, *w, *a;


    wc = cutoff / sqrt(2.0 * NC);
    fc = 1.0 / (2.0 * NC * wc * wc);

    for (iat = 0; iat < nat; iat++) {
        xi = rxyz[iat][0];
        yi = rxyz[iat][1];
        zi = rxyz[iat][2];
        n_sphere = 0;
        for (jat = 0; jat < nat; jat++) {
            for (ix = -ixyz; ix <= ixyz; ix++) {
                for (iy = -ixyz; iy <= ixyz; iy++) {
                    for (iz = -ixyz; iz <= ixyz; iz++) {
                        xj = rxyz[jat][0] + ix*lat[0][0] + iy*lat[1][0] + iz*lat[2][0];
                        yj = rxyz[jat][1] + ix*lat[0][1] + iy*lat[1][1] + iz*lat[2][1];
                        zj = rxyz[jat][2] + ix*lat[0][2] + iy*lat[1][2] + iz*lat[2][2];
                        d2 = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi) + (zj-zi)*(zj-zi);
                        if (d2 <= cutoff2) {
                            n_sphere++;
                            if (n_sphere > nx) {
                                fprintf(stderr, "Error: cutoff is too large.");
                                exit(1);
                            }
                            amp[n_sphere-1] = pow((1.0 - d2*fc), NC);
                            rxyz_sphere[n_sphere-1][0] = xj;
                            rxyz_sphere[n_sphere-1][1] = yj;
                            rxyz_sphere[n_sphere-1][2] = zj;
                            rcov_sphere[n_sphere-1] = rcov[jat];
                            if (jat == iat && ix == 0 && iy == 0 && iz == 0) {
                                ityp_sphere = 0;
                            } else {
                                ityp_sphere = types[jat];
                            }

                            for (il = 0; il < lseg; il++) {
                                if (il == 0){
                                    ind[il+lseg*(n_sphere-1)] = ityp_sphere * l + 0;
                                } else {
                                    ind[il+lseg*(n_sphere-1)] = ityp_sphere * l + 1;
                                }
                            }

                        }
                    }
                }
            }
        }
        n_sphere_min = MIN(n_sphere_min, n_sphere);
        n_sphere_max = MAX(n_sphere_max, n_sphere);


        // big overlap matrix

        nid = lseg * n_sphere;

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

        if ( (pvec = (double *) malloc(sizeof(double)*nid)) == NULL ){
            fprintf(stderr, "Memory could not be allocated.");
            exit(1);}


        creat_om( lseg, n_sphere, rxyz_sphere, rcov_sphere, amp, om );


        if ( (a = (double *) malloc(sizeof(double)*nid*nid)) == NULL) {
            fprintf(stderr, "Memory could not be allocated.");
            exit(1);}

        if ( (w = (double *) malloc(sizeof(double)*nid)) == NULL) {
            fprintf(stderr, "Memory could not be allocated.");
            exit(1);}

        for (i = 0; i < nid; i++)
            for (j = 0; j< nid; j++)
                a[i*nid + j] = om[j][i];



        lda = nid;
        lwork = -1;

        /// for (i = 0; i< nid*nid; i++)
        ///     printf("a[%d]= %g\n", i, a[i]);

        dsyev("V", "U", &nid, a, &lda, w, &wkopt, &lwork, &info);
        lwork = (int)wkopt;
        work = (double*) malloc(lwork*sizeof(double));
        dsyev("V", "U", &nid, a, &lda, w, work, &lwork, &info);
        if (info > 0 ) {
            fprintf(stderr, "Error: DSYEV 1");
            exit(1); }
        if (w[0] < -1E-12) {
            printf("w[0] = %g\n", w[0]);
            fprintf(stderr, "Error: Negative w");
            exit(1); }

        for (i = 0; i < nid; i++)
            pvec[i] = a[i + (nid-1)*nid];

        if (flag <= 0 ) {
            for (i = 0; i < nid; i++)
                lfp[iat][i] = w[nid-1-i];
            for (i = nid; i < nx; i++)
                lfp[iat][i] = 0.0;
        }

        free(w);
        free(a);
        free(work);


        if (flag >= 0 ){
            /* contract */
            nids = l*(ntyp+1);

            for (i = 0; i < nids; i++) {
                for (j = 0; j < nids; j++) {
                    omx[i][j] = 0.0;
                    omy[i][j] = 0.0;
                }
            }

            for (i = 0; i < nid; i++)
                for (j = 0; j < nid; j++)
                    omx[ind[i]][ind[j]] = omx[ind[i]][ind[j]] + pvec[i] * om[i][j] * pvec[j];


            if ( (a = (double *) malloc(sizeof(double)*nids*nids)) == NULL) {
                fprintf(stderr, "Memory could not be allocated.");
                exit(1);}

            if ( (w = (double *) malloc(sizeof(double)*nids)) == NULL) {
                fprintf(stderr, "Memory could not be allocated.");
                exit(1);}

            for (i = 0; i < nids; i++){
                for (j = 0; j< nids; j++){
                    a[i*nids + j] = omx[j][i];
                }
            }


            lda = nids;
            lwork = -1;
            dsyev("V", "U", &nids, a, &lda, w, &wkopt, &lwork, &info);
            lwork = (int)wkopt;
            work = (double*) malloc(lwork*sizeof(double));
            dsyev("V", "U", &nids, a, &lda, w, work, &lwork, &info);


            for (i = 0; i < nids; i++)
                sfp[iat][i] = w[nids-1-i];

            /*
               printf("iat = %d\n", iat);
               for (i = 0; i < nids; i++)
               printf(" %e ", sfp[iat][i]);
               printf("\n");
               */

            free(a);
            free(w);
        }

        free(work);
        for (i = 0; i < n_sphere; i++)
            free(om[i]);
        free(om);
        free(pvec);

    }
    printf("min  %d, max  %d\n", n_sphere_min, n_sphere_max);
}

void creat_om(int lseg, int n_sphere, double rxyz_sphere[][3], double rcov_sphere[],
        double amp[], double **om)
{

    int iat, jat, i, j;
    double xi, yi, zi, xj, yj, zj, xji, yji, zji;
    double d2, r, sji, stv;

    if (lseg == 1) { /* s orbital only */
        for (iat = 0; iat < n_sphere; iat++) {
            xi = rxyz_sphere[iat][0];
            yi = rxyz_sphere[iat][1];
            zi = rxyz_sphere[iat][2];

            for (jat = 0; jat < n_sphere; jat++) {
                xj = rxyz_sphere[jat][0];
                yj = rxyz_sphere[jat][1];
                zj = rxyz_sphere[jat][2];
                d2 = (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + (zi-zj)*(zi-zj);
                r = 0.5/(rcov_sphere[iat]*rcov_sphere[iat] + rcov_sphere[jat]*rcov_sphere[jat]);
                om[iat][jat] = pow(sqrt(4.0 * r * (rcov_sphere[iat]*rcov_sphere[jat])), 3)
                    * exp(-d2 * r) * amp[iat] * amp[jat];
            }
        }
    } else { /* both s and p orbital */
        for (iat = 0; iat < n_sphere; iat++) {
            xi = rxyz_sphere[iat][0];
            yi = rxyz_sphere[iat][1];
            zi = rxyz_sphere[iat][2];

            for (jat = 0; jat < n_sphere; jat++) {
                xj = rxyz_sphere[jat][0];
                yj = rxyz_sphere[jat][1];
                zj = rxyz_sphere[jat][2];
                xji = xj - xi;
                yji = yj - yi;
                zji = zj - zi;
                d2 = xji*xji + yji*yji + zji*zji;
                r = 0.5/(rcov_sphere[iat] * rcov_sphere[iat] + rcov_sphere[jat] * rcov_sphere[jat]);
                om[4*iat][4*jat] = pow(sqrt(4.0 * r * (rcov_sphere[iat]*rcov_sphere[jat])), 3)
                    * exp(-d2 * r) * amp[iat] * amp[jat];

                /* <si|pj> */
                sji = pow(sqrt(4.0*r*rcov_sphere[iat]*rcov_sphere[jat]), 3) * exp(-d2 * r);
                stv = sqrt(8.0) * rcov_sphere[jat] * r * sji;
                om[4*iat][4*jat+1] = stv * xji * amp[iat] * amp[jat];
                om[4*iat][4*jat+2] = stv * yji * amp[iat] * amp[jat];
                om[4*iat][4*jat+3] = stv * zji * amp[iat] * amp[jat];

                /* <pi|sj> */
                stv = sqrt(8.0) * rcov_sphere[iat] * r * sji * -1.0;
                om[4*iat+1][4*jat] = stv * xji * amp[iat] * amp[jat];
                om[4*iat+2][4*jat] = stv * yji * amp[iat] * amp[jat];
                om[4*iat+3][4*jat] = stv * zji * amp[iat] * amp[jat];

                /* <pi|pj> */
                stv = -8.0 * rcov_sphere[iat]*rcov_sphere[jat] * r * r * sji;
                om[4*iat+1][4*jat+1] = stv * ( xji * xji - 0.5/r ) * amp[iat] * amp[jat];
                om[4*iat+1][4*jat+2] = stv * ( yji * xji         ) * amp[iat] * amp[jat];
                om[4*iat+1][4*jat+3] = stv * ( zji * xji         ) * amp[iat] * amp[jat];
                om[4*iat+2][4*jat+1] = stv * ( xji * yji         ) * amp[iat] * amp[jat];
                om[4*iat+2][4*jat+2] = stv * ( yji * yji - 0.5/r ) * amp[iat] * amp[jat];
                om[4*iat+2][4*jat+3] = stv * ( zji * yji         ) * amp[iat] * amp[jat];
                om[4*iat+3][4*jat+1] = stv * ( xji * zji         ) * amp[iat] * amp[jat];
                om[4*iat+3][4*jat+2] = stv * ( yji * zji         ) * amp[iat] * amp[jat];
                om[4*iat+3][4*jat+3] = stv * ( zji * zji - 0.5/r ) * amp[iat] * amp[jat];
            }
        }
    }


    for (i = 0; i < lseg*n_sphere; i++){
        for (j = 0; j< lseg*n_sphere; j++){
            if (fabs(om[i][j] - om[j][i]) > 1e-6) {
                printf("OM SYMMETRY ERROR\n");
                printf("OM[%d][%d]=%g, OM[%d][%d]=%g\n", i, j, om[i][j], j, i, om[j][i]);
            }
        }
    }

}

int get_ixyz(double lat[3][3], double cutoff)
{
    int i, j, ixyz;
    int n = 3, lda = 3, info, lwork;
    double* work;
    double wkopt, w[3], a[9], lat2[3][3];

    lwork = -1;

    for (i=0; i < 3; i++)
        for (j=0; j < 3; j++)
            lat2[i][j] = lat[i][0]*lat[j][0] + lat[i][1]*lat[j][1] + lat[i][2]*lat[j][2];

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            a[i*3+j] = lat2[j][i];

    dsyev("V", "U", &n, a, &lda, w, &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    work = (double*) malloc(lwork*sizeof(double));
    dsyev("V", "U", &n, a, &lda, w, work, &lwork, &info);
    if (info > 0 ) {
        fprintf(stderr, "Error: DSYEV ixyz");
        exit(1);
    }

    ixyz = (int)(sqrt(1.0/w[0])*cutoff + 1);

    free(work);
    work = NULL;

    return ixyz;
}


