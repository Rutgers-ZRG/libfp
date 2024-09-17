/*******************************************************************************
 * Copyright (C) 2024 Zhu research group
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
#include <float.h>


extern void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda,
                double* w, double* work, int* lwork, int* info );

void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K,
            double *ALPHA, double *A, int *LDA, double *B, int *LDB,
            double *BETA, double *C, int *LDC);

double ddot_(int *N, double *DX, int *INCX, double *DY, int *INCY);


void get_fp(int flag, int ldfp, int log, int nat, int ntyp, int ixyz, int nx, int lseg, int l, double lat[3][3],
        double rxyz[][3], int types[], double rcov[], double cutoff, double **lfp, double **sfp, double ****dfp)
{
    int iat, jat, icat, ix, iy, iz, il, i, j, k, i1, i2, i3, i4;
    int iats, iorb, iiat;
    int n_sphere, ityp_sphere, nid, nids = l*(ntyp+1);
    int n_sphere_min = 1000000;
    int n_sphere_max = 0;
    int ind[nx*lseg];
    int indori[nx];
    double xi, yi, zi, xj, yj, zj, d2;
    double cutoff2 = cutoff*cutoff;
    double rxyz_sphere[nx][3];
    double rcov_sphere[nx];
    double alpha[nx];
    double amp[nx];
    double damp[nx];
    double fc, wc, amp_tmp;

    double **om, **aom, *pvec, **pvecs;
    double ****dom, ***dvdr;
    double omx[nids][nids], omy[nids][nids];

    int lda, ldb, ldc, info, lwork;
    double wkopt;
    double *work, *w, *a;

    char trans = 'N';  // No transpose
    double aa = 1.0, beta = 0.0;
    int incx = 1, incy = 1;


    int mm, nn, kk; 
    double dalpha = 1.0;
    double dbeta = 0.0;



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
                            amp_tmp = pow((1.0 - d2*fc), (NC-1));
                            amp[n_sphere-1] = amp_tmp * (1.0 - d2*fc);
                            damp[n_sphere-1] = -2. * fc * NC * amp_tmp;
                            // printf("ampf %d %g %g\n", n_sphere-1, amp_tmp, amp[n_sphere-1]);
                            rxyz_sphere[n_sphere-1][0] = xj;
                            rxyz_sphere[n_sphere-1][1] = yj;
                            rxyz_sphere[n_sphere-1][2] = zj;
                            // printf("rxzy, %g %g %g\n", xj, yj, zj);
                            rcov_sphere[n_sphere-1] = rcov[jat];
                            indori[n_sphere-1] = jat;
                            alpha[n_sphere-1] = 0.5/(rcov[jat]*rcov[jat]);
                            // printf("COV %g %g\n", rcov[jat], alpha[n_sphere-1]);
                            if (jat == iat && ix == 0 && iy == 0 && iz == 0) {
                                ityp_sphere = 0;
                                icat = n_sphere-1;
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
        if ( (pvecs = (double **) malloc(sizeof(double)*nid)) == NULL) {
            fprintf(stderr, "Memory could not be allocated.");
            exit(1);
        }
        for (i = 0; i < nid; i++){
            if ( (pvecs[i] = (double *) malloc(sizeof(double)*nid)) == NULL) {
                fprintf(stderr, "Memory could not be allocated.");
                exit(1);
            }
        }
        if ( (aom = (double **) malloc(sizeof(double)*nid)) == NULL) {
            fprintf(stderr, "Memory could not be allocated.");
            exit(1);
        }
        for (i = 0; i < nid; i++){
            if ( (aom[i] = (double *) malloc(sizeof(double)*nid)) == NULL) {
                fprintf(stderr, "Memory could not be allocated.");
                exit(1);
            }
        }

        if ( (pvec = (double *) malloc(sizeof(double)*nid)) == NULL ){
            fprintf(stderr, "Memory could not be allocated.");
            exit(1);}

        
        creat_om( lseg, n_sphere, rxyz_sphere, rcov_sphere, om );


        for (i = 0; i < n_sphere; i++){
            for (j = 0; j < n_sphere; j++){
                if (lseg == 1) {
                    aom[i][j] = om[i][j] * amp[i] * amp[j];
                } else {
                    aom[4*i][4*j] = om[4*i][4*j] * amp[i] * amp[j];
                    aom[4*i][4*j+1] = om[4*i][4*j+1] * amp[i] * amp[j];
                    aom[4*i][4*j+2] = om[4*i][4*j+2] * amp[i] * amp[j];
                    aom[4*i][4*j+3] = om[4*i][4*j+3] * amp[i] * amp[j];

                     /* <pi|sj> */
                    aom[4*i+1][4*j] = om[4*i+1][4*j] * amp[i] * amp[j];
                    aom[4*i+2][4*j] = om[4*i+2][4*j] * amp[i] * amp[j];
                    aom[4*i+3][4*j] = om[4*i+3][4*j] * amp[i] * amp[j];

                    /* <pi|pj> */
                    aom[4*i+1][4*j+1] = om[4*i+1][4*j+1] * amp[i] * amp[j];
                    aom[4*i+1][4*j+2] = om[4*i+1][4*j+2] * amp[i] * amp[j];
                    aom[4*i+1][4*j+3] = om[4*i+1][4*j+3] * amp[i] * amp[j];
                    aom[4*i+2][4*j+1] = om[4*i+2][4*j+1] * amp[i] * amp[j];
                    aom[4*i+2][4*j+2] = om[4*i+2][4*j+2] * amp[i] * amp[j];
                    aom[4*i+2][4*j+3] = om[4*i+2][4*j+3] * amp[i] * amp[j];
                    aom[4*i+3][4*j+1] = om[4*i+3][4*j+1] * amp[i] * amp[j];
                    aom[4*i+3][4*j+2] = om[4*i+3][4*j+2] * amp[i] * amp[j];
                    aom[4*i+3][4*j+3] = om[4*i+3][4*j+3] * amp[i] * amp[j];
                }
            }
        }

        if ( (a = (double *) malloc(sizeof(double)*nid*nid)) == NULL) {
            fprintf(stderr, "Memory could not be allocated.");
            exit(1);}
        
        memset(a, 0, sizeof(double)*nid*nid);

        if ( (w = (double *) malloc(sizeof(double)*nid)) == NULL) {
            fprintf(stderr, "Memory could not be allocated.");
            exit(1);}

        for (i = 0; i < nid; i++) {
            for (j = 0; j< nid; j++) {
                a[i*nid + j] = aom[j][i];
               
            }
        }


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
            fprintf(stderr, "Error: Negative w");
            exit(1); }

        for (i = 0; i < nid; i++)
            pvec[i] = a[i + (nid-1)*nid];

        for (i = 0; i < nid; i++)
            for (j = 0; j < nid; j++)
                pvecs[i][j] = a[i + j*nid];

        for (i = 0; i < nid; i++)
            lfp[iat][i] = w[nid-1-i];
        for (i = nid; i < nx; i++)
            lfp[iat][i] = 0.0;


        free(w);
        
        free(work);


        if (ldfp > 0) {

            if ( (dom = (double ****) malloc(sizeof(double)*nid)) == NULL) {
                fprintf(stderr, "Memory could not be allocated.");
                exit(1);
            }
            for (i = 0; i < nid; i++){
                if ( (dom[i] = (double ***) malloc(sizeof(double)*3)) == NULL) {
                fprintf(stderr, "Memory could not be allocated.");
                exit(1);
                }
            }
            for (i = 0; i < nid; i++){
                for (j = 0; j < 3; j++){
                    if ( (dom[i][j] = (double **) malloc(sizeof(double)*nid)) == NULL) {
                        fprintf(stderr, "Memory could not be allocated.");
                        exit(1);
                    }
                }
            }
            for (i = 0; i < nid; i++){
                for (j = 0; j < 3; j++){
                    for (k = 0; k < nid; k++){
                        if ( (dom[i][j][k] = (double *) malloc(sizeof(double)*nid)) == NULL) {
                            fprintf(stderr, "Memory could not be allocated.");
                            exit(1);
                        }
                    }
                }
            }
            for (i = 0; i < nid; i++) {
                for (i1 = 0; i1 < 3; i1++) {
                    for (i2 = 0; i2 < nid; i2++) {
                        for (i3 = 0; i3 < nid; i3++) {
                            dom[i][i1][i2][i3] = 0.;
                        }
                    }
                }
            }
            get_dom(n_sphere, icat, rxyz_sphere, alpha, amp, damp, om, dom);


        
            double *matt, *tmpA;
            int ik;
            mm = nid;
            nn = nid;
            kk = nid;
            lda = nid;
            ldb = nid;
            ldc = nid;
            if ( (matt = (double *) malloc(sizeof(double)*nid*nid)) == NULL) {
                fprintf(stderr, "Memory could not be allocated.");
                exit(1);}
            if ( (tmpA = (double *) malloc(sizeof(double)*nid*nid)) == NULL) {
                fprintf(stderr, "Memory could not be allocated.");
                exit(1);}
            for (iats = 0; iats < n_sphere; iats++) {
                iiat = indori[iats];
                for (ik = 0; ik < 3; ik++) {
                    k = 0;
                    for (j = 0; j < nid; j++) {
                        for (i = 0; i < nid; i++) {
                            matt[k++] = dom[iats][ik][i][j];
                        }
                    }
                    // Call dgemm_
                    dgemm_(&trans, &trans, &mm, &nn, &kk,
                            &dalpha, matt, &lda,
                            a, &ldb,
                            &dbeta, tmpA, &ldc); 
                    int incx = 1;
                    int incy = 1;
                    for (iorb = 0; iorb < nid; iorb++){
                        int iiorb = nid - iorb - 1;
                        double *vec = a + iiorb * nid;
                        double *tmp_col = tmpA + iiorb * nid;

                        double dot = ddot_(&nid, vec, &incx, tmp_col, &incy);
                        
                        // update dfp
                        dfp[iat][iiat][ik][iorb] += dot;
                    }

                }
            }
            free(tmpA);
            free(matt);
               
            //free dom
            for (i = 0; i < nid; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    for (k=0; k<nid; k++)
                    {
                        free(dom[i][j][k]);
                    }
                    free(dom[i][j]);
                }
                free(dom[i]);
            }
            free(dom);

        }

        free(a);


        if (flag > 0 ){
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
            dsyev_("V", "U", &nids, a, &lda, w, &wkopt, &lwork, &info);
            lwork = (int)wkopt;
            work = (double*) malloc(lwork*sizeof(double));
            dsyev_("V", "U", &nids, a, &lda, w, work, &lwork, &info);


            for (i = 0; i < nids; i++)
                sfp[iat][i] = w[nids-1-i];

        
            free(a);
            free(w);
            free(work);
        }

        
        work = NULL;
        for (i = 0; i < nid; i++)
        {
            free(om[i]);
            free(aom[i]);
            free(pvecs[i]);

        }
        

        free(om);
        free(aom);
        free(pvec);
        free(pvecs);


    }
    if (log > 0)
        printf("min  %d, max  %d\n", n_sphere_min, n_sphere_max);
}


void creat_om(int lseg, int n_sphere, double rxyz_sphere[][3], double rcov_sphere[],
         double **om)
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
                    * exp(-d2 * r);
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
                    * exp(-d2 * r);

                /* <si|pj> */
                sji = pow(sqrt(4.0*r*rcov_sphere[iat]*rcov_sphere[jat]), 3) * exp(-d2 * r);
                stv = sqrt(8.0) * rcov_sphere[jat] * r * sji;
                om[4*iat][4*jat+1] = stv * xji;
                om[4*iat][4*jat+2] = stv * yji;
                om[4*iat][4*jat+3] = stv * zji;

                /* <pi|sj> */
                stv = sqrt(8.0) * rcov_sphere[iat] * r * sji * -1.0;
                om[4*iat+1][4*jat] = stv * xji;
                om[4*iat+2][4*jat] = stv * yji;
                om[4*iat+3][4*jat] = stv * zji;

                /* <pi|pj> */
                stv = -8.0 * rcov_sphere[iat]*rcov_sphere[jat] * r * r * sji;
                om[4*iat+1][4*jat+1] = stv * ( xji * xji - 0.5/r );
                om[4*iat+1][4*jat+2] = stv * ( yji * xji         );
                om[4*iat+1][4*jat+3] = stv * ( zji * xji         );
                om[4*iat+2][4*jat+1] = stv * ( xji * yji         );
                om[4*iat+2][4*jat+2] = stv * ( yji * yji - 0.5/r );
                om[4*iat+2][4*jat+3] = stv * ( zji * yji         );
                om[4*iat+3][4*jat+1] = stv * ( xji * zji         );
                om[4*iat+3][4*jat+2] = stv * ( yji * zji         );
                om[4*iat+3][4*jat+3] = stv * ( zji * zji - 0.5/r );
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

    dsyev_("V", "U", &n, a, &lda, w, &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    work = (double*) malloc(lwork*sizeof(double));
    dsyev_("V", "U", &n, a, &lda, w, work, &lwork, &info);
    if (info > 0 ) {
        fprintf(stderr, "Error: DSYEV ixyz");
        exit(1);
    }

    ixyz = (int)(sqrt(1.0/w[0])*cutoff + 1);

    free(work);
    work = NULL;

    return ixyz;
}


void get_dom(int n_sphere, int icat, double rxyz_sphere[][3], double alpha[],
        double amp[], double damp[], double **om, double ****dom)
{   
    int iat, jat;
    double xi, yi, zi, xj, yj, zj, xij, yij, zij, r2, a1, a2, aa;
    double xic, yic, zic, xjc, yjc, zjc, pij, dipj, djpi;
    double dix, djx, dcx, diy, djy, dcy, diz, djz, dcz;

    // derivatives
    // <s|s>
    
    for (jat = 0; jat < n_sphere; jat++) {
        xj = rxyz_sphere[jat][0];
        yj = rxyz_sphere[jat][1];
        zj = rxyz_sphere[jat][2];
        

        for (iat = 0; iat < n_sphere; iat++) {
            xi = rxyz_sphere[iat][0];
            yi = rxyz_sphere[iat][1];
            zi = rxyz_sphere[iat][2];
            xij = xi - xj;
            yij = yi - yj;
            zij = zi - zj;
            r2 = xij*xij + yij*yij + zij*zij;
            a1 = alpha[iat] * alpha[jat];
            a2 = alpha[iat] + alpha[jat];
            
            aa = -2. * a1 / a2;
            xic = xi - rxyz_sphere[icat][0];
            yic = yi - rxyz_sphere[icat][1];
            zic = zi - rxyz_sphere[icat][2];
            xjc = xj - rxyz_sphere[icat][0];
            yjc = yj - rxyz_sphere[icat][1];
            zjc = zj - rxyz_sphere[icat][2];


            pij = amp[iat] * amp[jat];
            dipj = damp[iat] * amp[jat];
            djpi = damp[jat] * amp[iat];

            
            dix = pij * aa * om[iat][jat] * xij + dipj * om[iat][jat] * xic;
            djx = -pij * aa * om[iat][jat] * xij + djpi * om[iat][jat] * xjc;
            dcx = -dipj * om[iat][jat] * xic - djpi * om[iat][jat] * xjc;

            diy = pij * aa * om[iat][jat] * yij + dipj * om[iat][jat] * yic;
            djy = -pij * aa * om[iat][jat] * yij + djpi * om[iat][jat] * yjc;
            dcy = -dipj * om[iat][jat] * yic - djpi * om[iat][jat] * yjc;
            
            diz = pij * aa * om[iat][jat] * zij + dipj * om[iat][jat] * zic;
            djz = -pij * aa * om[iat][jat] * zij + djpi * om[iat][jat] * zjc;
            dcz = -dipj * om[iat][jat] * zic - djpi * om[iat][jat] * zjc;  

            dom[iat][0][iat][jat] += dix;
            dom[jat][0][iat][jat] += djx;
            dom[icat][0][iat][jat] += dcx;

            dom[iat][1][iat][jat] += diy;
            dom[jat][1][iat][jat] += djy;
            dom[icat][1][iat][jat] += dcy;

            dom[iat][2][iat][jat] += diz;
            dom[jat][2][iat][jat] += djz;
            dom[icat][2][iat][jat] += dcz;         


        }
    }

}
