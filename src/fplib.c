// Copyright (C) 2015 Li Zhu
// All rights reserved.

// This file is part of fplib.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fplib.h"

void get_fingerprint_periodic(int nat, int ntyp, int types[], double lat[3][3],
        double rxyz[][3], char *symb[], int natx, char *orb, double cutoff,
        double **sfp, double **lfp)
{
    int iat, jat, ix, iy, iz, i, j, k, ixyz;
    int n, lda, info, lwork;
    double wkopt;
    double* work;
    double* w;
    double* a;
    double factor_cutoff;
    double lat2[3][3];
    double rcov[nat];

    //Cell * cell;
    int lseg, l;


    if (strcmp(orb, "s") == 0){
        lseg = 1;
        l = 1;
    }else if (strcmp(orb, "sp") == 0){
        lseg = 4;
        l = 2;
    }else {
        fprintf(stderr, "Error: ORBITAL.");
        exit(1);
    }


    for (i = 0; i < nat; i++){
       rcov[i] = get_rcov( symb[types[i]-1] );
       //printf("rcov[%d] = %g\n", i, rcov[i]);
    }


    n = 3;
    lda = n;
    lwork = -1;
    w = (double *) malloc(n*sizeof(double));
    a = (double *) malloc(n*n*sizeof(double));

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
         fprintf(stderr, "Error: DSYEV 1");
         exit(1);
    }

    ixyz = (int)(sqrt(1.0/w[0])*cutoff + 1);

    free(w);
    free(a);
    free(work);

    // get sphere info

    //info = get_fp(cell, ixyz, natx, lseg, l, cutoff);
    info = get_fp(nat, ixyz, natx, lseg, l, lat, rxyz, types, rcov,  cutoff, lfp, sfp);

    //for (i = 0; i < nat; i++)
    //    for (j = 0; j < natx; j++)
    //        printf("%d  %d\n",i, j);
    //        lfp[i][j] = cell->lfp[i][j];
    //lfp[0][99] = cell->lfp[0][99];
    //lfp[5][9] = cell->lfp[5][9];


    //c_del_cell(cell);
    //printf("G9\n");

}

int get_fp(int nat, int ixyz, int nx, int lseg, int l, double lat[3][3],
        double rxyz[][3], int types[], double rcov[], double cutoff, double **lfp, double **sfp)
{
    int iat, jat, ix, iy, iz, il, i, j;
    int n_sphere, ityp_sphere, nid;
    int n_sphere_min = 1000000;
    int n_sphere_max = 0;
    int ind[nx];
    double xi, yi, zi, xj, yj, zj, d2;
    double cutoff2 = cutoff*cutoff;
    double rxyz_sphere[nx][3];
    double rcov_sphere[nx];
    double amp[nx];
    double fc, wc;

    double **om;

    int n, lda, info, lwork;
    double wkopt;
    double* work;
    double* w;
    double* a;


    wc = cutoff/sqrt(2.0*NC);
    fc = 1.0/(2.0*NC*wc*wc);

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
                                 return 1;
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
                            //for (il = 0; il < lseg; il++) {
                            //    if (il == 0){
                            //        ind[il+lseg*(n_sphere-1)] = ityp_sphere*il + 1;
                            //    } else {
                            //        ind[il+lseg*(n_sphere-1)] = ityp_sphere*il + 2;
                            //    }
                            //}

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

        //if ( (om_save = (double **) malloc(sizeof(double)*nid)) == NULL) {
        //    fprintf(stderr, "Memory could not be allocated.");
        //    exit(1); }
        //for ( i = 0; i < nid; i++ ) {
        //    if ( (om_save[i] = (double *) malloc(sizeof(double)*nid)) == NULL ) {
        //        fprintf(stderr, "Memory could not be allocated.");
        //        exit(1); } }


        creat_om( lseg, n_sphere, rxyz_sphere, rcov_sphere, amp, om );


        //for (i = 0; i < nid; i++){
        //    for (j = 0; j < nid; j++){
        //        om_save[i][j] = om[i][j];
        //    }
        //}

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
            lfp[iat][i] = w[nid-1-i];
        for (i = nid; i < nx; i++)
            lfp[iat][i] = 0.0;

        printf("%d, ", iat);
        for (i = 0; i < nx; i++)
            printf("  %g  ", lfp[iat][i]);
        printf("\n");

        // contract


        free(w);
        free(work);
        free(a);
        for (i = 0; i < n_sphere; i++)
            free(om[i]);
        free(om);

    }
    printf("min  %d, max  %d\n", n_sphere_min, n_sphere_max);
    return 0;
}

void creat_om(int lseg, int n_sphere, double rxyz_sphere[][3], double rcov_sphere[], 
        double amp[], double **om)
{

    int iat, jat, i, j, ii, jj;
    double xi, yi, zi, xj, yj, zj, xji, yji, zji;
    double d2, r, sji, stv;
    //double om4[n_sphere][lseg][n_sphere][lseg];


    if (lseg == 1) {
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
    } else {
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
                d2 = xji*xji + yji*yji+zji*zji;
                r = 0.5/(rcov_sphere[iat]*rcov_sphere[iat] + rcov_sphere[jat]*rcov_sphere[jat]);
                // <sj|si>
                //om4[iat][0][jat][0] = pow(sqrt(4.0 * r * (rcov_sphere[iat]*rcov_sphere[jat])), 3)
                //printf("%g %g %g %g %g %g %g\n", r, 4.0 * r * (rcov_sphere[iat]*rcov_sphere[jat]), 
                //        sqrt(4.0 * r * (rcov_sphere[iat]*rcov_sphere[jat])), -d2*r,
                //     exp(-d2 * r), amp[iat], amp[jat]);
                om[4*iat][4*jat] = pow(sqrt(4.0 * r * (rcov_sphere[iat]*rcov_sphere[jat])), 3)
                        * exp(-d2 * r) * amp[iat] * amp[jat];
                printf("w %g \n", om[4*iat][4*jat]);
                // <pj|si>
                sji = pow(sqrt(4.0*r*rcov_sphere[iat]*rcov_sphere[jat]), 3) * exp(-d2 * r);
                stv = sqrt(8.0) * rcov_sphere[jat] * r * sji;
                //om4[iat][0][jat][1] = stv * xji * amp[iat] * amp[jat];
                om[4*iat][4*jat+1] = stv * xji * amp[iat] * amp[jat];
                //om4[iat][0][jat][2] = stv * yji * amp[iat] * amp[jat];
                om[4*iat][4*jat+2] = stv * yji * amp[iat] * amp[jat];
                //om4[iat][0][jat][3] = stv * zji * amp[iat] * amp[jat];
                om[4*iat][4*jat+3] = stv * zji * amp[iat] * amp[jat];
            
                stv = sqrt(8.0) * rcov_sphere[iat] * r * sji;
                //om4[iat][1][jat][0] = stv * xji * amp[iat] * amp[jat];
                om[4*iat+1][4*jat] = stv * xji * amp[iat] * amp[jat];
                //om4[iat][2][jat][0] = stv * yji * amp[iat] * amp[jat];
                om[4*iat+1][4*jat] = stv * yji * amp[iat] * amp[jat];
                //om4[iat][3][jat][0] = stv * zji * amp[iat] * amp[jat];
                om[4*iat+1][4*jat] = stv * zji * amp[iat] * amp[jat];

                // <pj|pi>
                stv = -8.0 * rcov_sphere[iat]*rcov_sphere[jat] * r * r * sji;
                //om4[iat][1][jat][1] = stv * ( xji * xji - 0.5/r ) * amp[iat] * amp[jat];
                //om4[iat][1][jat][2] = stv * ( yji * xji         ) * amp[iat] * amp[jat];
                //om4[iat][1][jat][3] = stv * ( zji * xji         ) * amp[iat] * amp[jat];
                //om4[iat][2][jat][1] = stv * ( xji * yji         ) * amp[iat] * amp[jat];
                //om4[iat][2][jat][2] = stv * ( yji * yji - 0.5/r ) * amp[iat] * amp[jat];
                //om4[iat][2][jat][3] = stv * ( zji * yji         ) * amp[iat] * amp[jat];
                //om4[iat][3][jat][1] = stv * ( xji * zji         ) * amp[iat] * amp[jat];
                //om4[iat][3][jat][2] = stv * ( yji * zji         ) * amp[iat] * amp[jat];
                //om4[iat][3][jat][3] = stv * ( zji * zji - 0.5/r ) * amp[iat] * amp[jat];

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

        //for (i = 0; i < n_sphere; i++) 
        //    for ( ii = 0; ii < 4; ii++) 
        //        for ( j = 0; j < n_sphere; j++) 
        //            for ( jj = 0; jj < 4; jj++) 
        //                om[4*i + ii][4*j + jj] = om4[i][ii][j][jj];


    }

    
    for (i = 0; i < 1*n_sphere; i++)
        for (j = 0; j< 1*n_sphere; j++)
            if (fabs(om[i][j] - om[j][i]) > 1e-6) 
                printf(" om  %d %d %g %g\n", i, j, om[i][j], om[j][i]);




}

