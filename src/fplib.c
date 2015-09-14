// Copyright (C) 2015 Li Zhu 
// All rights reserved. 

// This file is part of fplib. 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fplib.h"
#include "cell.h"

extern void dsyev( char* jobz, char* uplo, int* n, double* a, int* lda,
                double* w, double* work, int* lwork, int* info );

int get_fingerprint_periodic(int nat, int ntyp, int types[], double lat[3][3],
        double rxyz[][3], char *symb[], int natx, char *orb, double cutoff,
        double sfp[][], double lfp[][], double e=0.0)
{
    int iat, jat, ix, iy, iz, i, j, k;
    int n, lda, info, lwork;
    double wkopt;
    double* work;
    double* w;
    double* a;
    double factor_cutoff;
    double lat2[3][3];
    double eiglat[3];
    double worklat[9];

    Cell cell;
    int lseg, l;

    
    if (strcmp(s, "s") == 0){
        lseg = 1;
        l = 1;
    }else if (strcmp(s, "sp") == 0){
        lseg = 4;
        l = 2;
    }else {
        fprintf(stderr, "Error: ORBITAL.");
        exit(1);
    }


    cell = c_new_cell(nat, ntyp, natx, lseg, l);
    c_set_cell(&cell, lat, rxyz, types, symb, e);

    n = 3;
    lda = n;
    lwork = -1;
    w = (double *) malloc(n*sizeof(double));
    a = (double *) malloc(n*n*sizeof(double));
  
    for (i=0; i < 3; i++) 
        for (j=0; j < 3; j++) 
            lat2[i][j] = lat[i][0]*lat[j][0] + lat[i][1]*lat[j][1] + lat[i][3]*lat[j][3];

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
   
    ixyz = (int)(sqrt(1.0/w[0])*cutoff + 1)

    free(w);
    free(a);
    free(work);

    // get sphere info
    
    info = get_fp(&cell, ixyz, cutoff);

}
 
int get_fp(Cell * cell, int ixyz, int nx, int lseg, int l, double cutoff){
    int iat, jat, ix, iy, iz, il, i, j;
    int n_sphere, ityp_sphere, nid; 
    int n_sphere_min = 1000000; 
    int n_sphere_max = 0;
    int ind[nx];
    double xi, yi, zi, xj, yj, zj, il;
    double cutoff2 = cutoff*cutoff;
    double rxyz_sphere[nx][3];
    double rcov_sphere[nx];
    double fc, wc;

    double **om, **om_save;

    int n, lda, info, lwork;
    double wkopt;
    double* work;
    double* w;
    double* a;


    wc = cutoff/sqrt(2.0*NC);
    fc = 1.0/(2.0*NC*wc*wc);

    for (iat = 0; iat < cell->n; iat++) {
        xi = cell->rxyz[iat][0];
        yi = cell->rxyz[iat][1];
        zi = cell->rxyz[iat][2];
        n_sphere = 0;
        for (jat = 0; jat < cell->n; jat++) {
            for (ix = -ixyz; ix <= ixyz; ix++) {
                for (iy = -ixyz; iy <= ixyz; iy++) {
                    for (iz = -ixyz; iz <= ixyz; iz++) {
                        xj = cell->rxyz[jat][0] + ix*cell->lat[0][0] 
                            + iy*cell->lat[1][0] + iz*cell->lat[2][0];
                        yj = cell->rxyz[jat][1] + ix*cell->lat[0][1]
                            + iy*cell->lat[1][1] + iz*cell->lat[2][1];
                        zj = cell->rxyz[jat][2] + ix*cell->lat[0][2]
                            + iy*cell->lat[1][2] + iz*cell->lat[2][2];
                        d2 = pow((xj-xi),2) + pow((yj-yi),2) + pow((zj-zi),2);
                        if (d2 <= cutoff2) {
                            n_sphere++;
                            if (n_sphere > nx) {
                                 fprintf(stderr, "Error: cutoff is too small.");
                                 return 1;
                            }
                            amp[n_sphere] = pow((1.0 - d2*fc), NC);
                            rxyz_sphere[n_sphere][0] = xj;
                            rxyz_sphere[n_sphere][1] = yj;
                            rxyz_sphere[n_sphere][2] = zj;
                            rcov_sphere[n_sphere] = cell->rcov[j];
                            if (jat == iat && ix == 0 && iy == 0 && iz == 0) {
                                ityp_sphere = 0;
                            } else {
                                ityp_sphere = cell->types[jat];
                            }
                            for (il = 0; il < lseg; il++) {
                                if (il == 0){
                                    ind[il+lseg*(nat_sphere-1)] = ityp_sphere*il + 1;
                                } else {
                                    ind[il+lseg*(nat_sphere-1)] = ityp_sphere*il + 2;
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

        nid = lseg*nat_sphere;

        if ( (om = (double **) malloc(sizeof(double)*nid)) == NULL) {
            printf(stderr, "Memory could not be allocated.");
            exit(1);
        }
        for (i = 0; i < nid; i++){
            if ( (om[i] = (double *) malloc(sizeof(double)*nid)) == NULL) {
                printf(stderr, "Memory could not be allocated.");
                exit(1);
            }
        }
        
        if ( (om_save = (double **) malloc(sizeof(double)*nid)) == NULL) {
            printf(stderr, "Memory could not be allocated.");
            exit(1); }
        for ( i = 0; i < nid; i++ ) {
            if ( (om_save[i] = (double *) malloc(sizeof(double)*nid)) == NULL ) {
                printf(stderr, "Memory could not be allocated.");
                exit(1); } }


        creat_om( lseg, n_sphere, rxyz_sphere, rcov_sphere, amp, &om );


        for (i = 0; i < nid; i++){
            for (j = 0; j < nid; j++){
                om_save[i][j] = om[i][j];
            }
        }

        if ( (a = (double *) malloc(sizeof(double)*nid*nid)) == NULL) {
            printf(stderr, "Memory could not be allocated.");
            exit(1);}
        
        if ( (w = (double *) malloc(sizeof(double)*nid)) == NULL) {
           printf(stderr, "Memory could not be allocated.");
           exit(1);}

        for (i = 0; i < nid; i++) 
            for (j = 0; j< nid; j++) 
                a[i*nid + j] = om[j][i];



        lda = nid;
        lwork = -1;

        dsyev("V", "U", &nid, a, &lda, w, &wkopt, &lwork, &info);
        lwork = (int)wkopt;
        work = (double*) malloc(lwork*sizeof(double));
        dsyev("V", "U", &nid, a, &lda, w, work, &lwork, &info);
        if (info > 0 ) {
            fprintf(stderr, "Error: DSYEV 1");
            exit(1); }
        if (w[0] < -1E-12) {
            fprintf(stderr, "Error: Negative w");
            exit(1); }
        
        for (i = 0; i < nid; i++) 
            cell->lfp[iat][i] = w[nid-1-i];
        for (i = nid; i < nx; i++)
            cell->lfp[iat][i] = 0.0;


        // contract

        
    }
}
