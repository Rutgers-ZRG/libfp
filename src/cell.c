// cell.c
// Copyright (C) 2015 Li Zhu
//

#include <stdlib.h>
#include <stdio.h>
#include "cell.h"
#include "mathfunc.h"

Cell c_new_cell(int n, int ntyp, int nx, int lseg, int l) {
    int i;
    Cell cell;
    cell.n = n;
    cell.ntyp = ntyp;
    if ((cell.types = (int *) malloc(sizeof(int) * size)) == NULL) {
        fprintf(stderr, "Memory could not be allocated.");
        exit(1);
    }
    if ((cell.rxyz = (double (*)[3]) malloc(sizeof(double[3])*n)) == NULL ) {
        fprintf(stderr, "Memory could not be allocated.");
        exit(1);
    }
    if ((cell.symb = (char **) malloc(sizeof(char *)*ntyp)) == NULL) {
        fprintf(stderr, "Memory could not be allocated.");
        exit(1);
    }
    if ((cell.rcov = (double *) malloc(sizeof(double)*n)) == NULL) {
        fprintf(stderr, "Memory could not be allocated.");
        exit(1);
    }

    if ((cell.sfp = (double **)malloc(sizeof(double)*n)) == NULL) {
        fprintf(stderr, "Memory could not be allocated.");
        exit(1);
    }
    for ( i = 0; i < n; i++ ) {
        if ((cell.sfp[i] == (double *)malloc(sizeof(double)*l*(ntyp+1))) == NULL) {
            fprintf(stderr, "Memory could not be allocated.");
            exit(1);
        }
    }

    if ((cell.lfp = (double **)malloc(sizeof(double)*n)) == NULL) {
        fprintf(stderr, "Memory could not be allocated.");
        exit(1);
    }
    for ( i = 0; i < n; i++ ) {
        if ((cell.lfp[i] == (double *)malloc(sizeof(double)*lseg)) == NULL) {
            fprintf(stderr, "Memory could not be allocated.");
            exit(1);
        }
    }
    return cell;
}

void c_del_cell(Cell * cell) {
    free(cell->types);
    free(cell->rxyz);
    free(cell->symb);
    free(cell->sfp);
    free(cell->lfp);
    free(cell->rcov);
}

void c_set_cell(Cell * cell, double lat[3][3], double rxyz[][3], int types[], 
        char *symb[], double e = 0) {
    int i, j;

    cell->e = e;
    cpmat3(cell->lat, lat);
    for (i = 0; i < cell->n; i++) {
        for (j = 0; j < 3; j++) {
            cell->rxyz[i][j] = rxyz[i][j];
        }
        cell->types[i] = types[i];
        cell->rocv[i] = get_rcov( symb[types[i]] );
    }

    for ( i = 0; i < cell->ntyp; i++) {
        cell->symb[i] = symb[i];
    }
}


