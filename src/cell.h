// cell.h 
// Copyright (C) 2015 Li Zhu
//

typedef struct {
    int n;
    int ntyp;
    double lat[3][3];
    int *types;
    double (*rxyz)[3];
    double **sfp;
    double **lfp;
    double *rcov;
    double e;
    char **symb;
} Cell;

Cell c_new_cell(int n);
void c_del_cell(Cell * cell);
void c_set_cell(Cell * cell, double lat[3][3], double rxyz[][3], int types[]);

