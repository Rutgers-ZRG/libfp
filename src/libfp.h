/*******************************************************************************
 * Copyright (C) 2015 Li Zhu 
 * All rights reserved. 
 * 
 * libfp.h
 * This file is part of libfp.
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

int get_major_version(void);
int get_minor_version(void);
int get_micro_version(void);

void get_fp_nonperiodic(int nid, int nat, int ntyp, int types[], double rxyz[][3], int znucl[], double fp[]);

void get_fp_periodic(int flag, int ldfp, int log, int lmax, int nat, int ntyp, int types[], double lat[3][3],
        double rxyz[][3], int znucl[], int natx, double cutoff, double **sfp, double **lfp, double ****dfp);

double get_fpdistance_periodic(int nat, int ntyp, int types[], int fp_len, 
        double **fp1, double **fp2, int f[]);

