/*******************************************************************************
 * 
 * apc.h
 * This file is part of fplib.
 *
 * SOLUTION OF THE LINEAR MIN-SUM ASSIGNMENT PROBLEM.  
 * HUNGARIAN METHOD. COMPLEXITY O(n^3).
 *
 * October 2008:
 * - Original FORTRAN code translated in C LANGUAGE by Andrea Tramontani 
 *   Andrea Tramontani 
 *   DEIS, University of Bologna 
 *   Viale Risorgimento, 2 
 *   40136 - Bologna (Italy)
 *
 * October 2015:
 * - motified by Li Zhu <z@zhuli.name>
 *
 * ****************************************************************************/

int apc(int n, double *a, double *z_p, int *f);

