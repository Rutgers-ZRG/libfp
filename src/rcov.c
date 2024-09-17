/*******************************************************************************
 * Copyright (C) 2015 Li Zhu 
 * All rights reserved. 
 * 
 * rcov.c
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

double get_rcov(int z) {
    double rcov;
    
    if      ( z == 1  )  rcov = 0.37;   /* "H"   */
    else if ( z == 2  )  rcov = 0.32;   /* "He"  */
    else if ( z == 3  )  rcov = 1.34;   /* "Li"  */
    else if ( z == 4  )  rcov = 0.90;   /* "Be"  */
    else if ( z == 5  )  rcov = 0.82;   /* "B"   */
    else if ( z == 6  )  rcov = 0.77;   /* "C"   */
    else if ( z == 7  )  rcov = 0.75;   /* "N"   */
    else if ( z == 8  )  rcov = 0.73;   /* "O"   */
    else if ( z == 9  )  rcov = 0.71;   /* "F"   */
    else if ( z == 10 )  rcov = 0.69;   /* "Ne"  */
    else if ( z == 11 )  rcov = 1.54;   /* "Na"  */
    else if ( z == 12 )  rcov = 1.30;   /* "Mg"  */
    else if ( z == 13 )  rcov = 1.18;   /* "Al"  */
    else if ( z == 14 )  rcov = 1.11;   /* "Si"  */
    else if ( z == 15 )  rcov = 1.06;   /* "P"   */
    else if ( z == 16 )  rcov = 1.02;   /* "S"   */
    else if ( z == 17 )  rcov = 0.99;   /* "Cl"  */
    else if ( z == 18 )  rcov = 0.97;   /* "Ar"  */
    else if ( z == 19 )  rcov = 1.96;   /* "K"   */
    else if ( z == 20 )  rcov = 1.74;   /* "Ca"  */
    else if ( z == 21 )  rcov = 1.44;   /* "Sc"  */
    else if ( z == 22 )  rcov = 1.36;   /* "Ti"  */
    else if ( z == 23 )  rcov = 1.25;   /* "V"   */
    else if ( z == 24 )  rcov = 1.27;   /* "Cr"  */
    else if ( z == 25 )  rcov = 1.39;   /* "Mn"  */
    else if ( z == 26 )  rcov = 1.25;   /* "Fe"  */
    else if ( z == 27 )  rcov = 1.26;   /* "Co"  */
    else if ( z == 28 )  rcov = 1.21;   /* "Ni"  */
    else if ( z == 29 )  rcov = 1.38;   /* "Cu"  */
    else if ( z == 30 )  rcov = 1.31;   /* "Zn"  */
    else if ( z == 31 )  rcov = 1.26;   /* "Ga"  */
    else if ( z == 32 )  rcov = 1.22;   /* "Ge"  */
    else if ( z == 33 )  rcov = 1.19;   /* "As"  */
    else if ( z == 34 )  rcov = 1.16;   /* "Se"  */
    else if ( z == 35 )  rcov = 1.14;   /* "Br"  */
    else if ( z == 36 )  rcov = 1.10;   /* "Kr"  */
    else if ( z == 37 )  rcov = 2.11;   /* "Rb"  */
    else if ( z == 38 )  rcov = 1.92;   /* "Sr"  */
    else if ( z == 39 )  rcov = 1.62;   /* "Y"   */
    else if ( z == 40 )  rcov = 1.48;   /* "Zr"  */
    else if ( z == 41 )  rcov = 1.37;   /* "Nb"  */
    else if ( z == 42 )  rcov = 1.45;   /* "Mo"  */
    else if ( z == 43 )  rcov = 1.56;   /* "Tc"  */
    else if ( z == 44 )  rcov = 1.26;   /* "Ru"  */
    else if ( z == 45 )  rcov = 1.35;   /* "Rh"  */
    else if ( z == 46 )  rcov = 1.31;   /* "Pd"  */
    else if ( z == 47 )  rcov = 1.53;   /* "Ag"  */
    else if ( z == 48 )  rcov = 1.48;   /* "Cd"  */
    else if ( z == 49 )  rcov = 1.44;   /* "In"  */
    else if ( z == 50 )  rcov = 1.41;   /* "Sn"  */
    else if ( z == 51 )  rcov = 1.38;   /* "Sb"  */
    else if ( z == 52 )  rcov = 1.35;   /* "Te"  */
    else if ( z == 53 )  rcov = 1.33;   /* "I"   */
    else if ( z == 54 )  rcov = 1.30;   /* "Xe"  */
    else if ( z == 55 )  rcov = 2.25;   /* "Cs"  */
    else if ( z == 56 )  rcov = 1.98;   /* "Ba"  */
    else if ( z == 57 )  rcov = 1.80;   /* "La"  */
    else if ( z == 58 )  rcov = 1.63;   /* "Ce"  */
    else if ( z == 59 )  rcov = 1.76;   /* "Pr"  */
    else if ( z == 60 )  rcov = 1.74;   /* "Nd"  */
    else if ( z == 61 )  rcov = 1.73;   /* "Pm"  */
    else if ( z == 62 )  rcov = 1.72;   /* "Sm"  */
    else if ( z == 63 )  rcov = 1.68;   /* "Eu"  */
    else if ( z == 64 )  rcov = 1.69;   /* "Gd"  */
    else if ( z == 56 )  rcov = 1.68;   /* "Tb"  */
    else if ( z == 66 )  rcov = 1.67;   /* "Dy"  */
    else if ( z == 67 )  rcov = 1.66;   /* "Ho"  */
    else if ( z == 68 )  rcov = 1.65;   /* "Er"  */
    else if ( z == 69 )  rcov = 1.64;   /* "Tm"  */
    else if ( z == 70 )  rcov = 1.70;   /* "Yb"  */
    else if ( z == 71 )  rcov = 1.60;   /* "Lu"  */
    else if ( z == 72 )  rcov = 1.50;   /* "Hf"  */
    else if ( z == 73 )  rcov = 1.38;   /* "Ta"  */
    else if ( z == 74 )  rcov = 1.46;   /* "W"   */
    else if ( z == 75 )  rcov = 1.59;   /* "Re"  */
    else if ( z == 76 )  rcov = 1.28;   /* "Os"  */
    else if ( z == 77 )  rcov = 1.37;   /* "Ir"  */
    else if ( z == 78 )  rcov = 1.28;   /* "Pt"  */
    else if ( z == 79 )  rcov = 1.44;   /* "Au"  */
    else if ( z == 80 )  rcov = 1.49;   /* "Hg"  */
    else if ( z == 81 )  rcov = 1.48;   /* "Tl"  */
    else if ( z == 82 )  rcov = 1.47;   /* "Pb"  */
    else if ( z == 83 )  rcov = 1.46;   /* "Bi"  */
    else if ( z == 84 )  rcov = 1.45;   /* "Po"  */
    else if ( z == 85 )  rcov = 1.47;   /* "At"  */
    else if ( z == 86 )  rcov = 1.42;   /* "Rn"  */
    else if ( z == 87 )  rcov = 2.23;   /* "Fr"  */
    else if ( z == 88 )  rcov = 2.01;   /* "Ra"  */
    else if ( z == 89 )  rcov = 1.86;   /* "Ac"  */
    else if ( z == 90 )  rcov = 1.75;   /* "Th"  */
    else if ( z == 91 )  rcov = 1.69;   /* "Pa"  */
    else if ( z == 92 )  rcov = 1.70;   /* "U"   */
    else if ( z == 93 )  rcov = 1.71;   /* "Np"  */
    else if ( z == 94 )  rcov = 1.72;   /* "Pu"  */
    else if ( z == 95 )  rcov = 1.66;   /* "Am"  */
    else if ( z == 96 )  rcov = 1.66;   /* "Cm"  */
    else if ( z == 97 )  rcov = 1.68;   /* "Bk"  */
    else if ( z == 98 )  rcov = 1.68;   /* "Cf"  */
    else if ( z == 99 )  rcov = 1.65;   /* "Es"  */
    else if ( z == 100)  rcov = 1.67;   /* "Fm"  */
    else if ( z == 101)  rcov = 1.73;   /* "Md"  */
    else if ( z == 102)  rcov = 1.76;   /* "No"  */
    else if ( z == 103)  rcov = 1.61;   /* "Lr"  */
    else if ( z == 104)  rcov = 1.57;   /* "Rf"  */
    else if ( z == 105)  rcov = 1.49;   /* "Db"  */
    else if ( z == 106)  rcov = 1.43;   /* "Sg"  */
    else if ( z == 107)  rcov = 1.41;   /* "Bh"  */
    else if ( z == 108)  rcov = 1.34;   /* "Hs"  */
    else if ( z == 109)  rcov = 1.29;   /* "Mt"  */
    else if ( z == 110)  rcov = 1.28;   /* "Ds"  */
    else if ( z == 111)  rcov = 1.21;   /* "Rg"  */
    else if ( z == 112)  rcov = 1.22;   /* "Cn"  */
    else rcov = 1.0;

    return rcov;
}
