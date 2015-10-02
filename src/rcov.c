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
#include <stdlib.h>
#include <string.h>

double get_rcov(char * s) {
    double rcov;
    
    if      (strcmp(s, "H" ) == 0)  rcov = 0.37;
    else if (strcmp(s, "He") == 0)  rcov = 0.32;
    else if (strcmp(s, "Li") == 0)  rcov = 1.34;
    else if (strcmp(s, "Be") == 0)  rcov = 0.90;
    else if (strcmp(s, "B" ) == 0)  rcov = 0.82;
    else if (strcmp(s, "C" ) == 0)  rcov = 0.77;
    else if (strcmp(s, "N" ) == 0)  rcov = 0.75;
    else if (strcmp(s, "O" ) == 0)  rcov = 0.73;
    else if (strcmp(s, "F" ) == 0)  rcov = 0.71;
    else if (strcmp(s, "Ne") == 0)  rcov = 0.69;
    else if (strcmp(s, "Na") == 0)  rcov = 1.54;
    else if (strcmp(s, "Mg") == 0)  rcov = 1.30;
    else if (strcmp(s, "Al") == 0)  rcov = 1.18;
    else if (strcmp(s, "Si") == 0)  rcov = 1.11;
    else if (strcmp(s, "P" ) == 0)  rcov = 1.06;
    else if (strcmp(s, "S" ) == 0)  rcov = 1.02;
    else if (strcmp(s, "Cl") == 0)  rcov = 0.99;
    else if (strcmp(s, "Ar") == 0)  rcov = 0.97;
    else if (strcmp(s, "K" ) == 0)  rcov = 1.96;
    else if (strcmp(s, "Ca") == 0)  rcov = 1.74;
    else if (strcmp(s, "Sc") == 0)  rcov = 1.44;
    else if (strcmp(s, "Ti") == 0)  rcov = 1.36;
    else if (strcmp(s, "V" ) == 0)  rcov = 1.25;
    else if (strcmp(s, "Cr") == 0)  rcov = 1.27;
    else if (strcmp(s, "Mn") == 0)  rcov = 1.39;
    else if (strcmp(s, "Fe") == 0)  rcov = 1.25;
    else if (strcmp(s, "Co") == 0)  rcov = 1.26;
    else if (strcmp(s, "Ni") == 0)  rcov = 1.21;
    else if (strcmp(s, "Cu") == 0)  rcov = 1.38;
    else if (strcmp(s, "Zn") == 0)  rcov = 1.31;
    else if (strcmp(s, "Ga") == 0)  rcov = 1.26;
    else if (strcmp(s, "Ge") == 0)  rcov = 1.22;
    else if (strcmp(s, "As") == 0)  rcov = 1.19;
    else if (strcmp(s, "Se") == 0)  rcov = 1.16;
    else if (strcmp(s, "Br") == 0)  rcov = 1.14;
    else if (strcmp(s, "Kr") == 0)  rcov = 1.10;
    else if (strcmp(s, "Rb") == 0)  rcov = 2.11;
    else if (strcmp(s, "Sr") == 0)  rcov = 1.92;
    else if (strcmp(s, "Y" ) == 0)  rcov = 1.62;
    else if (strcmp(s, "Zr") == 0)  rcov = 1.48;
    else if (strcmp(s, "Nb") == 0)  rcov = 1.37;
    else if (strcmp(s, "Mo") == 0)  rcov = 1.45;
    else if (strcmp(s, "Tc") == 0)  rcov = 1.56;
    else if (strcmp(s, "Ru") == 0)  rcov = 1.26;
    else if (strcmp(s, "Rh") == 0)  rcov = 1.35;
    else if (strcmp(s, "Pd") == 0)  rcov = 1.31;
    else if (strcmp(s, "Ag") == 0)  rcov = 1.53;
    else if (strcmp(s, "Cd") == 0)  rcov = 1.48;
    else if (strcmp(s, "In") == 0)  rcov = 1.44;
    else if (strcmp(s, "Sn") == 0)  rcov = 1.41;
    else if (strcmp(s, "Sb") == 0)  rcov = 1.38;
    else if (strcmp(s, "Te") == 0)  rcov = 1.35;
    else if (strcmp(s, "I" ) == 0)  rcov = 1.33;
    else if (strcmp(s, "Xe") == 0)  rcov = 1.30;
    else if (strcmp(s, "Cs") == 0)  rcov = 2.25;
    else if (strcmp(s, "Ba") == 0)  rcov = 1.98;
    else if (strcmp(s, "La") == 0)  rcov = 1.80;
    else if (strcmp(s, "Ce") == 0)  rcov = 1.63;
    else if (strcmp(s, "Pr") == 0)  rcov = 1.76;
    else if (strcmp(s, "Nd") == 0)  rcov = 1.74;
    else if (strcmp(s, "Pm") == 0)  rcov = 1.73;
    else if (strcmp(s, "Sm") == 0)  rcov = 1.72;
    else if (strcmp(s, "Eu") == 0)  rcov = 1.68;
    else if (strcmp(s, "Gd") == 0)  rcov = 1.69;
    else if (strcmp(s, "Tb") == 0)  rcov = 1.68;
    else if (strcmp(s, "Dy") == 0)  rcov = 1.67;
    else if (strcmp(s, "Ho") == 0)  rcov = 1.66;
    else if (strcmp(s, "Er") == 0)  rcov = 1.65;
    else if (strcmp(s, "Tm") == 0)  rcov = 1.64;
    else if (strcmp(s, "Yb") == 0)  rcov = 1.70;
    else if (strcmp(s, "Lu") == 0)  rcov = 1.60;
    else if (strcmp(s, "Hf") == 0)  rcov = 1.50;
    else if (strcmp(s, "Ta") == 0)  rcov = 1.38;
    else if (strcmp(s, "W" ) == 0)  rcov = 1.46;
    else if (strcmp(s, "Re") == 0)  rcov = 1.59;
    else if (strcmp(s, "Os") == 0)  rcov = 1.28;
    else if (strcmp(s, "Ir") == 0)  rcov = 1.37;
    else if (strcmp(s, "Pt") == 0)  rcov = 1.28;
    else if (strcmp(s, "Au") == 0)  rcov = 1.44;
    else if (strcmp(s, "Hg") == 0)  rcov = 1.49;
    else if (strcmp(s, "Tl") == 0)  rcov = 1.48;
    else if (strcmp(s, "Pb") == 0)  rcov = 1.47;
    else if (strcmp(s, "Bi") == 0)  rcov = 1.46;
    else if (strcmp(s, "Po") == 0)  rcov = 1.45;
    else if (strcmp(s, "At") == 0)  rcov = 1.47;
    else if (strcmp(s, "Rn") == 0)  rcov = 1.42;
    else if (strcmp(s, "Fr") == 0)  rcov = 2.23;
    else if (strcmp(s, "Ra") == 0)  rcov = 2.01;
    else if (strcmp(s, "Ac") == 0)  rcov = 1.86;
    else if (strcmp(s, "Th") == 0)  rcov = 1.75;
    else if (strcmp(s, "Pa") == 0)  rcov = 1.69;
    else if (strcmp(s, "U" ) == 0)  rcov = 1.70;
    else if (strcmp(s, "Np") == 0)  rcov = 1.71;
    else if (strcmp(s, "Pu") == 0)  rcov = 1.72;
    else if (strcmp(s, "Am") == 0)  rcov = 1.66;
    else if (strcmp(s, "Cm") == 0)  rcov = 1.66;
    else if (strcmp(s, "Bk") == 0)  rcov = 1.68;
    else if (strcmp(s, "Cf") == 0)  rcov = 1.68;
    else if (strcmp(s, "Es") == 0)  rcov = 1.65;
    else if (strcmp(s, "Fm") == 0)  rcov = 1.67;
    else if (strcmp(s, "Md") == 0)  rcov = 1.73;
    else if (strcmp(s, "No") == 0)  rcov = 1.76;
    else if (strcmp(s, "Lr") == 0)  rcov = 1.61;
    else if (strcmp(s, "Rf") == 0)  rcov = 1.57;
    else if (strcmp(s, "Db") == 0)  rcov = 1.49;
    else if (strcmp(s, "Sg") == 0)  rcov = 1.43;
    else if (strcmp(s, "Bh") == 0)  rcov = 1.41;
    else if (strcmp(s, "Hs") == 0)  rcov = 1.34;
    else if (strcmp(s, "Mt") == 0)  rcov = 1.29;
    else if (strcmp(s, "Ds") == 0)  rcov = 1.28;
    else if (strcmp(s, "Rg") == 0)  rcov = 1.21;
    else if (strcmp(s, "Cn") == 0)  rcov = 1.22;
    else rcov = 1.0;

    return rcov /  0.52917720859;
}
