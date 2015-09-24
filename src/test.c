#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fplib.h"


int main() {
    char ch, buf[72];
    FILE *fp;
    double sc, lat[3][3], (*pos)[3], (*rxyz)[3];
    double aa = 0.52917720859, t;
    int i, j, k, iconf, *typt;
    char sbuf[2]="Ab";
    char **symb, *orb="sp";

    int nconf, nat, ntyp, lseg=4, l=2, natx=100;
    double cutoff, **sfp, **lfp;
    int types[]={1,1,1,1,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3};

    scanf("%d", &nconf);
    scanf("%d", &nat);
    scanf("%d", &ntyp);
    scanf("%lf", &cutoff);

    sfp = (double **) malloc(sizeof(double)*nat);
    lfp = (double **) malloc(sizeof(double)*nat);

    for (i = 0; i < nat; i++){
        sfp[i] = (double *) malloc(sizeof(double)*(l+1)*(ntyp+1));
        lfp[i] = (double *) malloc(sizeof(double)*(natx*lseg));
    }

    symb = malloc(ntyp * sizeof(char*));
    for (i = 0; i < ntyp; i++)
        symb[i] = malloc(2*sizeof(char));
    typt = (int *) malloc(ntyp * sizeof(int));
    pos  = (double (*)[3]) malloc(nat * sizeof(double[3]));
    rxyz = (double (*)[3]) malloc(nat * sizeof(double[3]));


    fp = fopen("POSCARs", "r");

    if (fp == NULL) {
        perror("Error while opening the file.\n");
        exit(EXIT_FAILURE);
    }

    for (iconf = 0; iconf < nconf; iconf++) {

        fscanf(fp, "%s", buf);
        fscanf(fp, "%lf", &sc);
        for (i=0; i< 3; i++){
            fscanf(fp, "%lf %lf %lf", &lat[i][0], &lat[i][1], &lat[i][2]);
            lat[i][0] *= (sc/aa);
            lat[i][1] *= (sc/aa);
            lat[i][2] *= (sc/aa);
            printf("%15.9f %15.9f %15.9f\n", lat[i][0], lat[i][1], lat[i][2]);
        }

        for (i = 0; i < ntyp; i++){
            fscanf(fp, "%s", sbuf);
            strcpy(symb[i], sbuf);
        }

        for (i = 0; i < ntyp; i++)
            fscanf(fp, "%d", &typt[i]);

        fscanf(fp, "%s", buf);

        for (i = 0; i< nat; i++){
            fscanf(fp, "%le %le %le", &pos[i][0], &pos[i][1], &pos[i][2]);
            for (j = 0; j < 3; j++){
                if ( pos[i][j] >= 1.0 ){
                    pos[i][j] -= (int)(pos[i][j]);
                }
                if ( pos[i][j] < 0.) {
                    pos[i][j] -= (int)(pos[i][j]-1);
                }
            }
        }

        for (i = 0; i< nat; i++) {
            for (j = 0; j < 3; j++){
                t = 0.0;
                for (k = 0; k < 3; k++){ 
                    t += pos[i][k] * lat[k][j];
                }
                rxyz[i][j] = t;
            }

        }

       // for (i = 0; i< nat; i++)
       //         printf("%15.9f %15.9f %15.9f\n",  rxyz[i][0], rxyz[i][1], rxyz[i][2]);
        get_fingerprint_periodic(nat, ntyp, types, lat,
                rxyz, symb, natx, orb,  cutoff,
            sfp, lfp);


    }

    free(symb);
    free(pos);
    free(rxyz);
    return 0;
}




