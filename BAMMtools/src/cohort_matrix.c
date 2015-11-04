#include <R.h>

void cohort_matrix(int* tipStates, int* Nsamples, int* Ntips, double* mat);

void cohort_matrix(int* tipStates, int* Nsamples, int* Ntips, double* mat) {
    int i, j, k, ntips = *Ntips;    
    double f, n = (double)(*Nsamples);    
    for (i = 0; i < *Nsamples; i++) {
        for (j = 0; j < *Ntips; j++) {
            for (k = j+1; k < *Ntips; k++) {
                f = (double)(tipStates[j + i*ntips] == tipStates[k + i*ntips]);
                mat[j + k*ntips] = mat[k + j*ntips] += f / n; 
            }
        }
    }
}
