#include <stdlib.h>
#include <math.h>

void pathgraph_lammax(double *y, double *w, int *assign, int *dim, int *depth, double *lamval) {
    /* 
     Compute the smallest lambada value that set all elements in the path
     graph to zero. 
     */
    int i, j;
    double b_cumsum, lam_temp;
    int p = *dim;
    int d = *depth;
    // allocate memory for an array which will contain ||y_{s_i}||^2 for each index i
    double *b = malloc(d * sizeof(double));
    for (i = 0; i < d; i++) {
        b[i] = 0; // initialize each element to be 0.
    }
    for (i = 0; i < p; i++) {
        if (assign[i] != 0) {
            // Compute sum of squares of values of variables assigned to each node
            b[assign[i]-1] += y[i] * y[i];
        }
    }
    *lamval = sqrt(b[0]) / w[0];
    if (d > 1) {
        b_cumsum = b[0];
        for (j = 1; j < d; j++) {
            b_cumsum += b[j];
            lam_temp = sqrt(b_cumsum) / w[j];
            if (lam_temp > *lamval) *lamval = lam_temp;
        }
    }
}



