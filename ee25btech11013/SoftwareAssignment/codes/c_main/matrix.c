#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

double **memoAllocmatrix(int r, int c) {
    double **M = (double **)malloc(sizeof(double *) * r);
    if(M == NULL){
        printf("Memory alloc failed!!");
        exit(EXIT_FAILURE);
    }
    for(int i=0; i<r; i++) {
        M[i] = (double *)calloc(c, sizeof(double));  
    }
    return M;
}

void destroy_matrix(double **M, int r) {
    for(int i = 0; i < r; i++){
            free(M[i]);
    }    
    free(M);
}

double **transpose(double **A, int r, int c) {
    double **AT = memoAllocmatrix(c, r);
    for(int i=0;i<r;i++){
            for(int j=0;j<c;j++){
                AT[j][i] = A[i][j];
            }
    }
    return AT;
}

double **matrix_multiplication(double **A, int m, int n, double **B, int p) {
    double **C = memoAllocmatrix(m, p);

    for (int i=0;i<m;i++){
        for (int k=0;  k<n;k++){
            for (int j=0; j<p; j++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }   
    }
    return C;
}

double **matrix_copy(double **A, int r, int c) {
    double **B = memoAllocmatrix(r,c);
    for (int i=0;  i<r; i++){
        for (int j=0 ;j<c;j++){
            B[i][j] = A[i][j];
        }
    }
    return B;
}

double **to_double_matrix(unsigned char *data, int width, int height) {
    double **M = memoAllocmatrix(height, width);
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            M[row][col] = (double) data[row * width + col];
        }
    }

    return M;
}