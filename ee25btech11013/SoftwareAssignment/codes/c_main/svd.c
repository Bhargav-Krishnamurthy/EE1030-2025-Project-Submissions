#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "svd.h"

// sort singular values descending
typedef struct{ 
    double s; // s <- singular value
    int index; 
} SPair;

int cmp(const void *a, const void *b) {
    
    // sort singular values descending
    SPair *pa = (SPair *)a;
    SPair *pb = (SPair *)b;
    
    double da = pa->s;
    double db = pb->s;

    if (da < db) return 1; 
    if (da > db) return -1;
    return 0;
}

void qr_decomposition(double **A, double **Q, double **R, int n) {
    // R and Q zero
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            R[i][j] = 0.0;
            Q[i][j] = 0.0;
        }
    }


    for (int col=0; col<n; col++) {
        // copy col of A to a
        double *a = (double*)malloc(n * sizeof(double));
        for (int i=0;i<n;i++) {
            a[i] = A[i][col];
        }

        // gram-schmidt
        // proj of a on q
        // subtract from a
        for (int c=0; c<col; c++) {
            double mult = 0.0;

            for (int r=0;r<n;r++){
                mult = mult + (Q[r][c] * a[r]);
            }
            R[c][col] = mult;
            for (int r=0;r<n;r++) {
                a[r] = a[r] - (mult * Q[r][c]);
            }
        }
        double square = 0.0;

        for (int i=0;i<n;i++) {
            square += a[i]*a[i];
        }
        double norm = sqrt(square);

        /* normalise col to orthonormal basis */
        if (norm < 1e-14) {
            for (int i=0;i<n;i++) {
                Q[i][col] = 0.0;
            }
            R[col][col] = 0.0;
        } else {
            for (int i=0;i<n;i++) {
                Q[i][col] = a[i] / norm;
            }
            R[col][col] = norm;
        }
        free(a);
    }
}

void qr_iteration(double **A, int n, int iterations, double **eig_vec) {
    double **Q = memoAllocmatrix(n,n);
    double **R = memoAllocmatrix(n,n);

    // eig_vec = I
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            if (row == col) {
                eig_vec[row][col] = 1.0;   
            } else {
                eig_vec[row][col] = 0.0;  
            }
        }
    }


    for (int count=0; count<iterations; count++) {
        // A = Q * R
        qr_decomposition(A, Q, R, n);

        // Anext = R * Q
        double **Anext = matrix_multiplication(R, n, n, Q, n);

        // eig_vec = eig_vec * Q
        double **eig_vec_temp = matrix_multiplication(eig_vec, n, n, Q, n);

        // eig_vec = eig_vec_temp
        // A = Anext
        for (int i=0;i<n;i++){
            for (int j=0;j<n;j++){
                A[i][j] = Anext[i][j];
                eig_vec[i][j] = eig_vec_temp[i][j];
            }
        }

        destroy_matrix(eig_vec_temp, n);
        destroy_matrix(Anext, n);
    }

    destroy_matrix(Q, n);
    destroy_matrix(R, n);
}

void svd(double **A, int m, int n, const double *v, double sigma, double *u) {
    if (sigma < 1e-14) {
        for(int i=0;i<m;i++) {
            u[i] = 0.0;
        }
        return;
    }
    for (int i=0;i<m;i++){
        double x = 0.0;
        for (int j=0;j<n;j++) {
            x += v[j] * A[i][j];
        }
        u[i] = x / sigma;
    }
    
    double square = 0.0;
    for (int i=0;i<m;i++) {
        square += u[i]*u[i];
    }
    double norm = sqrt(square);
    // quick check to avoid div by zero
    if (norm > 1e-14) {
        for (int i=0;i<m;i++) {
            u[i] /= norm;
        }   
    }
}

void truncatedsvd(double **A, int m, int n, int k, double **A_out) {

    
    double **AT = transpose(A, m, n);
    double **ATA = matrix_multiplication(AT, n, m, A, n);

    // QR iteration
    double **ATA_copy = matrix_copy(ATA, n, n);
    double **V = memoAllocmatrix(n, n);
    qr_iteration(ATA_copy, n, QR_ITERATIONS, V);

    // Singular values
    int r = m;
    if(n<m){
        r = n;
    }
    SPair *sv_pairs = malloc(r * sizeof(SPair));
    for (int i=0; i<r; i++) {
        double eig = ATA_copy[i][i];
        if (eig < 1e-14) {
            eig = 0.0;
        }
        sv_pairs[i].index = i;
        sv_pairs[i].s = sqrt(eig);
    }
    qsort(sv_pairs, r, sizeof(SPair), cmp);

    
    double **V_ordered = memoAllocmatrix(n, r);
    double **U_sorted = memoAllocmatrix(m, r);
    double *vcol = malloc(n * sizeof(double));
    double *ucol = malloc(m * sizeof(double));
    for (int kk=0; kk<r; kk++) {
        int index = sv_pairs[kk].index;
        for(int i=0; i<n; i++){
            V_ordered[i][kk] = V[i][index];
            vcol[i] = V[i][index];
        }
        svd(A, m, n, vcol, sv_pairs[kk].s, ucol);
        for (int i=0; i<m; i++) {
            U_sorted[i][kk] = ucol[i];
        }
    }
    free(vcol);
    free(ucol);

    
    if (k > r){
        k = r;
    }
    double **U_k = memoAllocmatrix(m, k);
    double **Sigma_k = memoAllocmatrix(k, k);
    double **V_kT = memoAllocmatrix(k, n);

    for (int i=0; i<m; i++)
        for (int j=0; j<k; j++)
            U_k[i][j] = U_sorted[i][j];

    for (int i=0; i<k; i++)
        Sigma_k[i][i] = sv_pairs[i].s;

    for (int i=0; i<k; i++)
        for (int j=0; j<n; j++)
            V_kT[i][j] = V_ordered[j][i];

    // A_k = U_k * Sigma_k * V_kT
    double **temp = matrix_multiplication(U_k, m, k, Sigma_k, k);
    double **A_k = matrix_multiplication(temp, m, k, V_kT, n);


    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            A_out[i][j] = A_k[i][j];


    destroy_matrix(U_k, m);
    destroy_matrix(Sigma_k, k);
    destroy_matrix(V_kT, k);
    destroy_matrix(temp, m);
    destroy_matrix(A_k, m);
    destroy_matrix(V_ordered, n);
    destroy_matrix(U_sorted, m);
    free(sv_pairs);
    destroy_matrix(V, n);
    destroy_matrix(ATA_copy, n);
    destroy_matrix(ATA, n);
    destroy_matrix(AT, n);
}

