#ifndef SVD_H
#define SVD_H

#define QR_ITERATIONS 10

void qr_decomposition(double **A, double **Q, double **R, int n);
void qr_iteration(double **A, int n, int iterations, double **eig_vec);

void svd(double **A, int m, int n, const double *v, double sigma, double *u);
void truncatedsvd(double **A, int m, int n, int k, double **A_out);

#endif
