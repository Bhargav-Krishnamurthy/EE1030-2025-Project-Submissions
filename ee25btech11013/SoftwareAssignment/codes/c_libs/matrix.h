#ifndef MATRIX_H
#define MATRIX_H

double **memoAllocmatrix(int r, int c);
void destroy_matrix(double **M, int r);
double **transpose(double **A, int r, int c);
double **matrix_multiplication(double **A, int m, int n, double **B, int p);
double **matrix_copy(double **A, int r, int c);

double **to_double_matrix(unsigned char *data, int width, int height);

#endif
