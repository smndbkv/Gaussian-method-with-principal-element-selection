#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>

#define EPS 1e-16
#define MAX_PRINT 10
void init_matrix(double *a, int n, int s, double *b, double *x_exact);
void print_matrix(double *a, int n, int m, int r = MAX_PRINT);
void print(double *a, int n, int m = -1, int r = -1);
void get_block(int n, int m, int i, int j, double *a, double *&c, int &v, int &h);
void set_block(int n, int m, int i, int j, double *a, double *c, int v, int h);
double f(int s, int n, int i, int j);
double dot_product(double *a, double *b, int n);
double get_r1(double *a, double *x, double *b, int n);
double get_r2(double *x, double *x_exact, int n);
int gauss_method(int n, int m, double *a, double *b, double *x);

#endif