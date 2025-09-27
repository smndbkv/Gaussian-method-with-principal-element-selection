#include <iostream>
#include <stdio.h>
#include "matrix.h"

int main(int argc, char **argv)
{
    int n = 0, m = 0, r = 0, s = 0;
    double t1 = 0, t2 = 0, r1 = 0, r2 = 0;
    // char *filename = nullptr;
    //  double time = 0;
    if (!((argc == 5 || argc == 6) &&
          sscanf(argv[1], "%d", &n) == 1 &&
          sscanf(argv[2], "%d", &m) == 1 &&
          sscanf(argv[3], "%d", &r) == 1 &&
          sscanf(argv[4], "%d", &s) == 1 &&
          s >= 0 && s <= 4))
    {
        printf("Usage: %s n m r s <filename>\n", argv[0]);
        return 0;
    }
    // if (argc == 6)
    // {
    //     filename = argv[5];
    // }
    double *a = new double[n * n], *b = new double[n], *x = new double[n], *x_exact = new double[n];

    init_matrix(a, n, s, b, x_exact);
    printf("Initial matrix:\n");
    print_matrix(a, n, n, r);
    printf("Initial right side:\n");
    print_matrix(b, 1, n, r);
    printf("Initial exact solution:\n");
    print_matrix(x_exact, 1, n, r);

    t1 = clock();
    // gauss_method(n, m, a, b, x);
    for (int i = 0; i < n; i++)
    {
        x[i] = (i + 1) % 2;
    }
    t1 = (clock() - t1) / CLOCKS_PER_SEC;

    t2 = clock();
    r1 = get_r1(a, x, b, n);
    r2 = get_r2(x, x_exact, n);
    t2 = (clock() - t2) / CLOCKS_PER_SEC;

    printf(
        "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n",
        argv[0], 11, r1, r2, t1, t2, s, n, m);

    delete[] a;
    delete[] b;
    delete[] x;
    delete[] x_exact;
    return 0;
}