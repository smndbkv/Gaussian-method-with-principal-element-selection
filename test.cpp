#include "matrix.h"
#include <iostream>
#include <stdio.h>

int main(int argc, char **argv)
{
    int n = 0, m = 0, r = 0, s = 0;
    double t1 = 0, t2 = 0, r1 = 0, r2 = 0;
    char *file_name = nullptr;

    if (!((argc == 5 || argc == 6) && sscanf(argv[1], "%d", &n) == 1 && sscanf(argv[2], "%d", &m) == 1 && sscanf(argv[3], "%d", &r) == 1 && sscanf(argv[4], "%d", &s) == 1 && s >= 0 && s <= 4))
    {
        printf("Usage: %s n m r s <filename>\n", argv[0]);
        return -1;
    }
    if (argc == 6)
    {
        file_name = argv[5];
    }
    else if (argc == 5 && s == 0)
    {
        printf("Usage: %s n m r 0 filename\n", argv[0]);
        return -2;
    }
    double *a = new double[n * n], *b = new double[n], *x = new double[n],
           *x_exact = new double[n];

    io_status st = init_matrix(a, n, s, b, x_exact, file_name);
    switch (st)
    {
    case SUCCESS:
        break;
    case ERROR_OPEN:
        printf("Cannot open %s\n", file_name);
        break;
    case ERROR_READ:
        printf("Cannot read %s\n", file_name);
        break;
    }
    if (st != io_status::SUCCESS)
    {
        delete[] a;
        delete[] b;
        delete[] x;
        delete[] x_exact;
        return -3;
    }

    print(a, n);
    double *block = new double[m * m], *block_inv = new double[m * m];

    int v, h;
    get_block(a, n, m, 0, 0, block, v, h);
    printf("v = %d, h = %d\n", v, h);
    print(block, v, h);
    if (inverse(block, v, block_inv))
    {
        printf("qwertewertrewqertrewq\n");
    }
    set_block(a, n, m, 0, 0, block_inv, v, h);
    print(a, n);

    delete[] block;
    delete[] block_inv;
    delete[] a;
    delete[] b;
    delete[] x;
    delete[] x_exact;
    return 0;
}