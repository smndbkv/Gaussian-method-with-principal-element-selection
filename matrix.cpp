#include "matrix.h"

void init_matrix(double *a, int n, int s, double *b, double *x_exact)
{
    if (s != 0)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                a[i * n + j] = f(s, n, i + 1, j + 1);
            }
        }
        int p = (n - 1) / 2;
        for (int i = 0; i < n; i++)
        {
            b[i] = 0;
            x_exact[i] = (i + 1) % 2;
            for (int j = 0; j <= p; j++)
            {
                b[i] += a[i * n + 2 * j];
            }
        }
    }
}

void get_block(int n, int m, int i, int j, double *a, double *&c, int &v, int &h)
{
    int k = n / m, l = n % m;
    v = (i < k ? m : l);
    h = (j < k ? m : l);
    c = a + i * n * m + j * m;
}

void set_block(int n, int m, int i, int j, double *a, double *c, int v, int h)
{
    for (int r = 0; r < v; r++)
    {
        for (int s = 0; s < h; s++)
        {
            a[(i * m + r) * n + (j * m + s)] = c[r * h + s];
        }
    }
}

void print_matrix(double *a, int n, int m, int r)
{
    int nn = std::min(n, r), mm = std::min(m, r);
    for (int i = 0; i < nn; i++)
    {
        for (int j = 0; j < mm; j++)
        {
            printf(" %10.3e", a[i * m + j]);
        }
        printf("\n");
    }
}
void print(double *a, int n, int m, int r)
{
    if (m == -1)
    {
        m = n;
    }
    if (r == -1)
    {
        r = n;
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            printf("%lf ", a[i * r + j]);
        }
        printf("\n");
    }
}
double f(int s, int n, int i, int j)
{
    switch (s)
    {
    case 1:
        return n - std::max(i, j) + 1;
        break;
    case 2:
        return std::max(i, j);
        break;
    case 3:
        return std::abs(i - j);
        break;
    case 4:
        return 1. / (i + j - 1);
        break;
    }
    return 0;
}

double dot_product(double *a, double *b, int n)
{
    double s = 0;
    for (int i = 0; i < n; i++)
    {
        s += a[i] * b[i];
    }
    return s;
}
double get_r1(double *a, double *x, double *b, int n)
{
    double s1 = 0, s2 = 0;
    for (int i = 0; i < n; i++)
    {
        s1 += std::abs(dot_product(a + i * n, x, n) - b[i]);
        s2 += std::abs(b[i]);
    }
    if (std::abs(s2) < EPS)
    {
        return -1;
    }
    return s1 / s2;
}

double get_r2(double *x, double *x_exact, int n)
{
    double s1 = 0, s2 = 0;
    for (int i = 0; i < n; i++)
    {
        s1 += std::abs(x[i] - x_exact[i]);
        s2 += std::abs(x_exact[i]);
    }
    if (std::abs(s2) < EPS)
    {
        return -1;
    }
    return s1 / s2;
}

// int gauss_method(int n, int m, double *a, double *b, double *x)
// {
//     return -1;
// }
