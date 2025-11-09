#include "matrix.h"
#include <iostream>
#include <stdio.h>
#include "args.h"

#include <fenv.h>
int feenableexcept(int excepts);
int fedisableexcept(int excepts);
int fegetexcept(void);

void *thread_func(void *ptr)
{
  args *Arg = (args *)ptr;
  int m = Arg->m, q = Arg->q;
  gauss_status *err = Arg->err;
  double *c = new double[m * m], *g = new double[m * m], *d = new double[m * m], *f = new double[m * m];
  double t = Arg->get_cpu_time();
  err[q] = gauss_method(Arg, c, g, d, f);
  t = Arg->get_cpu_time() - t;
  Arg->time = t;
  delete[] c;
  delete[] g;
  delete[] d;
  delete[] f;
  return 0;
}

int main(int argc, char **argv)
{
  int n = 0, m = 0, p = 0, r = 0, s = 0;
  double t1 = 0, t2 = 0, r1 = 0, r2 = 0;
  char *file_name = nullptr;

  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

  if (!((argc == 6 || argc == 7) && sscanf(argv[1], "%d", &n) == 1 && sscanf(argv[2], "%d", &m) == 1 && sscanf(argv[3], "%d", &p) == 1 && sscanf(argv[4], "%d", &r) == 1 && sscanf(argv[5], "%d", &s) == 1 && s >= 0 && s <= 4))
  {
    printf("Usage: %s n m p r s <filename>\n", argv[0]);
    return -1;
  }
  if (argc == 7)
  {
    file_name = argv[6];
  }
  else if (argc == 6 && s == 0)
  {
    printf("Usage: %s n m r 0 filename\n", argv[0]);
    return -2;
  }
  if (m > n)
  {
    m = n;
  }
  if (p > n / m)
  {
    p = n / m;
  }
  double *a = new double[n * n], *b = new double[n], *x = new double[n],
         *x_exact = new double[n];

  io_status st = init_matrix(a, n, s, b, x_exact, file_name);
  switch (st)
  {
  case io_status::SUCCESS:
    break;
  case io_status::ERROR_OPEN:
    printf("Cannot open %s\n", file_name);
    break;
  case io_status::ERROR_READ:
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
  printf("Initial matrix:\n");
  print_matrix(a, n, n, r);
  printf("Initial right side:\n");
  print_matrix(b, 1, n, r);
  printf("Initial exact solution:\n");
  print_matrix(x_exact, 1, n, r);

  pthread_t *tid = new pthread_t[p];
  pthread_barrier_t barrier;
  pthread_barrier_init(&barrier, NULL, p);
  gauss_status *err = new gauss_status[p];
  for (int q = 0; q < p; q++)
    err[q] = gauss_status::DONE;
  args *Arg = new args[p];
  int *per = new int[n / m + 1];
  double *min_norm = new double[p];
  int *min_i = new int[p];
  int *min_j = new int[p];
  double **block = new double *[p];
  for (int q = 0; q < p; q++)
  {
    block[q] = new double[m * m];
  }
  for (int q = 0; q < p; q++)
  {
    Arg[q].a = a;
    Arg[q].n = n;
    Arg[q].m = m;
    Arg[q].b = b;
    Arg[q].x = x;
    Arg[q].err = err;
    Arg[q].q = q;
    Arg[q].p = p;
    Arg[q].per = per;
    Arg[q].min_norm = min_norm;
    Arg[q].min_i = min_i;
    Arg[q].min_j = min_j;
    Arg[q].barrier = &barrier;
    Arg[q].block = block;
  }

  t1 = clock();
  struct timespec start_t1, end_t1;
  clock_gettime(CLOCK_MONOTONIC, &start_t1);

  for (int q = 1; q < p; q++)
  {
    if (pthread_create(tid + q, 0, thread_func, Arg + q))
    {
      printf("Cannot create thread %d\n", q);
    }
  }
  thread_func(Arg + 0);

  for (int k = 1; k < p; k++)
  {
    pthread_join(tid[k], 0);
  }
  t1 = (clock() - t1) / CLOCKS_PER_SEC;
  clock_gettime(CLOCK_MONOTONIC, &end_t1);

  t1 = (end_t1.tv_sec - start_t1.tv_sec) * 1e9;
  t1 = (t1 + (end_t1.tv_nsec - start_t1.tv_nsec)) * 1e-9;
  for (int q = 0; q < p; q++)
  {
    switch (err[q])
    {
    case gauss_status::DONE:
      break;
    case gauss_status::ZERO_MATRIX:
      printf("Initial matrix is zero\n");
      break;
    case gauss_status::NOT_APPLICABLE:
      printf("This solution method is not applicable to this matrix\n");
      break;
    }
    if (err[q] != gauss_status::DONE)
    {
      printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = "
             "%d M = %d\n",
             argv[0], 11, -1., -1., 0., 0., s, n, m);
      delete[] a;
      delete[] b;
      delete[] x;
      delete[] x_exact;
      delete[] tid;
      delete[] err;
      delete[] Arg;
      delete[] per;
      delete[] min_norm;
      delete[] min_i;
      delete[] min_j;
      for (int q = 0; q < p; q++)
      {
        delete[] block[q];
      }
      delete[] block;
      return -4;
    }
  }
  // gauss_status g_st = gauss_method(n, m, a, b, x, c, g, d, f, p);
  // switch (g_st)
  // {
  // case gauss_status::DONE:
  //   break;
  // case gauss_status::ZERO_MATRIX:
  //   printf("Initial matrix is zero\n");
  //   break;
  // case gauss_status::NOT_APPLICABLE:
  //   printf("This solution method is not applicable to this matrix\n");
  //   break;
  // }
  // if (g_st != gauss_status::DONE)
  // {
  //   printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = "
  //          "%d M = %d\n",
  //          argv[0], 11, -1., -1., 0., 0., s, n, m);
  //   delete[] a;
  //   delete[] b;
  //   delete[] x;
  //   delete[] x_exact;
  //   delete[] c;
  //   delete[] g;
  //   delete[] d;
  //   delete[] f;
  //   delete[] p;
  //   return -4;
  // }

  printf("Solution vector:\n");
  print_matrix(x, 1, n, r);

  // инизиализируем еще раз чтобы посчитать невязки
  init_matrix(a, n, s, b, x_exact, file_name);

  t2 = clock();
  r1 = get_r1(a, x, b, n);
  r2 = get_r2(x, x_exact, n);
  t2 = (clock() - t2) / CLOCKS_PER_SEC;

  for (int q = 0; q < p; q++)
  {
    printf("CPU time thread %.2lf\n", Arg[q].time);
  }
  printf(
      "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
      argv[0], 11, r1, r2, t1, t2, s, n, m, p);

  delete[] a;
  delete[] b;
  delete[] x;
  delete[] x_exact;
  delete[] tid;
  delete[] err;
  delete[] Arg;
  delete[] per;
  delete[] min_norm;
  delete[] min_i;
  delete[] min_j;
  for (int q = 0; q < p; q++)
  {
    delete[] block[q];
  }
  delete[] block;
  return 0;
}