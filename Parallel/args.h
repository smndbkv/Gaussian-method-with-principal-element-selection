#ifndef ARGS_H
#define ARGS_H
#include <pthread.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "io_status.h"
class args
{
public:
    double *a = nullptr; // матрица
    int n = 0;           // её размер
    int m = 0;           // размер блока
    double *b = nullptr;
    double *x = nullptr;
    gauss_status *err = nullptr;
    int p = 0; // число потоков
    int q = 0; // номер потока
    int *per = nullptr;
    double *min_norm = nullptr;
    int *min_i = nullptr;
    int *min_j = nullptr;
    pthread_barrier_t *barrier;
    double **block = nullptr;
    double time = 0;
    void reduce_sum(double *a = nullptr, int n = 0)
    {
        static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
        static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
        static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
        static int t_in = 0;
        static int t_out = 0;
        static double *r = nullptr;
        int i;

        if (p <= 1)
            return;

        pthread_mutex_lock(&m);

        if (r == nullptr)
        {
            r = a;
        }
        else
        {
            for (i = 0; i < n; i++)
            {
                r[i] += a[i];
            }
        }

        t_in++;
        if (t_in >= p)
        {
            t_out = 0;
            pthread_cond_broadcast(&c_in);
        }
        else
        {
            while (t_in < p)
            {
                pthread_cond_wait(&c_in, &m);
            }
        }

        if (a != r)
        {
            for (i = 0; i < n; i++)
            {
                a[i] = r[i];
            }
        }

        t_out++;
        if (t_out >= p)
        {
            t_in = 0;
            r = nullptr;
            pthread_cond_broadcast(&c_out);
        }
        else
        {
            while (t_out < p)
            {
                pthread_cond_wait(&c_out, &m);
            }
        }

        pthread_mutex_unlock(&m);
        return;
    }
    double get_full_time()
    {
        struct timeval t;
        gettimeofday(&t, 0);
        return t.tv_sec + t.tv_usec / 1e6;
    }
    double get_cpu_time()
    {
        struct rusage t;
        getrusage(RUSAGE_THREAD, &t);
        return t.ru_utime.tv_sec + t.ru_utime.tv_usec / 1e6;
    }
};

#endif