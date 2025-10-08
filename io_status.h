#ifndef IO_STATUS_H
#define IO_STATUS_H
#define EPS (1e-16)
#define LEN (1234)
#define MAX(a, b) ((a > b) ? a : b)
#define ABS(x) ((x) > 0 ? (x) : (-x))
#define CMP(x, y) \
	((ABS(((x) - (y))) <= ABS(((0.5) * (EPS) * ((x) + (y))))) ? 1 : 0)
#define MAX_PRINT 10
#define SWAP(x, y)        \
	{                     \
		double _ = x;     \
		x = (double &&)y; \
		y = _;            \
	}
typedef enum io_status_
{
	SUCCESS,
	ERROR_OPEN,
	ERROR_READ
} io_status;
#endif
