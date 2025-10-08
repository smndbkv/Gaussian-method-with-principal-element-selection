FLAGS = -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
FLAGS_GPROF = -mfpmath=sse -fstack-protector-all -pg -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
FLAGS_GPROF_O3 = -O3 -mfpmath=sse -fstack-protector-all -pg -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format

all: a.out

a.out: main.cpp matrix.cpp matrix.h io_status.h
	g++ $(FLAGS) $^ -o a.out

gprof: main.cpp matrix.cpp matrix.h io_status.h
	g++ $(FLAGS_GPROF) $^ -o gp.out

gprof_O3: main.cpp matrix.cpp matrix.h io_status.h
	g++ $(FLAGS_GPROF_O3) $^ -o gpO3.out

zip:
	zip Dubkov_SA.zip *.h *.cpp Makefile