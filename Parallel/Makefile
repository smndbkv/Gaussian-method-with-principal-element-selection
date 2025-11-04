FLAGS = -pthread -std=c++17 -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
FLAGS_GPROF = -pthread -mfpmath=sse -fstack-protector-all -pg -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
FLAGS_GPROF_O2 = -pthread -O2 -mfpmath=sse -fstack-protector-all -pg -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
FLAGS_GPROF_O3 = -pthread -O3 -mfpmath=sse -fstack-protector-all -pg -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
N = 4
M = 1
P = 1
R = 10
TEST_FILES = /home/semyon/BOGACHEVTASK/5_sem/matrix_test/*.txt



all: a.out

a.out: main.cpp matrix.h matrix.cpp io_status.h args.h
	g++ $(FLAGS) $^ -o a.out

gprof: main.cpp matrix.cpp matrix.h io_status.h
	g++ $(FLAGS_GPROF) $^ -o gp.out


gprof_O2: main.cpp matrix.cpp matrix.h io_status.h
	g++ $(FLAGS_GPROF_O2) $^ -o gpO2.out

gprof_O3: main.cpp matrix.cpp matrix.h io_status.h
	g++ $(FLAGS_GPROF_O3) $^ -o gpO3.out

zip:
	zip Dubkov_SA.zip *.h *.cpp Makefile -x "test.cpp"


test: a.out
	@echo "Running tests..."
	@for file in $(TEST_FILES); do \
		echo "Testing with $$file:"; \
		./a.out $(N) $(M) $(P) $(R) 0 $$file; \
		echo "---"; \
	done

.PHONY: test clean