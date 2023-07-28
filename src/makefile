CC=g++
SRC=s21_matrix.cpp
OBJ=s21_matrix.o
CFLAGS= -g -Wall -Werror -Wextra -std=c++17
TESTFLAGS=-lgtest

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)    # MacOS
    OPEN_CMD = open
    MEMCHECK_CMD = leaks -atExit -- ./test.out | grep LEAK:
else
    OPEN_CMD = xdg-open
    MEMCHECK_CMD = valgrind --tool=memcheck --leak-check=yes ./test.out
endif

all: gcov_report

s21_matrix_oop.a:
	$(CC) $(CFLAGS) $(SRC) -c
	ar -rcs s21_matrix_oop.a $(OBJ)

test: s21_matrix_oop.a
	$(CC) $(CFLAGS) $(TESTFLAGS) test.cc s21_matrix_oop.a -o test.out 
	./test.out

memcheck: test
	$(MEMCHECK_CMD)

gcov_report: test
	$(CC) test.cc -c
	$(CC) --coverage  s21_matrix.cpp  test.o -o test.out $(TESTFLAGS)
	./test.out
	lcov -t "test" -o test.info -c -d ./
	genhtml -o report test.info
	$(OPEN_CMD) report/index.html

clean:
	rm -rf *.out *.o s21_matrix_oop.a *.gcda *.gcno *.info 
	-rm -rf report

clang:
	clang-format -i *.cpp *.cpp *.h
	clang-format -n *.cpp *.h *.cc
