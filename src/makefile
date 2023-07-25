CC=g++
SRC=s21_matrix.cpp
OBJ=s21_matrix.o
CFLAGS= -g -Wall -Werror -Wextra -std=c++17
TESTFLAGS=-lgtest

all: gcov_report


s21_matrix_oop.a:
	$(CC) $(CFLAGS) $(SRC) -c
	ar -rcs s21_matrix_oop.a $(OBJ)

test: s21_matrix_oop.a
	$(CC) $(CFLAGS) $(TESTFLAGS) test.cc s21_matrix_oop.a -o test.out 
	./test.out

gcov_report:
	$(CC) test.cc -c
	$(CC) --coverage  s21_matrix.cpp  test.o -o test.out $(TESTFLAGS)
	./test.out
	lcov -t "test" -o test.info -c -d ./
	genhtml -o report test.info
	open report/index.html


clean:
	rm -rf *.out *.o s21_matrix_oop.a *.gcda *.gcno *.info 
	-rm -rf report

clang:
	clang-format -i *.cpp *.cpp *.h
	clang-format -n *.cpp *.h *.cc