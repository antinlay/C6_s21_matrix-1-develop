GCC = gcc -Wall -Werror -Wextra -g
SOURCE = s21_matrix.c
OBJ = $(SOURCE:.c=.o)
GCOVFLAGS =-L. --coverage
HTML = lcov -t "test" -o rep.info -c -d .
OS = $(shell uname)

ifeq ($(OS), Darwin)
	CHECK_H=-lcheck
else
	CHECK_H=-lcheck -pthread -lrt -lm
endif

all: gcov_report

clean:
	rm -rf *.o *.a *.so *.gcda *.gcno *.gch rep.info *.html *.css test report *.txt *.dSYM

test: s21_matrix.a
	$(GCC) matrix_test.c s21_matrix.a -o test $(CHECK_H) $(GCOVFLAGS)
	./test

s21_matrix.a: clean $(OBJ)
	ar rcs s21_matrix.a $(OBJ)
	cp s21_matrix.a libs21_matrix.a
	ranlib s21_matrix.a
	ranlib libs21_matrix.a

gcov_report: test
	$(GCC) matrix_test.c $(SOURCE) -o test $(CHECK_H) $(GCOVFLAGS)
	./test
	$(HTML)
	genhtml -o report rep.info
	open report/index.html

check:
	cppcheck --enable=all --suppress=missingIncludeSystem --inconclusive --check-config $(SOURCE) *.h
	cp ../materials/linters/CPPLINT.cfg CPPLINT.cfg
	python3 ../materials/linters/cpplint.py --extensions=c $(SOURCE) *.h
	rm -rf CPPLINT.cfg
	make test
ifeq ($(OS), Darwin)
	leaks --atExit -- test
else
	CK_FORK=no valgrind --vgdb=no --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=RESULT_VALGRIND.txt ./test
endif
