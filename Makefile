all: clean s21_decimal.a

s21_decimal.a: clean s21_decimal.o
	ar rcs s21_decimal.a s21_decimal.o
	ar rcs s21_decimal s21_decimal.o

s21_decimal_test.o: clean s21_decimal.c s21_decimal.h
	gcc -fprofile-arcs -ftest-coverage -lgcov -fPIC -Wall -Werror -Wextra -c s21_decimal.c -o s21_decimal.o

s21_decimal.o: clean s21_decimal.c s21_decimal.h
	gcc -Wall -Werror -Wextra -c s21_decimal.c -o s21_decimal.o

s21_decimal_test.a: clean s21_decimal_test.o
	ar rcs s21_decimal.a s21_decimal.o
	ar rcs s21_decimal s21_decimal.o

clean:
	rm -rf *.a
	rm -rf *.o
	rm -rf *.out
	rm -rf *.gcda
	rm -rf *.gcno
	rm -rf *.gcov
	rm -rf *.html
	rm -rf s21_decimal


test: clean all
	gcc test.c -L. s21_decimal.a -o tests.out -lcheck_pic -pthread -lrt -lm -lsubunit

gcov_report: clean s21_decimal_test.a
	gcc test.c -fprofile-arcs -lgcov -L. s21_decimal.a -o tests.out -lcheck_pic -pthread -lrt -lm -lsubunit
	./tests.out
	gcovr --html-details example_html.details.html



