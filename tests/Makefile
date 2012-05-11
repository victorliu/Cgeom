CFLAGS = -O2 -Wall

all: test_la test_poly
test_la: test_common.h test_la.c
	gcc $(CFLAGS) test_common.c test_la.c ../geom_la.c -I.. -I. -o test_la
test_poly: test_common.h test_poly.c ../geom_poly.c
	gcc $(CFLAGS) test_common.c test_poly.c ../geom_la.c ../geom_poly.c -I.. -I. -o test_poly
