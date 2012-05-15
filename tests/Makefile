CFLAGS = -O2 -Wall
#CFLAGS = -O0 -ggdb -Wall
DUMALIB = /d/dev/libs/output/duma/libduma.a

all: test_la test_poly test_poly_bound test_shapes test_bvh test_bvh_mem test_shapeset
test_la: test_common.h test_la.c ../geom_la.c
	gcc $(CFLAGS) test_common.c test_la.c ../geom_la.c -I.. -I. -o test_la
test_poly: test_common.h test_poly.c ../geom_la.c ../geom_poly.c
	gcc $(CFLAGS) test_common.c test_poly.c ../geom_la.c ../geom_poly.c -I.. -I. -o test_poly
test_poly_bound: test_common.h test_poly_bound.c ../geom_la.c ../geom_poly.c
	gcc $(CFLAGS) test_common.c test_poly_bound.c ../geom_la.c ../geom_poly.c -I.. -I. -o test_poly_bound
test_shapes: test_common.h test_shapes.c ../geom_la.c ../geom_poly.c ../geom_shapes.c
	gcc $(CFLAGS) test_common.c test_shapes.c ../geom_la.c ../geom_poly.c ../geom_shapes.c -I.. -I. $(DUMALIB) -o test_shapes
test_bvh: test_common.h test_bvh.c ../geom_bvh.c
	gcc $(CFLAGS) test_common.c test_bvh.c ../geom_bvh.c -I.. -I. $(DUMALIB) -o test_bvh
test_bvh_mem: test_common.h test_bvh_mem.c ../geom_bvh.c
	gcc $(CFLAGS) test_common.c test_bvh_mem.c ../geom_bvh.c -I.. -I. $(DUMALIB) -o test_bvh_mem
test_shapeset: test_common.h test_shapeset.c ../geom_bvh.c ../geom_shapes.c ../geom_la.c ../geom_shapeset.c ../geom_poly.c
	gcc $(CFLAGS) test_common.c test_shapeset.c ../geom_bvh.c ../geom_shapes.c ../geom_la.c ../geom_shapeset.c ../geom_poly.c -I.. -I. $(DUMALIB) -o test_shapeset
