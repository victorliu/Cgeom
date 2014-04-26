CFLAGS = -O2 -Wall
#CFLAGS = -O0 -ggdb -Wall
DUMALIB = /d/dev/libs/output/duma/libduma.a

LIBS = -lm

all: test_la test_poly test_poly_bound test_shapes test_bvh test_bvh_mem test_shapeset
test_la: test_common.h test_la.c ../geom_la.c
	gcc $(CFLAGS) test_common.c test_la.c ../geom_la.c -I.. -I. -o test_la $(LIBS)
test_poly: test_common.h test_poly.c ../geom_la.c ../geom_poly.c
	gcc $(CFLAGS) test_common.c test_poly.c ../geom_la.c ../geom_poly.c -I.. -I. -o test_poly $(LIBS)
test_poly_bound: test_common.h test_poly_bound.c ../geom_la.c ../geom_poly.c
	gcc $(CFLAGS) test_common.c test_poly_bound.c ../geom_la.c ../geom_poly.c -I.. -I. -o test_poly_bound $(LIBS)
test_shapes: test_common.h test_shapes.c ../geom_la.c ../geom_poly.c ../geom_shapes.c
	gcc $(CFLAGS) test_common.c test_shapes.c ../geom_la.c ../geom_poly.c ../geom_shapes.c ../geom_predicates.c -I.. -I. $(DUMALIB) -o test_shapes $(LIBS)
test_bvh: test_common.h test_bvh.c ../geom_bvh.c
	gcc $(CFLAGS) test_common.c test_bvh.c ../geom_bvh.c -I.. -I. $(DUMALIB) -o test_bvh $(LIBS)
test_bvh_mem: test_common.h test_bvh_mem.c ../geom_bvh.c
	gcc $(CFLAGS) test_common.c test_bvh_mem.c ../geom_bvh.c -I.. -I. $(DUMALIB) -o test_bvh_mem $(LIBS)
test_shapeset: test_common.h test_shapeset.c ../geom_bvh.c ../geom_shapes.c ../geom_la.c ../geom_shapeset.c ../geom_poly.c
	gcc $(CFLAGS) test_common.c test_shapeset.c ../geom_bvh.c ../geom_shapes.c ../geom_predicates.c ../geom_la.c ../geom_shapeset.c ../geom_poly.c -I.. -I. $(DUMALIB) -o test_shapeset $(LIBS)
test_arc: test_common.h test_arc.c ../geom_arc.c
	gcc $(CFLAGS) -I../.. test_common.c test_arc.c ../geom_arc.c ../geom_sphereavg.c ../geom_la.c -I.. -I. -o test_arc $(LIBS)
test_arc_unparam: test_common.h test_arc_unparam.c ../geom_arc.c
	gcc $(CFLAGS) -I.. -I../.. -I. test_common.c test_arc_unparam.c ../geom_arc.c ../geom_sphereavg.c ../geom_la.c -o test_arc_unparam $(LIBS)
test_circle_fit: test_common.h test_circle_fit.c ../geom_circum.c
	gcc $(CFLAGS) test_common.c test_circle_fit.c ../geom_circum.c ../geom_predicates.c -I.. -I. -o test_circle_fit $(LIBS)
test_polyline_thicken: test_common.h test_polyline_thicken.c ../geom_arclinegraph.c
	gcc $(CFLAGS) -I../.. test_common.c test_polyline_thicken.c ../geom_arclinegraph.c ../geom_la.c -I.. -I. -o test_polyline_thicken $(LIBS)
