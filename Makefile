CC = gcc
CFLAGS = -Wall -I.. -O0 -ggdb

OBJS = \
	geom_la.o \
	geom_poly.o \
	geom_predicates.o \
	geom_circum.o \
	geom_shapes.o \
	geom_bvh.o \
	geom_shapeset.o \
	geom_sphereavg.o \
	geom_arc.o

all: libgeom.a
libgeom.a: $(OBJS)
	ar rcs libgeom.a $(OBJS)

geom_la.o: geom_la.c geom_la.h
	$(CC) -c $(CFLAGS) geom_la.c -o geom_la.o
geom_poly.o: geom_poly.c geom_poly.h geom_la.h
	$(CC) -c $(CFLAGS) geom_poly.c -o geom_poly.o
geom_predicates.o: geom_predicates.c geom_predicates.h
	$(CC) -c $(CFLAGS) geom_predicates.c -o geom_predicates.o
geom_circum.o: geom_circum.c geom_circum.h
	$(CC) -c $(CFLAGS) geom_circum.c -o geom_circum.o
geom_shapes.o: geom_shapes.c geom_shapes.h geom_poly.h geom_la.h geom_predicates.h
	$(CC) -c $(CFLAGS) geom_shapes.c -o geom_shapes.o
geom_bvh.o: geom_bvh.c geom_bvh.h
	$(CC) -c $(CFLAGS) geom_bvh.c -o geom_bvh.o
geom_shapeset.o: geom_shapeset.c geom_bvh.h geom_shapes.h
	$(CC) -c $(CFLAGS) geom_shapeset.c -o geom_shapeset.o
geom_sphereavg.o: geom_la.h geom_sphereavg.h
	$(CC) -c $(CFLAGS) geom_sphereavg.c -o geom_sphereavg.o
geom_arc.o: geom_arc.c geom_arc.h
	$(CC) -c $(CFLAGS) geom_arc.c -o geom_arc.o

clean:
	rm -f *.o libgeom.a
