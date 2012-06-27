#ifndef GEOM_SHAPESET_H_INCLUDED
#define GEOM_SHAPESET_H_INCLUDED

#include <Cgeom/geom_shapes.h>

// A shapeset is a collection of shapes. This object only stores pointers
// to shapes without managing memory.

typedef struct geom_shapeset2d_struct* geom_shapeset2d;
typedef struct geom_shapeset3d_struct* geom_shapeset3d;

#define GEOM_SHAPESET2D_FLAG_UNBOUNDED 1
#define GEOM_SHAPESET3D_FLAG_UNBOUNDED 1

geom_shapeset2d geom_shapeset2d_new();
geom_shapeset3d geom_shapeset3d_new();

void geom_shapeset2d_destroy(geom_shapeset2d ss);
void geom_shapeset3d_destroy(geom_shapeset3d ss);

int geom_shapeset2d_set_lattice(geom_shapeset2d ss, const double lattice[4]);
int geom_shapeset3d_set_lattice(geom_shapeset3d ss, const double lattice[9]);

// returns index of shape
int geom_shapeset2d_add(geom_shapeset2d ss, geom_shape2d *s);
int geom_shapeset3d_add(geom_shapeset3d ss, geom_shape3d *s);

void geom_shapeset2d_finalize(geom_shapeset2d ss);
void geom_shapeset3d_finalize(geom_shapeset3d ss);

unsigned int geom_shapeset2d_size(geom_shapeset2d ss);
unsigned int geom_shapeset3d_size(geom_shapeset3d ss);

geom_shape2d* geom_shapeset2d_index(geom_shapeset2d ss, int index);
geom_shape3d* geom_shapeset3d_index(geom_shapeset3d ss, int index);

int geom_shapeset2d_index_aabb(geom_shapeset2d ss, int index, geom_aabb2d *box);
int geom_shapeset3d_index_aabb(geom_shapeset3d ss, int index, geom_aabb3d *box);

// Returns item of largest index
int geom_shapeset2d_query_pt(geom_shapeset2d ss, const double p[2]);
int geom_shapeset3d_query_pt(geom_shapeset3d ss, const double p[3]);

int geom_shapeset2d_foreach(
	geom_shapeset2d ss,
	int (*func)(geom_shape2d *s, const geom_aabb2d *box, unsigned int flags, void *data),
	void *data
);
int geom_shapeset3d_foreach(
	geom_shapeset3d ss,
	int (*func)(geom_shape3d *s, const geom_aabb3d *box, unsigned int flags, void *data),
	void *data
);

#endif // GEOM_SHAPESET_H_INCLUDED
