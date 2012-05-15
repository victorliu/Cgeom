#ifndef GEOM_SHAPESET_H_INCLUDED
#define GEOM_SHAPESET_H_INCLUDED

#include <geom_shapes.h>

// A shapeset is a collection of shapes. This object only stores pointers
// to shapes without managing memory.

typedef struct geom_shapeset2d_struct* geom_shapeset2d;
typedef struct geom_shapeset3d_struct* geom_shapeset3d;

geom_shapeset2d geom_shapeset2d_new(const double lattice[4]);
geom_shapeset3d geom_shapeset3d_new(const double lattice[6]);

void geom_shapeset2d_destroy(geom_shapeset2d ss);
void geom_shapeset3d_destroy(geom_shapeset3d ss);

// returns index of shape
int geom_shapeset2d_add(geom_shapeset2d ss, geom_shape2d *s);
int geom_shapeset3d_add(geom_shapeset3d ss, geom_shape3d *s);

void geom_shapeset2d_finalize(geom_shapeset2d ss);
void geom_shapeset3d_finalize(geom_shapeset3d ss);

unsigned int geom_shapeset2d_size(geom_shapeset2d ss);
unsigned int geom_shapeset3d_size(geom_shapeset3d ss);

geom_shape2d* geom_shapeset2d_index(geom_shapeset2d ss, int index);
geom_shape3d* geom_shapeset3d_index(geom_shapeset3d ss, int index);

// Returns item of largest index
int geom_shapeset2d_query_pt(geom_shapeset2d ss, const double p[2]);
int geom_shapeset3d_query_pt(geom_shapeset3d ss, const double p[3]);

#endif // GEOM_SHAPESET_H_INCLUDED
