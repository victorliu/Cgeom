#ifndef GEOM_BVH_H_INCLUDED
#define GEOM_BVH_H_INCLUDED

// Bounding volume hierarchy structures
//
// The underlying implementation of the BVH is an R-tree. We assume
// that the structures are _static_; that is, they are created and
// never modified. This assumption allows us to generate more
// efficient trees by bulk-loading the trees up front.
//  The bulk loading method we use is Sort-Tile-Recursive (STR) due
// to its simplicity and reasonable effectiveness.

typedef struct geom_bvh2d_struct* geom_bvh2d;
typedef struct geom_bvh3d_struct* geom_bvh3d;

// Construct a new BVH from an iterator function.
// The function should fill in the box info (c,h), and an optional tag
// The data parameter is passed to the iterator.
geom_bvh2d geom_bvh2d_new(unsigned int n, int (*shape_iterator)(double c[2], double h[2], int *tag, void *data), void *data);
geom_bvh3d geom_bvh3d_new(unsigned int n, int (*shape_iterator)(double c[2], double h[3], int *tag, void *data), void *data);

void geom_bvh2d_destroy(geom_bvh2d bvh);
void geom_bvh3d_destroy(geom_bvh3d bvh);

// p is the query point
// In 2D, p[2] is {x,y}. In 3D, p[3] is {x,y,z}
// The query function is passed leaf boxes which contain the point p, along with the tag.
// The function should return nonzero to continue the query.
// A zero return value terminates the query.
int geom_bvh2d_query_pt(
	geom_bvh2d bvh,
	const double p[2],
	int (*query_func)(int tag, const double c[2], const double h[2], void *data),
	void *data
);
int geom_bvh3d_query_pt(
	geom_bvh3d bvh,
	const double p[3],
	int (*query_func)(int tag, const double c[3], const double h[3], void *data),
	void *data
);

// Same as query_pt, but returns all leaf boxes which intersect the given box (c,h).
int geom_bvh2d_query_box(
	geom_bvh2d bvh,
	const double c[2], const double h[2],
	int (*query_func)(int tag, const double c[2], const double h[2], void *data),
	void *data
);
int geom_bvh3d_query_box(
	geom_bvh3d bvh,
	const double c[3], const double h[3],
	int (*query_func)(int tag, const double c[3], const double h[3], void *data),
	void *data
);

// Performs a full traversal of all the boxes in the BVH.
// The flag internal is set to 0 for leaf nodes, and 1 for non-leaf nodes.
int geom_bvh2d_traverse(
	geom_bvh2d bvh,
	int (*func)(int tag, const double c[2], const double h[2], int leaf, void *data),
	void *data
);
int geom_bvh3d_traverse(
	geom_bvh3d bvh,
	int (*func)(int tag, const double c[3], const double h[3], int leaf, void *data),
	void *data
);

#endif // GEOM_BVH_H_INCLUDED
