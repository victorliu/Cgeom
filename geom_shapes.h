#ifndef GEOM_SHAPES_H_INCLUDED
#define GEOM_SHAPES_H_INCLUDED

// Shape representations
// =====================
//  Each shape is defined relative to its own local coordinate frame.
// This is prefered for numerical robustness and ease of translation.
//
// All matrices are column major:
//   [ A[0]  A[2] ]
//   [ A[1]  A[3] ]
//
// The shape structs themselves should be allocated dynamically due to
// the variable sized arrays at the end for certain types of shapes.

typedef enum{
	GEOM_SHAPE2D_ELLIPSE,
	GEOM_SHAPE2D_POLYGON
} geom_shape2d_type;

// An ellipse may not be an actual ellipse if the axes are not orthogonal
// In this case, the ellipse is the locus of points p such that
// [e[0]/|e[0]|^2 . (p-o)]^2 + [e[1]/|e[1]|^2 . (p-o)]^2 = 1
typedef struct{
	double A[4]; // semi axes of ellipse in columns
	double B[4]; // inverse of A
	// The block is defined by the set
	//   norm(B . {x,y})_2 = 1
	//   B = inv(A)
} geom_shape2d_ellipse;

typedef struct{
	unsigned int nv;
	double v[2]; // vertices relative to local origin (xy pairs)
	// v is a variable sized array of size 2*nv
} geom_shape2d_polygon;

typedef struct{
	geom_shape2d_type type;
	int tag; // user data
	double org[2]; // origin point
	union{
		geom_shape2d_ellipse ellipse;
		geom_shape2d_polygon polygon;
	} s;
} geom_shape2d;


typedef enum{
	GEOM_SHAPE3D_BLOCK, // parallelopiped, axes not necessarily perpendicular
	GEOM_SHAPE3D_POLY, // convex polyhedron; intersection of a set of halfspaces
	GEOM_SHAPE3D_ELLIPSOID, // axes not necessarily perpendicular
	GEOM_SHAPE3D_FRUSTUM, // encompasses cylinder and cone
	GEOM_SHAPE3D_EXTRUSION // a 2D shape extruded along a 3rd axis
} geom_shape3d_type;

typedef struct{
	double A[9]; // semi axes in columns
	double B[9];
	// The block is defined by the set
	//   norm(B . {x,y,z})_inf = 1
	//   B = inv(A)
} geom_shape3d_block;

typedef struct{
	unsigned int np;
	double p[4]; // plane normals (normal + offset)
	// p is a variable sized array of size 4*np
	// Halfspaces defined by:
	//   p[4*i+0] * x + p[4*i+1] * y + p[4*i+2] * z <= d
	// If p is a column-major matrix, then we have the region defined by
	//    norm(p^T . {x,y,z,1})_inf <= 0
} geom_shape3d_poly;

typedef struct{
	double Q[9]; // orthogonal matrix
	double len;
	double r_base;
	double r_tip;
	// The last column of Q is the direction of the axis of the frustum.
	// The axis starts at the local origin and ends at the tip cap.
	// Let the matrix T be defined as
	//   T = Q.diag(r_base,r_base,len)
	// Then
	//   inv(T) = inv(diag(r_base,r_base,len)) . Q^T
	// Let r = inv(T).{x,y,z}
	// The frustum is defined by the set
	//   max(
	//     2*norm_inf(r.z - 0.5),
	//     norm_2(r.xy) / (1 + (r_tip-1)*r.z)
	//   ) <= 1
} geom_shape3d_frustum;

typedef struct{
	double A[9]; // semi axes in columns
	double B[9];
	// The ellipsoid is defined by the set
	//   (B . {x,y,z})^2 = 1
	//   B = inv(A)
} geom_shape3d_ellipsoid;

typedef struct{
	double Q[9]; // orthogonal matrix
	double len;
	geom_shape2d s2;
	// The first two columsn of Q define the "xy-plane"
	// in which the 2D shape s2 lives. The third column is
	// the extrusion axis of length len.
} geom_shape3d_extrusion;

typedef struct{
	geom_shape3d_type type;
	int tag; // user data
	double org[3]; // origin point
	union{
		geom_shape3d_block     block;
		geom_shape3d_poly      poly;
		geom_shape3d_ellipsoid ellipsoid;
		geom_shape3d_frustum   frustum;
		geom_shape3d_extrusion extrusion;
	} s;
} geom_shape3d;


// Axis Aligned Bounding Boxes
// Represented by a center point c and half widths h.
// This representation is preferred so that small boxes at large
// translations are represented to high precision, and so that
// translation does not affect the precision of the size of the box.
typedef struct{
	double c[2], h[2];
} geom_aabb2d;

typedef struct{
	double c[3], h[3];
} geom_aabb3d;


// These functions do not check input arguments for NULL.
// All operations return false if not supported

// Expands the box b1 to include the box b2.
void geom_aabb3d_union(geom_aabb3d *b1, const geom_aabb3d *b2);
void geom_aabb2d_union(geom_aabb2d *b1, const geom_aabb2d *b2);

// Expands the box to include the point.
void geom_aabb3d_union_pt(geom_aabb3d *b1, const double p[3]);
void geom_aabb2d_union_pt(geom_aabb2d *b1, const double p[2]);

// Determines if p lies within the box. Returns 0 if no, 1 if yes.
int geom_aabb3d_contains(const geom_aabb3d *b, const double p[3]);
int geom_aabb2d_contains(const geom_aabb2d *b, const double p[2]);

// Determines if two boxes intersect. Returns 0 if no, 1 if yes.
int geom_aabb3d_intersects(const geom_aabb3d *a, const geom_aabb3d *b);
int geom_aabb2d_intersects(const geom_aabb2d *a, const geom_aabb2d *b);

// Determines if p lies within the shape. Returns 0 if no, 1 if yes.
int geom_shape3d_contains(const geom_shape3d *s, const double p[3]);
int geom_shape2d_contains(const geom_shape2d *s, const double p[2]);

// Compute the bounding box of the shape.
int geom_shape3d_get_aabb(const geom_shape3d *s, geom_aabb3d *b);
int geom_shape2d_get_aabb(const geom_shape2d *s, geom_aabb2d *b);

// Returns an approximate outward normal vector to the shape at the point
// given by p. p should be "near" the surface of the shape, although
// any p should produce some n.
int geom_shape3d_normal(const geom_shape3d *s, const double p[3], double n[3]);
int geom_shape2d_normal(const geom_shape2d *s, const double p[2], double n[2]);

/*
// Determines if two shapes intersect. Returns 0 if no, 1 if yes.
int geom_shape3d_intersects(const geom_shape3 *s, const geom_shape3 *t);
int geom_shape2d_intersects(const geom_shape2 *s, const geom_shape2 *t);

// Computes the overlapping volume/area between a shape and
// the given box.
int geom_shape3d_aabb_overlap(const geom_shape3 *s, const geom_aabb3 *b);
int geom_shape2d_aabb_overlap(const geom_shape2 *s, const geom_aabb2 *b);

// Computes the overlapping volume/area between a shape and
// the given simplex. The orientation of the simplex determines
// the sign of the returned volume.
double geom_shape3d_simplex_overlap(const geom_shape3 *s, const double t[12]);
double geom_shape2d_simplex_overlap(const geom_shape2 *s, const double t[6]);

// Determines if a shape intersects a given line segment defined by the point
// and vector. Returns the number of intersections (up to 2) in t. The values
// in t are the offsets along v, and are always in the range [0,1].
int geom_shape3d_segment_intersect(const geom_shape3 *s,
	const double a[3], const double v[3], double t[2]);
int geom_shape2d_segment_intersect(const geom_shape2 *s,
	const double a[2], const double v[2], double t[2]);

// Returns the fourier transform (real and imag part in ft) of the shape
// at the k-point 2*pi*f. There is no normalization factor to the Fourier
// integral.
int geom_shape3d_fourier_transform(const geom_shape3d *s, const double f[3], double ft[2]);
int geom_shape2d_fourier_transform(const geom_shape2d *s, const double f[2], double ft[2]);
*/
#include <stdio.h>

// Output a 3D shape description to a POVRay block. The content string is output
// after the shape description to allow setting material properties.
int geom_shape3d_output_POVRay(const geom_shape3d *s, FILE *fp, const char *content);
int geom_aabb3d_output_POVRay(const geom_aabb3d *b, FILE *fp, const char *content);

#define GEOM_SHAPE2D_OUTPUT_POSTSCRIPT_FILL        0x1 // apply a fill to the figure
#define GEOM_SHAPE2D_OUTPUT_POSTSCRIPT_NOSTROKE    0x2 // do not stroke the figure
#define GEOM_SHAPE2D_OUTPUT_POSTSCRIPT_NOTRANSLATE 0x4 // draw using only local coords

// Output a 2D shape description to PostScript commands. The options
// are defined above.
int geom_shape2d_output_postscript(const geom_shape2d *s, FILE *fp, int opts);
int geom_aabb2d_output_postscript(const geom_aabb2d *b, FILE *fp, int opts);

#endif // GEOM_SHAPES_H_INCLUDED
