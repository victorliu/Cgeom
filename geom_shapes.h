#ifndef GEOM_SHAPES_H_INCLUDED
#define GEOM_SHAPES_H_INCLUDED

// Shape representations
// =====================
//  Each shape is defined relative to its own local coordinate frame.
// This is prefered for numerical robustness and ease of translation.

// For the purposes of mutual intersection tests, there are only the following types:
//   convex, ellipsoid, frustum
typedef enum{
	SHAPE3_SPHERE,
	SHAPE3_BLOCK, // axes not necessarily perpendicular
	SHAPE3_FRUSTUM, // encompasses cylinder and cone
	SHAPE3_ELLIPSOID, // axes not necessarily perpendicular
	SHAPE3_TETRAHEDRON,
	SHAPE3_CONVEX // convex hull of a set of points
} geom_shape3d_type;

typedef struct{
	double r;
} geom_shape3d_sphere;

typedef struct{
	double Q[9]; // orthogonal matrix
	double len;
	double r_base;
	double r_tip_base;
	// The last column of Q is the direction of the axis of the frustum.
	// The axis starts at the local origin and ends at the tip cap.
	// Let the matrix T be defined as
	//   T = Q.diag(r_base,r_base,len)
	// Then
	//   inv(T) = inv(diag(len,r_base,r_base)) . Q^T
	// Let r = inv(T).{x,y,z}
	// The frustum is defined by the set
	//   max(
	//     2*norm_inf(r.z - 0.5),
	//     norm_2(r.xy)^2 / (1 + (r_tip_base-1)*r.z)
	//   ) <= 1
} geom_shape3d_frustum;

typedef struct{
	double A[9]; // semi axes in columns
	double B[9];
	// The block is defined by the set
	//   norm(B . {x,y,z})_inf = 1
	//   B = inv(A)
} geom_shape3d_block;

typedef struct{
	double A[9]; // semi axes in columns
	double B[9];
	// The ellipsoid is defined by the set
	//   (B . {x,y,z})^2 = 1
	//   B = inv(A)
} geom_shape3d_ellipsoid;

typedef struct{
	double p[9]; // other points
} geom_shape3d_tetrahedron;

typedef struct{
	unsigned int np;
	double *p; // plane points (triplets)
	double *n; // plane normals (triplets)
} geom_shape3d_convex;

typedef struct{
	geom_shape3d_type type;
	double o[3]; // origin point
	union{
		geom_shape3d_sphere sphere;
		geom_shape3d_frustum frustum;
		geom_shape3d_block block;
		geom_shape3d_ellipsoid ellipsoid;
		geom_shape3d_tetrahedron tetrahedron;
		geom_shape3d_convex convex;
	} s;
} geom_shape3d;



// For the purposes of mutual intersection tests, there are only 2 types:
//   ellipse, arcspline
typedef enum{
	SHAPE2_CIRCLE,
	SHAPE2_ELLIPSE,
	SHAPE2_QUAD,
	SHAPE2_POLYGON,
	SHAPE2_ARCSPLINE
} geom_shape2d_type;

typedef struct{
	double r;
} geom_shape2d_circle;

// An ellipse may not be an actual ellipse if the axes are not orthogonal
// In this case, the ellipse is the locus of points p such that
// [e[0]/|e[0]|^2 . (p-o)]^2 + [e[1]/|e[1]|^2 . (p-o)]^2 = 1
typedef struct{
	double A[4]; // semi axes of quad in ellipse
	double B[4];
	// The block is defined by the set
	//   norm(B . {x,y,z})_inf = 1
	//   B = inv(A)
} geom_shape2d_ellipse;

typedef struct{
	double A[4]; // semi axes of quad in columns
	double B[4];
	// The block is defined by the set
	//   norm(B . {x,y,z})_inf = 1
	// B = inv(A)
} geom_shape2d_quad;

typedef struct{
	unsigned int n;
	double *p; // points relative to local origin (pairs)
} geom_shape2d_polygon;

typedef struct{
	unsigned int n;
	double *p; // points relative to local origin (pairs)
	double *h; // heights h[i] between p[i] and p[i+1]
} geom_shape2d_arcspline;

typedef struct{
	geom_shape2d_type type;
	double o[2]; // origin point
	union{
		geom_shape2d_circle circle;
		geom_shape2d_ellipse ellipse;
		geom_shape2d_quad quad;
		geom_shape2d_polygon polygon;
		geom_shape2d_arcspline arcspline;
	} s;
} geom_shape2d;

typedef struct{
	double c[3]; // center
	double h[3]; // extents
} geom_aabb3d;

typedef struct{
	double c[2];
	double h[2];
} geom_aabb2d;

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
int geom_shape2d_contains(const geom_shape2d *s, const double p[2]); // arcspline not implemented

// Compute the bounding box of the shape.
int geom_shape3d_get_aabb(const geom_shape3d *s, geom_aabb3d *b);
int geom_shape2d_get_aabb(const geom_shape2d *s, geom_aabb2d *b); // arcspline not implemented
/*
// Determines if two shapes intersect. Returns 0 if no, 1 if yes.
int geom_shape3_intersects(const geom_shape3 *s, const geom_shape3 *t); // not implemented
int geom_shape2_intersects(const geom_shape2 *s, const geom_shape2 *t); // not implemented

// Computes the overlapping volume/area between a shape and
// the given box.
int geom_shape3_aabb_overlap(const geom_shape3 *s, const geom_aabb3 *b); // not implemented
int geom_shape2_aabb_overlap(const geom_shape2 *s, const geom_aabb2 *b); // not implemented

// Computes the overlapping volume/area between a shape and
// the given simplex. The orientation of the simplex determines
// the sign of the returned volume.
double geom_shape3_simplex_overlap(const geom_shape3 *s, const double t[3*4]); // not implemented
double geom_shape2_simplex_overlap(const geom_shape2 *s, const double t[2*3]); // not implemented

// Returns an approximate outward normal vector to the shape at the point
// given by p. p should be "near" the surface of the shape, although
// any p should produce some n.
int geom_shape3_normal(const geom_shape3 *s, const double p[3], double n[3]); // not implemented
int geom_shape2_normal(const geom_shape2 *s, const double p[2], double n[2]); // not implemented

// Determines if a shape intersects a given line segment defined by the point
// and vector. Returns the number of intersections (up to 2) in t. The values
// in t are the offsets along v, and are always in the range [0,1].
int geom_shape3_segment_intersect(const geom_shape3 *s,
	const double a[3], const double v[3], double t[2]); // not implemented
int geom_shape2_segment_intersect(const geom_shape2 *s,
	const double a[2], const double v[2], double t[2]); // not implemented
*/
// Returns the fourier transform (real and imag part in ft) of the shape
// at the k-point 2*pi*f. There is no normalization factor to the Fourier
// integral.
int geom_shape3d_fourier_transform(const geom_shape3d *s, const double f[3], double ft[2]);
int geom_shape2d_fourier_transform(const geom_shape2d *s, const double f[2], double ft[2]);

#include <stdio.h>

// Output a 3D shape description to a POVRay block. The content string is output
// after the shape description to allow setting material properties.
int geom_shape3d_output_POVRay(const shape3d *s, FILE *fp, const char *content);
int geom_aabb3d_output_POVRay(const aabb3d *b, FILE *fp, const char *content);

#define SHAPE_OUTPUT_POSTSCRIPT_FILL        1 // apply a fill to the figure
#define SHAPE_OUTPUT_POSTSCRIPT_NOSTROKE    2 // do not stroke the figure
#define SHAPE_OUTPUT_POSTSCRIPT_NOTRANSLATE 4 // draw using only local coords

// Output a 2D shape description to PostScript commands. The options
// are defined above.
int geom_shape2d_output_postscript(const geom_shape2d *s, FILE *fp, int opts);
int geom_aabb2d_output_postscript(const geom_aabb2d *b, FILE *fp, int opts);

#endif // GEOM_SHAPES_H_INCLUDED
