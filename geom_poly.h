#ifndef GEOM_POLY_H_INCLUDED
#define GEOM_POLY_H_INCLUDED

float  geom_polygon_area2f(unsigned int nv, const float  *v);
double geom_polygon_area2d(unsigned int nv, const double *v);

// Returns 1 if inside, 0 otherwise
int geom_polygon_inside2f(unsigned int nv, const float  *v, const float  p[2]);
int geom_polygon_inside2d(unsigned int nv, const double *v, const double p[2]);

// Returns a normal vector for a point on (near) the polygon's boundary
// For points far from the boundary, no guarantee is made about how reasonable
// the result is.
void geom_polygon_normal2f(unsigned int nv, const float  *v, const float  p[2], float  n[2]);
void geom_polygon_normal2d(unsigned int nv, const double *v, const double p[2], double n[2]);

// A convex region defined as the intersection of n halfspaces.
// Each halfspace is defined by:
//   p[3*i+0]*r[0] + p[3*i+1]*r[1] <= p[3*i+2]
// Therefore, treating p as a column major matrix with three rows
// and np columns, the convex region is defined by:
//   |p^T . r|_inf <= 0
int geom_convex_inside2f(unsigned int np, const float  *p, const float  r[2]);
int geom_convex_inside2d(unsigned int np, const double *p, const double r[2]);

void geom_convex_normal2f(unsigned int np, const float  *p, const float  r[2], float  n[2]);
void geom_convex_normal2d(unsigned int np, const double *p, const double r[2], double n[2]);

// All the above considerations are analogous in 3D.
int geom_convex_inside3f(unsigned int np, const float  *p, const float  r[3]);
int geom_convex_inside3d(unsigned int np, const double *p, const double r[3]);

void geom_convex_normal3f(unsigned int np, const float  *p, const float  r[3], float  n[3]);
void geom_convex_normal3d(unsigned int np, const double *p, const double r[3], double n[3]);

#endif // GEOM_POLY_H_INCLUDED
