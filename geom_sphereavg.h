#ifndef GEOM_SPHEREAVG_H_INCLUDED
#define GEOM_SPHEREAVG_H_INCLUDED

/* Perform weighted averages on a sphere.
 * Given n vectors on the unit sphere in v as successive triplets,
* and n weights in w, computes the weighted average vector in avg.
 */
void geom_sphereavg3d(int n, const double *v, const double *w, double tol, double *avg);

/* Returns sin(a*x)/sin(x) for -1 <= a <= 1
 * This function is typically used to implement the slerp function,
 * and must be robust for small values of x. Note that this function
 * is not appropriate when x is near +/-Pi, since the slerp function
 * becomes ill-defined. When x is near 0, then the slerp function
 * can still be regularized in a reasonable way.
 */
double geom_sin_ratio(double a, double x);

/* Returns asin(rs*sx)/asin(sx) for -1 <= sx <= 1, -1 <= rs <= 1
 * This function is the inverse of the slerp function when
 * given values of the sin-ratio and sin(x), returning a.
 */
double geom_asin_ratio(double rs, double sx);

/* returns c that is the circular interpolation of a and b
 * (not necessarily normalized, but equal length). s is assumed to be in [0,1]
 */
void geom_slerp2d(const double a[2], const double b[2], double s, double c[2]);


double geom_unslerp2d(const double a[2], const double b[2], const double c[2]);

#endif /* GEOM_SPHEREAVG_H_INCLUDED */
