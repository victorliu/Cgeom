#include <stdlib.h>
#include <math.h>
#include <Cgeom/geom_predicates.h>

/*****************************************************************************/
/*                                                                           */
/*  tetcircumcenter()   Find the circumcenter of a tetrahedron.              */
/*                                                                           */
/*  The result is returned both in terms of xyz coordinates and xi-eta-zeta  */
/*  coordinates, relative to the tetrahedron's point `a' (that is, `a' is    */
/*  the origin of both coordinate systems).  Hence, the xyz coordinates      */
/*  returned are NOT absolute; one must add the coordinates of `a' to        */
/*  find the absolute coordinates of the circumcircle.  However, this means  */
/*  that the result is frequently more accurate than would be possible if    */
/*  absolute coordinates were returned, due to limited floating-point        */
/*  precision.  In general, the circumradius can be computed much more       */
/*  accurately.                                                              */
/*                                                                           */
/*  The xi-eta-zeta coordinate system is defined in terms of the             */
/*  tetrahedron.  Point `a' is the origin of the coordinate system.          */
/*  The edge `ab' extends one unit along the xi axis.  The edge `ac'         */
/*  extends one unit along the eta axis.  The edge `ad' extends one unit     */
/*  along the zeta axis.  These coordinate values are useful for linear      */
/*  interpolation.                                                           */
/*                                                                           */
/*  If `xi' is NULL on input, the xi-eta-zeta coordinates will not be        */
/*  computed.                                                                */
/*                                                                           */
/*****************************************************************************/

void geom_circum_tet3d(
	const double a[3],
	const double b[3],
	const double c[3],
	const double d[3],
	double circumcenter[3],
	double *xi,
	double *eta,
	double *zeta
){
  double xba, yba, zba, xca, yca, zca, xda, yda, zda;
  double balength, calength, dalength;
  double xcrosscd, ycrosscd, zcrosscd;
  double xcrossdb, ycrossdb, zcrossdb;
  double xcrossbc, ycrossbc, zcrossbc;
  double denominator;
  double xcirca, ycirca, zcirca;

  /* Use coordinates relative to point `a' of the tetrahedron. */
  xba = b[0] - a[0];
  yba = b[1] - a[1];
  zba = b[2] - a[2];
  xca = c[0] - a[0];
  yca = c[1] - a[1];
  zca = c[2] - a[2];
  xda = d[0] - a[0];
  yda = d[1] - a[1];
  zda = d[2] - a[2];
  /* Squares of lengths of the edges incident to `a'. */
  balength = xba * xba + yba * yba + zba * zba;
  calength = xca * xca + yca * yca + zca * zca;
  dalength = xda * xda + yda * yda + zda * zda;
  /* Cross products of these edges. */
  xcrosscd = yca * zda - yda * zca;
  ycrosscd = zca * xda - zda * xca;
  zcrosscd = xca * yda - xda * yca;
  xcrossdb = yda * zba - yba * zda;
  ycrossdb = zda * xba - zba * xda;
  zcrossdb = xda * yba - xba * yda;
  xcrossbc = yba * zca - yca * zba;
  ycrossbc = zba * xca - zca * xba;
  zcrossbc = xba * yca - xca * yba;

  /* Calculate the denominator of the formulae. */
#ifdef EXACT
  /* Use orient3d() from http://www.cs.cmu.edu/~quake/robust.html     */
  /*   to ensure a correctly signed (and reasonably accurate) result, */
  /*   avoiding any possibility of division by zero.                  */
  denominator = 0.5 / geom_orient3d(b, c, d, a);
#else
  /* Take your chances with floating-point roundoff. */
  denominator = 0.5 / (xba * xcrosscd + yba * ycrosscd + zba * zcrosscd);
#endif

  /* Calculate offset (from `a') of circumcenter. */
  xcirca = (balength * xcrosscd + calength * xcrossdb + dalength * xcrossbc) *
           denominator;
  ycirca = (balength * ycrosscd + calength * ycrossdb + dalength * ycrossbc) *
           denominator;
  zcirca = (balength * zcrosscd + calength * zcrossdb + dalength * zcrossbc) *
           denominator;
  circumcenter[0] = a[0]+xcirca;
  circumcenter[1] = a[1]+ycirca;
  circumcenter[2] = a[2]+zcirca;

  if (xi != (double *) NULL) {
    /* To interpolate a linear function at the circumcenter, define a    */
    /*   coordinate system with a xi-axis directed from `a' to `b',      */
    /*   an eta-axis directed from `a' to `c', and a zeta-axis directed  */
    /*   from `a' to `d'.  The values for xi, eta, and zeta are computed */
    /*   by Cramer's Rule for solving systems of linear equations.       */
    *xi = (xcirca * xcrosscd + ycirca * ycrosscd + zcirca * zcrosscd) *
          (2.0 * denominator);
    *eta = (xcirca * xcrossdb + ycirca * ycrossdb + zcirca * zcrossdb) *
           (2.0 * denominator);
    *zeta = (xcirca * xcrossbc + ycirca * ycrossbc + zcirca * zcrossbc) *
            (2.0 * denominator);
  }
}

/*****************************************************************************/
/*                                                                           */
/*  tricircumcenter()   Find the circumcenter of a triangle.                 */
/*                                                                           */
/*  The result is returned both in terms of x-y coordinates and xi-eta       */
/*  coordinates, relative to the triangle's point `a' (that is, `a' is       */
/*  the origin of both coordinate systems).  Hence, the x-y coordinates      */
/*  returned are NOT absolute; one must add the coordinates of `a' to        */
/*  find the absolute coordinates of the circumcircle.  However, this means  */
/*  that the result is frequently more accurate than would be possible if    */
/*  absolute coordinates were returned, due to limited floating-point        */
/*  precision.  In general, the circumradius can be computed much more       */
/*  accurately.                                                              */
/*                                                                           */
/*  The xi-eta coordinate system is defined in terms of the triangle.        */
/*  Point `a' is the origin of the coordinate system.  The edge `ab' extends */
/*  one unit along the xi axis.  The edge `ac' extends one unit along the    */
/*  eta axis.  These coordinate values are useful for linear interpolation.  */
/*                                                                           */
/*  If `xi' is NULL on input, the xi-eta coordinates will not be computed.   */
/*                                                                           */
/*****************************************************************************/

void geom_circum_tri2d(
	const double a[2],
	const double b[2],
	const double c[2],
	double circumcenter[2],
	double *xi,
	double *eta
){
  double xba, yba, xca, yca;
  double balength, calength;
  double denominator;
  double xcirca, ycirca;
  
  /* Use coordinates relative to point `a' of the triangle. */
  xba = b[0] - a[0];
  yba = b[1] - a[1];
  xca = c[0] - a[0];
  yca = c[1] - a[1];
  /* Squares of lengths of the edges incident to `a'. */
  balength = xba * xba + yba * yba;
  calength = xca * xca + yca * yca;

  /* Calculate the denominator of the formulae. */
#ifdef EXACT
  /* Use orient2d() from http://www.cs.cmu.edu/~quake/robust.html     */
  /*   to ensure a correctly signed (and reasonably accurate) result, */
  /*   avoiding any possibility of division by zero.                  */
  denominator = 0.5 / geom_orient2d(b, c, a);
#else
  /* Take your chances with floating-point roundoff. */
  denominator = 0.5 / (xba * yca - yba * xca);
#endif

  /* Calculate offset (from `a') of circumcenter. */
  xcirca = (yca * balength - yba * calength) * denominator;  
  ycirca = (xba * calength - xca * balength) * denominator;  
  circumcenter[0] = a[0]+xcirca;
  circumcenter[1] = a[1]+ycirca;

  if (xi != (double *) NULL) {
    /* To interpolate a linear function at the circumcenter, define a     */
    /*   coordinate system with a xi-axis directed from `a' to `b' and    */
    /*   an eta-axis directed from `a' to `c'.  The values for xi and eta */
    /*   are computed by Cramer's Rule for solving systems of linear      */
    /*   equations.                                                       */
    *xi = (xcirca * yca - ycirca * xca) * (2.0 * denominator);
    *eta = (ycirca * xba - xcirca * yba) * (2.0 * denominator);
  }
}

/*****************************************************************************/
/*                                                                           */
/*  tricircumcenter3d()   Find the circumcenter of a triangle in 3D.         */
/*                                                                           */
/*  The result is returned both in terms of xyz coordinates and xi-eta       */
/*  coordinates, relative to the triangle's point `a' (that is, `a' is       */
/*  the origin of both coordinate systems).  Hence, the xyz coordinates      */
/*  returned are NOT absolute; one must add the coordinates of `a' to        */
/*  find the absolute coordinates of the circumcircle.  However, this means  */
/*  that the result is frequently more accurate than would be possible if    */
/*  absolute coordinates were returned, due to limited floating-point        */
/*  precision.  In general, the circumradius can be computed much more       */
/*  accurately.                                                              */
/*                                                                           */
/*  The xi-eta coordinate system is defined in terms of the triangle.        */
/*  Point `a' is the origin of the coordinate system.  The edge `ab' extends */
/*  one unit along the xi axis.  The edge `ac' extends one unit along the    */
/*  eta axis.  These coordinate values are useful for linear interpolation.  */
/*                                                                           */
/*  If `xi' is NULL on input, the xi-eta coordinates will not be computed.   */
/*                                                                           */
/*****************************************************************************/

void geom_circum_tri3d(
	const double a[3],
	const double b[3],
	const double c[3],
	double circumcenter[3],
	double *xi,
	double *eta
){
  double xba, yba, zba, xca, yca, zca;
  double balength, calength;
  double xcrossbc, ycrossbc, zcrossbc;
  double denominator;
  double xcirca, ycirca, zcirca;

  /* Use coordinates relative to point `a' of the triangle. */
  xba = b[0] - a[0];
  yba = b[1] - a[1];
  zba = b[2] - a[2];
  xca = c[0] - a[0];
  yca = c[1] - a[1];
  zca = c[2] - a[2];
  /* Squares of lengths of the edges incident to `a'. */
  balength = xba * xba + yba * yba + zba * zba;
  calength = xca * xca + yca * yca + zca * zca;
  
  /* Cross product of these edges. */
#ifdef EXACT
  /* Use orient2d() from http://www.cs.cmu.edu/~quake/robust.html     */
  /*   to ensure a correctly signed (and reasonably accurate) result, */
  /*   avoiding any possibility of division by zero.                  */
  xcrossbc = geom_orient2d(&b[1], &c[1], &a[1]);
  {
   double a2[2] = { a[2], a[0] };
   double b2[2] = { b[2], b[0] };
   double c2[2] = { c[2], c[0] };
   ycrossbc = geom_orient2d(b2, c2, a2);
  }
  zcrossbc = geom_orient2d(b, c, a);
#else
  /* Take your chances with floating-point roundoff. */
  xcrossbc = yba * zca - yca * zba;
  ycrossbc = zba * xca - zca * xba;
  zcrossbc = xba * yca - xca * yba;
#endif

  /* Calculate the denominator of the formulae. */
  denominator = 0.5 / (xcrossbc * xcrossbc + ycrossbc * ycrossbc +
                       zcrossbc * zcrossbc);

  /* Calculate offset (from `a') of circumcenter. */
  xcirca = ((balength * yca - calength * yba) * zcrossbc -
            (balength * zca - calength * zba) * ycrossbc) * denominator;
  ycirca = ((balength * zca - calength * zba) * xcrossbc -
            (balength * xca - calength * xba) * zcrossbc) * denominator;
  zcirca = ((balength * xca - calength * xba) * ycrossbc -
            (balength * yca - calength * yba) * xcrossbc) * denominator;
  circumcenter[0] = a[0]+xcirca;
  circumcenter[1] = a[1]+ycirca;
  circumcenter[2] = a[2]+zcirca;

  if (xi != (double *) NULL) {
    /* To interpolate a linear function at the circumcenter, define a     */
    /*   coordinate system with a xi-axis directed from `a' to `b' and    */
    /*   an eta-axis directed from `a' to `c'.  The values for xi and eta */
    /*   are computed by Cramer's Rule for solving systems of linear      */
    /*   equations.                                                       */

    /* There are three ways to do this calculation - using xcrossbc, */
    /*   ycrossbc, or zcrossbc.  Choose whichever has the largest    */
    /*   magnitude, to improve stability and avoid division by zero. */
    if (((xcrossbc >= ycrossbc) ^ (-xcrossbc > ycrossbc)) &&
        ((xcrossbc >= zcrossbc) ^ (-xcrossbc > zcrossbc))) {
      *xi = (ycirca * zca - zcirca * yca) / xcrossbc;
      *eta = (zcirca * yba - ycirca * zba) / xcrossbc;
    } else if ((ycrossbc >= zcrossbc) ^ (-ycrossbc > zcrossbc)) {
      *xi = (zcirca * xca - xcirca * zca) / ycrossbc;
      *eta = (xcirca * zba - zcirca * xba) / ycrossbc;
    } else {
      *xi = (xcirca * yca - ycirca * xca) / zcrossbc;
      *eta = (ycirca * xba - xcirca * yba) / zcrossbc;
    }
  }
}

double geom_circum_fit2d(
	int n,
	const double *a, /* length 2*n of (x,y) pairs */
	double c[2], double *r
){
	/* We should have
	 *   (x - cx)^2 + (y - cy)^2 = r^2
	 * Let cx = mx + vx, cy = my + vy where
	 *   (mx,my) = avg(x,y)
	 * Then, we have
	 *   (x - mx - vx)^2 + (y - my - vy)^2 = r^2
	 *   (x - mx)^2 - 2(x-mx)vx + (y - my)^2 - 2(y-my)vy = r^2 - vx^2 - vy^2
	 * This is linear in (q, vx, vy) where q = r^2 - vx^2 - vy^2:
	 *   [ 1  2(x-mx)  2(y-my) ] [ q;vx;vy ] = (x - mx)^2 + (y - my)^2
	 * which we write as A*u = b
	 *
	 * Let's partition A = [ 1 B ] and u = [ q;v ] where 1 = ones;
	 * Note that with our choice of M, the columns of B have mean zero, so
	 *   1' * B = 0
	 * Therefore, forming the normal equations, we have
	 *   [ n    0  ] [ q ] = [ 1'*b ]
	 *   [ 0  B'*B ] [ v ]   [ B'*b ]
	 * We then only need to solve q = avg(b), and B'*B*v = B'*b
	 */
	int i;
	double m[2] = { 0,0 };
	double BB[3] = { 0,0,0 };
	double b[2] = { 0,0 };
	const double in = 1./(double)n;
	*r = 0;
	for(i = 0; i < n; ++i){
		m[0] += a[2*i+0];
		m[1] += a[2*i+1];
	}
	m[0] *= in; m[1] *= in;
	for(i = 0; i < n; ++i){
		const double t[2] = { a[2*i+0]-m[0], a[2*i+1]-m[1] };
		const double bi = t[0]*t[0] + t[1]*t[1];
		const double Bi0 = 2*t[0];
		const double Bi1 = 2*t[1];
		b[0] += Bi0 * bi;
		b[1] += Bi1 * bi;
		*r += bi;
		BB[0] += Bi0*Bi0;
		BB[1] += Bi0*Bi1;
		BB[2] += Bi1*Bi1;
	}
	*r *= in;
	
	/* At this point BB contains B'*B, *r contains q */
	/* Solve the system
	 *   [ BB[0] BB[1] ] [ vx ] = [ b[0] ]
	 *   [ BB[1] BB[2] ] [ vy ]   [ b[1] ]
	 */
	{
		double idet = 1. / (BB[0]*BB[2] - BB[1]*BB[1]);
		double v[2] = {
			idet * (BB[2] * b[0] - BB[1] * b[1]),
			idet * (BB[0] * b[1] - BB[1] * b[0])
		};
		/* r = sqrt(q + v'*v) */
		*r = sqrt(*r + v[0]*v[0] + v[1]*v[1]);
		c[0] = m[0]+v[0]; c[1] = m[1]+v[1];
		return sqrt(hypot(
			BB[0]*v[0] + BB[1]*v[1] - b[0],
			BB[1]*v[0] + BB[2]*v[1] - b[1]
		)) * in;
	}
}

