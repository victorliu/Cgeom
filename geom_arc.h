/* A circular arc is parameterized by the start and endpoints */
/* and the bulge factor g. Positive bulges make the arc bulge */
/* to the right when looking at the endpoint from the start.  */

/* Computes g when arc is specified by endpoints a and b and
 * point c within arc.
 */
double geom_arc_g_from_pt(const double *a, const double *b, const double *c);

/* returns the length of the arc */
double geom_arc_length(const double a[2], const double b[2], double g);

/* Computes a tight bounding circle for the arc */
void geom_arc_bound_circle(
	const double a[2], const double b[2], double g,
	double c[2], double *r
);

/* Computes the extremal point on the arc that is most in the */
/* direction of v (not necessarily normalized). */
void geom_arc_extremum(
	const double a[2], const double b[2], double g,
	const double v[2], double p[2]
);

/* Computes a tight bounding rectangle for the arc */
void geom_arc_bound_rect(
	const double a[2], const double b[2], double g,
	double xb[2], double yb[2]
);

/* Computes the point p and optionally tangent vector t */
/* of the arc at arc length parameter s in [0,1]. Note  */
/* that t is not normalized. */
void geom_arc_param(
	const double a[2], const double b[2], double g,
	double s, double p[2], double t[2]
);
/* Given an arc and a point p that is approximately on the arc,
 * returns the arc-length parameter in [0,1] of p. This function
 * is only guaranteed to give reasonable answers for outputs of
 * geom_arc_param.
 */
double geom_arc_unparam(
	const double a[2], const double b[2], double g,
	const double p[2]
);

/* Returns the circle parameterization of the arc. */
void geom_arc_circle(
	const double a[2], const double b[2], double g,
	double c[2], double *r, double theta[2]
);

/* Computes a cubic Bezier approximation to the arc. */
/* Returns the two inner control points. */
void geom_arc_bezier(
	const double a[2], const double b[2], double g,
	double ap[2], double bp[2]
);

/* Given an arc, splits it into x- and y-monotone segments.
 * Note that there can be up to 4 such segments.
 * The number of segments is returned, and pg contains
 * tuples of {x,y,g} where x,y are the coordinates
 * of the endpoints of the segments (without duplicates)
 * and g is the bulge factor between the i-th and i+1-th point.
 * Therefore, pg can possibly be 3*6-1=17 in length.
 */
int geom_arc_split_monotone(
	const double a[2], const double b[2], double g,
	double *pg
);

/* Given two arcs, computes their intersections points, if any.
 * There may be up to two intersection points. The number of
 * intersections is returned.
 */
int geom_arc_intersect(
	const double a1[2], const double b1[2], double g1,
	const double a2[2], const double b2[2], double g2,
	double p[4]
);

/* Compute the offset curve to the arc.
 * Note that g stays the same.
 * Otherwise, the arc is extended by d on both ends.
 */
void geom_arc_offset(
	const double a[2], const double b[2], double g,
	double d, double ao[2], double bo[2]
);

void geom_arc_extend(
	const double a[2], const double b[2], double g,
	const double d[2], double ao[2], double bo[2], double *go
);
