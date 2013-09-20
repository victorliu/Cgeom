/* A circular arc is parameterized by the start and endpoints */
/* and the bulge factor g. Positive bulges make the arc bulge */
/* to the right when looking at the endpoint from the start.  */

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


