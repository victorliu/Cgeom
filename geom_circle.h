/* Returns number of intersections between two circles (c1,r1) and
 * (c2,r2). The intersection points are (p[0],p[1]) and (p[2],p[3]).
 */
int geom_circle_circle_intersect(
	const double c1[2], double r1,
	const double c2[2], double r2,
	double p[4]
);

/* Returns number of intersections between a circle (c,r) and
 * a line (a,b). The intersection points are the parameters t:
 *   p = a + t*(b-a)
 */
int geom_circle_line_intersect(
	const double c[2], double r,
	const double a[2], const double b[2],
	double t[2]
);
