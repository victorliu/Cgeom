#include <Cgeom/geom_la.h>

int geom_circle_circle_intersect(
	const double c1[2], double r1,
	const double c2[2], double r2,
	double p[4]
){
	const double d[2] = { c2[0]-c1[0], c2[1]-c1[1] };
	const double dl = geom_norm2d(d);
	const double rsum = r1+r2;
	const double rdiff = fabs(r1-r2);
	if(dl > rsum){ return 0; }
	if(dl == rsum){
		double t = r1/dl;
		p[0] = (1-t)*c1[0] + t*c2[0];
		p[1] = (1-t)*c1[1] + t*c2[1];
		return 0;
	}
	if(dl < rdiff){ return 0; }
	if(dl == rdiff){
		if(r1 > r2){
			double t = r2/dl;
			p[0] = (1+t)c2[0] - t*c1[0];
			p[1] = (1+t)c2[1] - t*c1[1];
			return 1;
		}else{
			double t = r1/dl;
			p[0] = (1+t)c1[0] - t*c2[0];
			p[1] = (1+t)c1[1] - t*c2[1];
			return 1;
		}
	}
	/* At this point, we need to find p s.t.
	 *   |p-c1|^2 == r1^2 and |p-c2|^2 == r2^2
	 * Parameterize p = c1 + t*(c2-c1) + s*rot90(c2-c1)
	 *   t = (r1^2 - r2^2 + d^2) / [2 d^2]
	 *   s = sqrt(r1^2/d^2 - t^2)
	 * At this point we have established that d < r1+r2 and d > |r1-r2|
	 */
	const double q1 = r1/dl;
	const double q2 = r2/dl;
	double t = 0.5*(1+(q1+q2)*(q1-q2));
	double s = (q1-t)*(q1+t);
	if(0 == s){
		p[0] = (1-t)*c[0] + t*c2[0];
		p[1] = (1-t)*c[1] + t*c2[1];
		return 1;
	}else if(s < 0){
		return 0;
	}else{
		s = sqrt(s);
		p[0] = (1-t)*c[0] + t*c2[0] - s*(c2[1]-c1[1]);
		p[1] = (1-t)*c[1] + t*c2[1] + s*(c2[0]-c1[0]);
		p[2] = (1-t)*c[0] + t*c2[0] + s*(c2[1]-c1[1]);
		p[3] = (1-t)*c[1] + t*c2[1] - s*(c2[0]-c1[0]);
		return 2;
	}
}

int geom_circle_line_intersect(
	const double c[2], double r,
	const double a[2], const double b[2],
	double t[2]
){
	/* Let the intersection points be parameterized by
	 *   p = a + t*(b-a)
	 * Then we must have
	 *   [ (a-c) + t*(b-a) ]^2 = r^2
	 *   ca^2 - r^2 + 2*t*ca*ab + t^2*ab^2 = 0
	 */
	const double ab[2] = { b[0]-a[0], b[1]-a[1] };
	const double ca[2] = { a[0]-c[0], a[1]-c[1] };
	const double ab2 = ab[0]*ab[0] + ab[1]*ab[1];
	const double p = (ca[0]*ab[0] + ca[1]*ab[1]) / ab2;
	const double q = (ca[0]*ca[0] + ca[1]*ca[1] - r*r) / ab2;
	/* Solve t^2 + 2*t*p + q = 0 */
	double disc = p*p-q;
	if(disc == 0){
		t[0] = -p;
		return 1;
	}else if(disc < 0){
		return 0;
	}else{
		disc = sqrt(disc);
		t[0] = -p + disc;
		t[1] = -p - disc;
		return 2;
	}
}

