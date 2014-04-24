#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <Cgeom/geom_sphereavg.h>

#define G_QUARTER 0.4142135623730950488

/* returns 2 [1-sqrt(1-x^2)] / x^2 */
static double q2ratio(double x){
	if(0 == x){ return 1; }
	else if(x < 0){ x = -x; }
	{
		const double x2 = x*x;
		if(x < 3e-3){ return (0.125*x2 + 0.25)*x2 + 1.; }
		else{
			return 2.*(1. - sqrt((1.+x)*(1.-x)))/x2;
		}
	}
}

/* If neg == 0
 *   Solves a*x^2 + 2*b*x + a = 0 for x
 * Else
 *   Solves a*x^2 + 2*b*x - a = 0 for x
 * Always returns the root of smaller magnitude; the other root
 * is simply +/- 1/x (plus sign for neg == 0)
 * If the result is complex, the real part is returned.
 *
 */
static double qsolve(int neg, double a, double b){
	if(0 == a){
		return 0;
	}
	if(0 == b){
		return 1;
	}
	double r;
	if(fabs(b) > fabs(a)){
		r = a/b;
		if(neg){
			return r / (1. + sqrt(1. + r*r));
		}else{
			return -r / (1. + sqrt(1. - r*r));
		}
	}else{
		r = b/a;
		if(!neg){ /* complex root */
			return -r;
		}else{
			if(r > 0){
				return 1. / (r + sqrt(1.+r*r));
			}else{
				return 1. / (r - sqrt(1.+r*r));
			}
		}
	}
}


/* returns asin(x)/x */
static double asinc(double x){
	if(0 == x){ return 1; }
	else if(x < 0){ x = -x; }
	const double x2 = x*x;
	if(x < 1e-4){
		return (1./6.)*x2 + 1.;
	}else if(x < 3e-3){
		return ((3./40.)*x2 + (1./6.))*x2 + 1.;
	}else{
		return asin(x)/x;
	}
}

/* Given g = tan(x), computes r[0] = tan(f*x) and r[1] = tan((1-f)*x)
 * for f in [0,1]
 */
static void tanfrac(double g, double f, double *g0, double *g1){
	if(0 == f){
		*g0 = 0; *g1 = g;
	}else if(1 == f){
		*g0 = g; *g1 = 0;
	}else{
		/* Let a = r[0], b = r[1]. We have the relation that
		 *   a = (g-b)/(1+g*b)
		 * or
		 *   a + b + g a b - g = 0
		 */
		double a, b;
		int flip = 0;
		if(f > 0.5){
			flip = 1;
			f = 1.-f;
		}
		/* Solve a = tan(f*x) */
		a = tan(f*atan(g));
		b = (g-a) / (1.+g*a);
		if(flip){
			*g0 = b; *g1 = a;
		}else{
			*g0 = a; *g1 = b;
		}
	}
}

double geom_arc_g_from_pt(const double *a, const double *b, const double *c){
	const double mc[2] = {
		c[0] - (0.5*a[0]+0.5*b[0]),
		c[1] - (0.5*a[1]+0.5*b[1])
	};
	const double am[2] = {
		0.5*(b[0]-a[0]),
		0.5*(b[1]-a[1])
	};
	const double iam2 = 1./(am[0]*am[0]+am[1]*am[1]);
	const double alpha = (mc[0]*am[0] + mc[1]*am[1]) * iam2;
	const double beta  = (mc[0]*am[1] - mc[1]*am[0]) * iam2;
	const double g = qsolve(1, 2.*beta, 1. - (alpha*alpha + beta*beta));
	if(beta >= 0){
		if(g >= 0){ return g; }
		else{ return -1./g; }
	}else{
		if(g <= 0){ return g; }
		else{ return -1./g; }
	}
}

double geom_arc_length(const double a[2], const double b[2], double g){
	const double t = 0.5 * hypot(b[0]-a[0], b[1]-a[1]);
	if(0 == g){
		return 2*t;
	}else if(fabs(g) < 0.125){
		const double alpha = 2*g / (1+g*g);
		return t*asinc(alpha);
	}else{
		const double r = t * ((1+g*g) / (2*g));
		const double theta = atan2(2*g, (1+g)*(1-g));
		return fabs(r*theta);
	}
}

void geom_arc_param(
	const double a[2], const double b[2], double g,
	double s, double p[2], double t[2]
){
	if(0 == g){
		p[0] = (1.-s) * a[0] + s * b[0];
		p[1] = (1.-s) * a[1] + s * b[1];
		if(NULL != t){
			t[0] = b[0] - a[0];
			t[1] = b[1] - a[1];
		}
	}else if(fabs(g) > G_QUARTER){
		/* Use slerp */
		const double m[2] = { /* midpoint of a and b */
			0.5*a[0] + 0.5*b[0],
			0.5*a[1] + 0.5*b[1]
		};
		const double dt2 = (1.+g)*(1.-g)/(4.*g);
		const double t2 = hypot(b[0]-a[0], b[1]-a[1]);
		const double r = 0.25*t2*fabs(1./g + g);
		const double c[2] = {
			m[0] + dt2 * (a[1]-b[1]),
			m[1] + dt2 * (b[0]-a[0])
		};
		const double u[2] = { a[0]-c[0], a[1]-c[1] };
		const double v[2] = { b[0]-c[0], b[1]-c[1] };
		const double qu = atan2(u[1], u[0]);
		const double qv = atan2(v[1], v[0]);
		double q;
		if(g > 0){
			q = qv - qu;
			if(q < 0){ q += 2*M_PI; }
			if(q >= 2*M_PI){ q -= 2*M_PI; }
			q = s*q + qu;
		}else{
			q = qu - qv;
			if(q < 0){ q += 2*M_PI; }
			if(q >= 2*M_PI){ q -= 2*M_PI; }
			q = qu - s*q;
		}
		p[0] = r*cos(q);
		p[1] = r*sin(q);
		/*geom_slerp2d(u, v, s, p); apparently inaccurate */
		if(NULL != t){
			if(g > 0){
				t[0] = -p[1]; t[1] =  p[0];
			}else{
				t[0] =  p[1]; t[1] = -p[0];
			}
		}
		if(0 == s){
			p[0] = a[0]; p[1] = a[1];
		}else if(1 == s){
			p[0] = b[0]; p[1] = b[1];
		}else{
			p[0] += c[0]; p[1] += c[1];
		}
	}else{
		const double mv[2] = { /* vector from midpoint to V */
			0.5 * g * (b[1]-a[1]),
			0.5 * g * (a[0]-b[0])
		};
		const double m[2] = { /* midpoint of a and b */
			0.5*a[0] + 0.5*b[0],
			0.5*a[1] + 0.5*b[1]
		};
		const double theta = atan2(2*g, (1+g)*(1-g));
		const double sr0 = geom_sin_ratio(   s, theta);
		const double sr1 = geom_sin_ratio(1.-s, theta);
		if(0 == s){
			p[0] = a[0]; p[1] = a[1];
		}else if(1 == s){
			p[0] = b[0]; p[1] = b[1];
		}else{
			const double lv = 4./(1.+g*g) * sr0 * sr1;
			const double sr = geom_sin_ratio(2.*s - 1., theta);
			const double la = 0.5 * (1. - sr) - 0.5*lv;
			const double lb = 0.5 * (1. + sr) - 0.5*lv;
			p[0] = la*a[0] + lb*b[0] + lv*(m[0]+mv[0]);
			p[1] = la*a[1] + lb*b[1] + lv*(m[1]+mv[1]);
		}
		if(NULL != t){ /* compute tangent vector */
			/* P' = D1 * mb + D2 * mv
			 * D1 = cos((2*s-1)*theta}
			 * D2 = 2/(1+g*g) * [ cos(s*theta) * sr1 - cos((1-s)*theta) * sr0 ] 
			 */
			const double mb[2] = { b[0] - m[0], b[1] - m[1] };
			const double d1 = cos((2.*s-1.) * theta);
			const double d2 = 2./(1.+g*g) * (cos(s*theta) * sr1 - cos((1.-s)*theta) * sr0);
			t[0] = d1*mb[0] + d2*mv[0];
			t[1] = d1*mb[1] + d2*mv[1];
		}
	}
}

double geom_arc_unparam(
	const double a[2], const double b[2], double g,
	const double p[2]
){
	if(0 == g){
		/* Linear interpolation on larger coordinate */
		const double dx = b[0]-a[0];
		const double dy = b[1]-a[1];
		if(fabs(dx) > fabs(dy)){
			return (p[0]-a[0]) / dx;
		}else{
			return (p[1]-a[1]) / dy;
		}
	}else if(fabs(g) > G_QUARTER){
		/* Compute from circle, un-slerp */
		const double m[2] = { /* midpoint of a and b */
			0.5*a[0] + 0.5*b[0],
			0.5*a[1] + 0.5*b[1]
		};
		const double dt2 = (1.+g)*(1.-g)/(4.*g);
		const double c[2] = {
			m[0] + dt2 * (a[1]-b[1]),
			m[1] + dt2 * (b[0]-a[0])
		};
		const double u[2] = { a[0]-c[0], a[1]-c[1] };
		const double v[2] = { b[0]-c[0], b[1]-c[1] };
		const double w[2] = { p[0]-c[0], p[1]-c[1] };
		const double qu = atan2(u[1], u[0]) / (2*M_PI);
		const double qv = atan2(v[1], v[0]) / (2*M_PI);
		double qw = atan2(w[1], w[0]) / (2*M_PI);
/*printf("qu = %g, qv = %g, qw = %g\n", qu, qv, qw);*/
		double q;
		if(g > 0){
			q = qv - qu;
			qw -= qu;
			if(q < 0){ q += 1; }
			if(q >= 1){ q -= 1; }
			if(qw < 0){ qw += 1; }
			if(qw >= 1){ qw -= 1; }
			return qw/q;
		}else{
			q = qu - qv;
			qw = qu - qw;
/*printf("  qu = %g, qv = %g, qw = %g, q = %g\n", qu, qv, qw, q);*/
			if(q < 0){ q += 1; }
			if(q >= 1){ q -= 1; }
			if(qw < 0){ qw += 1; }
			if(qw >= 1){ qw -= 1; }
			return qw/q;
		}
	}else{
		/* Compute barycentric coordinates of p wrt a,b,v */
		const double m[2] = { /* midpoint of a and b */
			0.5*a[0] + 0.5*b[0],
			0.5*a[1] + 0.5*b[1]
		};
		const double mp[2] = { p[0] - m[0], p[1] - m[1] };
		const double mb[2] = { b[0] - m[0], b[1] - m[1] };
		const double t2 = mb[0]*mb[0] + mb[1]*mb[1];
		const double dot = mb[0]*mp[0] + mb[1]*mp[1];
		const double lAB = dot/t2;
		/* At this point lAB is sin((2s-1)theta)/sin(theta).
		 * and the sign of dot determines which side of the middle of the arc
		 * we are on.
		 */
		const double s2 = geom_asin_ratio(lAB, 2*g/(1+g*g));
		return 0.5+0.5*s2;
	}
}

void geom_arc_bound_circle(
	const double a[2], const double b[2], double g,
	double c[2], double *r
){
	/* set c to midpoint m for now */
	c[0] = 0.5*a[0] + 0.5*b[0];
	c[1] = 0.5*a[1] + 0.5*b[1];
	if(fabs(g) <= 1){
		/* c should just be m */
		*r = 0.5*hypot(b[0]-a[0], b[1]-a[1]);
	}else{
		const double t2 = hypot(b[0]-a[0], b[1]-a[1]);
		const double dt2 = (1.+g)*(1.-g)/(4.*g);
		c[0] += dt2 * (a[1]-b[1]);
		c[1] += dt2 * (b[0]-a[0]);
		/* *r = t * (1+g*g)/(2*g); */
		/* Since g > 1, we can do better: */
		*r = 0.25 * t2 * (1./g + g);
	}
}

void geom_arc_extremum(
	const double a[2], const double b[2], double g,
	const double v[2], double p[2]
){
	const double vlen = hypot(v[0], v[1]);
	double ta[2], tb[2];
	double dot = v[0]*(b[0]-a[0]) + v[1]*(b[1]-a[1]);
	
	int bint = 0;
	if(0 != g){
		geom_arc_param(a, b, g, 0, p, ta);
		geom_arc_param(a, b, g, 1, p, tb);
		double dota = ta[0]*v[0]+ta[1]*v[1];
		double dotb = tb[0]*v[0]+tb[1]*v[1];
		if(fabs(g) < 1){
			bint = (dota > 0 && dotb < 0);
		}else{
			bint = (dota > 0 || dotb < 0);
		}
	}
	if(bint){
		const double t2 = hypot(b[0]-a[0], b[1]-a[1]);
		const double mv[2] = { /* vector from midpoint to V */
			0.5 * g * (b[1]-a[1]),
			0.5 * g * (a[0]-b[0])
		};
		const double m[2] = { /* midpoint of a and b */
			0.5*a[0] + 0.5*b[0],
			0.5*a[1] + 0.5*b[1]
		};
		const double gg1 = 1.+g*g;
		if(fabs(g) < 0.125){
			const double alpha = dot / (vlen * t2);
			double alphap = alpha * gg1/fabs(2.*g);
			if(alphap < -1){ alphap = -1; }
			if(alphap >  1){ alphap =  1; }
			const double betap = q2ratio(alpha) * alphap*alphap / gg1;
			const double la = 0.5*betap - 0.5*alphap;
			const double lb = 0.5*betap + 0.5*alphap;
			const double lv = 1.-betap;
			p[0] = la*a[0] + lb*b[0] + lv*(m[0]+mv[0]);
			p[1] = la*a[1] + lb*b[1] + lv*(m[1]+mv[1]);
		}else{
			const double dt2 = (1.+g)*(1.-g)/(4.*g);
			const double c[2] = {
				m[0] + dt2 * (a[1]-b[1]),
				m[1] + dt2 * (b[0]-a[0])
			};
			const double r_v = fabs(t2 * gg1 / (4.*g * vlen));
			p[0] = c[0] + r_v * v[0];
			p[1] = c[1] + r_v * v[1];
		}
	}else{ // same as segment extremum
		if(dot > 0){
			p[0] = b[0];
			p[1] = b[1];
		}else{
			p[0] = a[0];
			p[1] = a[1];
		}
	}
}

void geom_arc_bound_rect(
	const double a[2], const double b[2], double g,
	double xb[2], double yb[2]
){
	double v[2], p[2];
	if(a[0] < b[0]){
		xb[0] = a[0];
		xb[1] = b[0];
	}else{
		xb[0] = b[0];
		xb[1] = a[0];
	}
	if(a[1] < b[1]){
		yb[0] = a[1];
		yb[1] = b[1];
	}else{
		yb[0] = b[1];
		yb[1] = a[1];
	}
	if(0 == g){ return; }
	
	
	v[0] = -1; v[1] = 0;
	geom_arc_extremum(a, b, g, v, p);
	if(p[0] < xb[0]){ xb[0] = p[0]; }
	
	v[0] = 1; v[1] = 0;
	geom_arc_extremum(a, b, g, v, p);
	if(p[0] > xb[1]){ xb[1] = p[0]; }
	
	v[0] = 0; v[1] = -1;
	geom_arc_extremum(a, b, g, v, p);
	if(p[1] < yb[0]){ yb[0] = p[1]; }
	
	v[0] = 0; v[1] = 1;
	geom_arc_extremum(a, b, g, v, p);
	if(p[1] > yb[1]){ yb[1] = p[1]; }
}

void geom_arc_circle(
	const double a[2], const double b[2], double g,
	double c[2], double *r, double theta[2]
){
	const double t2 = hypot(b[0]-a[0], b[1]-a[1]);
	const double dt2 = (1.+g)*(1.-g)/(4.*g);
	c[0] = (0.5*a[0] + 0.5*b[0]) + dt2 * (a[1]-b[1]);
	c[1] = (0.5*a[1] + 0.5*b[1]) + dt2 * (b[0]-a[0]);
	*r = 0.25 * t2 * fabs(1./g + g);
	theta[0] = atan2(a[1]-c[1], a[0]-c[0]);
	theta[1] = atan2(b[1]-c[1], b[0]-c[0]);
}

void geom_arc_bezier(
	const double a[2], const double b[2], double g,
	double ap[2], double bp[2]
){
	const double t2 = hypot(b[0]-a[0], b[1]-a[1]);
	const double v[2] = {
		0.5*(b[0]-a[0]),
		0.5*(b[1]-a[1])
	};
	const double n[2] = { (b[1]-a[1])/t2, (a[0]-b[0])/t2 };
	const double alpha = (2./3.) * (1.+g) * (1.-g);
	const double beta = (2./3.) * g * t2;
	ap[0] = a[0] + alpha*v[0] + beta*n[0];
	ap[1] = a[1] + alpha*v[1] + beta*n[1];
	bp[0] = b[0] - alpha*v[0] + beta*n[0];
	bp[1] = b[1] - alpha*v[1] + beta*n[1];
}

double geom_arc_subdivide(
	const double a[2], const double b[2], double g,
	double p[2]
){
	double gp = qsolve(1, g, 1.);
	if(NULL != p){
		p[0] = 0.5*a[0] + 0.5*b[0] + 0.5*g*(b[1]-a[1]);
		p[1] = 0.5*a[0] + 0.5*b[0] + 0.5*g*(a[0]-b[0]);
	}
	if(g > 0 && gp < 0){ return -1./gp; }
	else if(g < 0 && gp > 0){ return -1./gp; }
	else{ return gp; }
}

/* Classify the quadrants:
 *
 *    1 | 0
 *   ---+---
 *    2 | 3
 * If cw is false, then quadrants are CW side inclusive and CCW side exclusive
 * otherwise, the quadrants are CW side exclusive and CCW side inclusive.
 */
static int quadrant_classify(int cw, const double v[2]){
	if(!cw){
		if(v[0] > 0 && v[1] >= 0){
			return 0;
		}else if(v[0] <= 0 && v[1] > 0){
			return 1;
		}else if(v[0] < 0 && v[1] <= 0){
			return 2;
		}else{ /* v[0] >= 0 && v[1] < 0 */
			return 3;
		}
	}else{
		if(v[0] >= 0 && v[1] > 0){
			return 0;
		}else if(v[0] < 0 && v[1] >= 0){
			return 1;
		}else if(v[0] <= 0 && v[1] < 0){
			return 2;
		}else{ /* v[0] > 0 && v[1] <= 0 */
			return 3;
		}
	}
}

int geom_arc_split_monotone(
	const double a[2], const double b[2], double g,
	double *pg
){
	if(0 == g){
		pg[0] = a[0]; pg[1] = a[1]; pg[2] = 0;
		pg[3] = b[0]; pg[4] = b[1];
		return 1;
	}
	if(fabs(g) > G_QUARTER){ /* more than a quarter circle */
		/* We must have 2, 3, 4, or 5 segments */
		/* find center */
		const int dirlist[4][2] = {
			{ 0,1 }, { -1,0 }, { 0,-1 }, { 1,0 }
		};
		const double gg = (1+g)*(1-g) / (4*g);
		const double c[2] = {
			(0.5*a[0]+0.5*b[0]) + (a[1]-b[1])*gg,
			(0.5*a[1]+0.5*b[1]) + (b[0]-a[0])*gg
		};
		const double r = 0.25*(1./g + g) * hypot(b[0]-a[0],b[1]-a[1]);
		double t0[2], t1[2], p[2];
		int q0, q1, n = 0;
		const int cw = (g > 0 ? 0 : 1);
		const int inc = (cw ? 3 : 1);
		geom_arc_param(a, b, g, 0, p, t0);
		geom_arc_param(a, b, g, 1, p, t1);
//fprintf(stderr, "inc = %d, g = %g\n", inc, g);
		q0 = quadrant_classify( cw, t0);
		q1 = quadrant_classify(!cw, t1);
		
		pg[0] = a[0]; pg[1] = a[1];
		while(q0 != q1 || n < 1){
			const int ax = (q0+4-cw+3)%4;
			double t, gs;
//fprintf(stderr, "q0 = %d, q1 = %d, ax = %d\n", q0, q1, ax);
			++n;
			pg[3*n+0] = c[0] + dirlist[ax][0] * r;
			pg[3*n+1] = c[1] + dirlist[ax][1] * r;
			t = 0.5*hypot(pg[3*n+0]-pg[3*(n-1)+0], pg[3*n+1]-pg[3*(n-1)+1]);
			gs = qsolve(0, t, -r);
			if((g > 0 && gs < 0) || (g < 0 && gs > 0)){ gs = -1./gs; }
//fprintf(stderr, "gs = %g\n", gs);
			pg[3*(n-1)+2] = gs;
			q0 = (q0+inc)%4;
		}
		++n;
		pg[3*n+0] = b[0];
		pg[3*n+1] = b[1];
		{
			double t = 0.5*hypot(pg[3*n+0]-pg[3*(n-1)+0], pg[3*n+1]-pg[3*(n-1)+1]);
			pg[3*(n-1)+2] = qsolve(0, t, -r);
		}
//fprintf(stderr, "n = %d\n", n);
		return n;
	}else{ /* less than a quarter circle */
		/* We must have 1 or 2 segments */
		int cs = -1;
		double u[2], v[2], p[2], s;
		/* First obtain the tangents u and v */
		geom_arc_param(a, b, g, 0, p, u);
		geom_arc_param(a, b, g, 1, p, v);
		/* If any coordinate of u and v have different signs,
		 * then we have to make a split.
		 */
		if(u[0] > 0 && v[0] < 0){ cs = 0; }
		else if(u[0] < 0 && v[0] > 0){ cs = 1; }
		else if(u[1] > 0 && v[1] < 0){ cs = 2; }
		else if(u[1] < 0 && v[1] > 0){ cs = 3; }

		if(cs >= 0){
			static const double d[4][2] = {
				{ 1, 0 }, {-1, 0 }, { 0, 1 }, { 0,-1 }
			};
			geom_arc_extremum(a, b, g, d[cs], &pg[3]);
			s = geom_arc_unparam(a, b, g, &pg[3]);
			tanfrac(g, s, &pg[2], &pg[5]);
			pg[0] = a[0]; pg[1] = a[1];
			pg[6] = b[0]; pg[7] = b[1];
			return 2;
		}else{
			pg[0] = a[0]; pg[1] = a[1]; pg[2] = g;
			pg[3] = b[0]; pg[4] = b[1];
			return 1;
		}
	}
}

