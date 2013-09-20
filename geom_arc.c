#include <stddef.h>
#include <math.h>
#include <stdio.h>

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

/* returns sin(a*x)/sin(x) for -1 <= a <= 1 */
static double sinratio(double a, double x){
	double sign = 1;
	if(0 == x){ return a; }
	else if(x < 0){ x = -x; }
	
	if(0 == a){ return 0; }
	else if(a < 0){ a = -a; sign = -1; }
	if(1 == a){ return sign; }
	if(M_PI_2 == x){ return sign*sin(a*M_PI_2); }
	if(x < 1e-8){ return sign*a; }
	{
		const double x2 = x*x;
		const double a2 = a*a;
		if(x < 1e-4){ return sign*a*(1. + (1./6.)*(1.-a2)*x2); }
		else{
			return sign*sin(a*x) / sin(x);
		}
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
		const double sr0 = sinratio(   s, theta);
		const double sr1 = sinratio(1.-s, theta);
		const double lv = 4./(1.+g*g) * sr0 * sr1;
		const double sr = sinratio(2.*s - 1., theta);
		const double la = 0.5 * (1. - sr) - 0.5*lv;
		const double lb = 0.5 * (1. + sr) - 0.5*lv;
		p[0] = la*a[0] + lb*b[0] + lv*(m[0]+mv[0]);
		p[1] = la*a[1] + lb*b[1] + lv*(m[1]+mv[1]);
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
			double alphap = alpha * gg1/(2.*g);
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
