#include <stdlib.h>
#include <math.h>
#include <Cgeom/geom_la.h>
#include <Cgeom/geom_predicates.h>
#include <float.h>

float  geom_polygon_area2f(unsigned int n, const float  *v){
	unsigned int p, q;
	float area = 0.f;
	for(p=n-1, q=0; q < n; p = q++){
		area += v[2*p+0]*v[2*q+1] - v[2*q+0]*v[2*p+1];
	}
	return 0.5f*area;
}
double geom_polygon_area2d(unsigned int n, const double *v){
	unsigned int p, q;
	double area = 0.;
	for(p=n-1, q=0; q < n; p = q++){
		area += v[2*p+0]*v[2*q+1] - v[2*q+0]*v[2*p+1];
	}
	return 0.5*area;
}

/* License for point-in-polygon code:

Copyright (c) 1970-2003, Wm. Randolph Franklin

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

	1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
	2. Redistributions in binary form must reproduce the above copyright notice in the documentation and/or other materials provided with the distribution.
	3. The name of W. Randolph Franklin may not be used to endorse or promote products derived from this Software without specific prior written permission. 

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 

*/

/* Returns 1 if inside, 0 otherwise */
int geom_polygon_inside2f(unsigned int nv, const float  *v, const float  p[2]){
	unsigned int i, j;
	int c = 0;
	for(i = 0, j = nv-1; i < nv; j = i++){
		if(
			((v[2*i+1] > p[1]) != (v[2*j+1]>p[1])) &&
			(p[0] < (v[2*j+0]-v[2*i+0]) * (p[1]-v[2*i+1]) / (v[2*j+1]-v[2*i+1]) + v[2*i+0])
		){
			c = !c;
		}
	}
	return c;
}
int geom_polygon_inside2d(unsigned int nv, const double *v, const double p[2]){
	unsigned int i, j;
	int c = 0;
	for(i = 0, j = nv-1; i < nv; j = i++){
		if(
			((v[2*i+1] > p[1]) != (v[2*j+1]>p[1])) &&
			(p[0] < (v[2*j+0]-v[2*i+0]) * (p[1]-v[2*i+1]) / (v[2*j+1]-v[2*i+1]) + v[2*i+0])
		){
			c = !c;
		}
	}
	return c;
}

/* Returns a normal vector for a point on (near) the polygon's boundary */
void geom_polygon_normal2f(unsigned int nv, const float  *v, const float  p[2], float  n[2]){
	unsigned int i, j;
	float mindist = FLT_MAX;
	unsigned int idist = -1;
	n[0] = 0.f; n[1] = 0.f;
	for(j = 0, i = nv-1; j < nv; i = j++){
		/* compute distance from r to segment */
		float u[2] = {
			v[2*j+0] - v[2*i+0],
			v[2*j+1] - v[2*i+1]
		};
		float vp[2] = {
			p[0] - v[2*i+0],
			p[1] - v[2*i+1]
		};
		const float ulen = geom_norm2f(u);
		const float u2 = u[0]*u[0] + u[1]*u[1];
		const float vpu = vp[0]*u[0] + vp[1]*u[1];
		const float t = vpu/u2;
		const float off[2] = {vp[0] - t*u[0], vp[1] - t*u[1] };
		u[0] /= ulen; u[1] /= ulen;
		float dist = geom_norm2f(off);
		int curind = -1;
		if(t < 0.f){
			dist = geom_norm2f(vp);
			curind = i;
		}else if(t > 1.f){
			vp[0] = p[0] - v[2*j+0];
			vp[1] = p[1] - v[2*j+1];
			dist = geom_norm2f(vp);
			curind = j;
		}
		if(-1 != idist && curind == idist){
			n[0] = 0.5f*( u[1] + n[0]);
			n[1] = 0.5f*(-u[0] + n[1]);
			idist = -1;
		}else if(dist < mindist){
			mindist = dist;
			n[0] = u[1];
			n[1] = -u[0];
			idist = curind;
		}
	}
	geom_normalize2f(n);
}
void geom_polygon_normal2d(unsigned int nv, const double *v, const double p[2], double n[2]){
	unsigned int i, j;
	double mindist = DBL_MAX;
	unsigned int idist = -1;
	n[0] = 0.0; n[1] = 0.0;
	for(j = 0, i = nv-1; j < nv; i = j++){
		/* compute distance from r to segment */
		double u[2] = {
			v[2*j+0] - v[2*i+0],
			v[2*j+1] - v[2*i+1]
		};
		double vp[2] = {
			p[0] - v[2*i+0],
			p[1] - v[2*i+1]
		};
		const double ulen = geom_norm2d(u);
		const double u2 = u[0]*u[0] + u[1]*u[1];
		const double vpu = vp[0]*u[0] + vp[1]*u[1];
		const double t = vpu/u2;
		const double off[2] = {vp[0] - t*u[0], vp[1] - t*u[1] };
		u[0] /= ulen; u[1] /= ulen;
		double dist = geom_norm2d(off);
		int curind = -1;
		if(t < 0.0){
			dist = geom_norm2d(vp);
			curind = i;
		}else if(t > 1.0){
			vp[0] = p[0] - v[2*j+0];
			vp[1] = p[1] - v[2*j+1];
			dist = geom_norm2d(vp);
			curind = j;
		}
		if(-1 != idist && curind == idist){
			n[0] = 0.5f*( u[1] + n[0]);
			n[1] = 0.5f*(-u[0] + n[1]);
			idist = -1;
		}else if(dist < mindist){
			mindist = dist;
			n[0] = u[1];
			n[1] = -u[0];
			idist = curind;
		}
	}
	geom_normalize2d(n);
}

int geom_convex_inside2f(unsigned int np, const float  *p, const float  r[2]){
	unsigned int i;
	for(i = 0; i < np; ++i){
		if(p[3*i+0] * r[0] + p[3*i+1] * r[1] > p[3*i+2]){ return 0; }
	}
	return 1;
}
int geom_convex_inside2d(unsigned int np, const double *p, const double r[2]){
	unsigned int i;
	for(i = 0; i < np; ++i){
		if(p[3*i+0] * r[0] + p[3*i+1] * r[1] > p[3*i+2]){ return 0; }
	}
	return 1;
}

void geom_convex_normal2f(unsigned int np, const float  *p, const float  r[2], float  n[2]){
	unsigned int i;
	float maxdist = -FLT_MAX;
	n[0] = 0.f; n[1] = 0.f;
	for(i = 0; i < np; ++i){
		float d = p[3*i+0] * r[0] + p[3*i+1] * r[1] - p[3*i+2];
		if(d >= maxdist){
			maxdist = d;
			n[0] = p[3*i+0];
			n[1] = p[3*i+1];
		}
	}
	geom_normalize2f(n);
}
void geom_convex_normal2d(unsigned int np, const double *p, const double r[2], double n[2]){
	unsigned int i;
	double maxdist = -DBL_MAX;
	n[0] = 0.; n[1] = 0.;
	for(i = 0; i < np; ++i){
		double d = p[3*i+0] * r[0] + p[3*i+1] * r[1] - p[3*i+2];
		if(d >= maxdist){
			maxdist = d;
			n[0] = p[3*i+0];
			n[1] = p[3*i+1];
		}
	}
	geom_normalize2d(n);
}

int geom_convex_inside3f(unsigned int np, const float  *p, const float  r[3]){
	unsigned int i;
	for(i = 0; i < np; ++i){
		if(p[4*i+0] * r[0] + p[4*i+1] * r[1] + p[4*i+2] * r[2] > p[4*i+3]){ return 0; }
	}
	return 1;
}
int geom_convex_inside3d(unsigned int np, const double *p, const double r[3]){
	unsigned int i;
	for(i = 0; i < np; ++i){
		if(p[4*i+0] * r[0] + p[4*i+1] * r[1] + p[4*i+2] * r[2] > p[4*i+3]){ return 0; }
	}
	return 1;
}

void geom_convex_normal3f(unsigned int np, const float  *p, const float  r[3], float  n[3]){
	unsigned int i;
	float maxdist = -FLT_MAX;
	n[0] = 0.f; n[1] = 0.f; n[2] = 0.f;
	for(i = 0; i < np; ++i){
		float d = p[4*i+0] * r[0] + p[4*i+1] * r[1] + p[4*i+2] * r[2] - p[4*i+3];
		if(d >= maxdist){
			maxdist = d;
			n[0] = p[4*i+0];
			n[1] = p[4*i+1];
			n[2] = p[4*i+2];
		}
	}
	geom_normalize3f(n);
}
void geom_convex_normal3d(unsigned int np, const double *p, const double r[3], double n[3]){
	unsigned int i;
	double maxdist = -DBL_MAX;
	n[0] = 0.; n[1] = 0.; n[2] = 0.;
	for(i = 0; i < np; ++i){
		double d = p[4*i+0] * r[0] + p[4*i+1] * r[1] + p[4*i+2] * r[2] - p[4*i+3];
		if(d >= maxdist){
			maxdist = d;
			n[0] = p[4*i+0];
			n[1] = p[4*i+1];
			n[2] = p[4*i+2];
		}
	}
	geom_normalize3d(n);
}

/* begin stuff for geom_convex_bound3d */

typedef struct Settings_t {
	/* For regularization. Minimum value of abs(D_ii) in the kkt D factor */
	double kkt_reg;
	double resid_tol;
	double eps;
	int max_iters;
	int refine_steps;
	/* Better start obviates the need for s_init and z_init */
	double s_init;
	double z_init;
} Settings;

typedef struct Workspace_t {
	double *s_inv;   /* size  n   */
	double *s_inv_z; /* size  n   */
	double *rhs;     /* size 2n+3 */
	double *x;       /* size 2n+3 */
	double *lhs_aff; /* size 2n+3 */
	double *lhs_cc;  /* size 2n+3 */
	double *buffer;  /* size 2n+3 */
	double *d_inv;   /* size 2n+3 */
	double *L;       /* size 4n+3 */
} Workspace;

/* computes y = A * x */
static void dmmv(unsigned int m, unsigned int n, const double *A, unsigned int ldA, const double *x, double *y){
	unsigned int i, j;
	for(i = 0; i < m; ++i){
		y[i] = 0;
		for(j = 0; j < n; ++j){
			y[i] -= A[i+j*ldA] * x[j];
		}
	}
}
/* computes y = A^T * x */
static void dmmTv(unsigned int m, unsigned int n, const double *A, unsigned int ldA, const double *x, double *y){
	unsigned int i, j;
	for(j = 0; j < n; ++j){
		y[j] = 0;
		for(i = 0; i < m; ++i){
			y[j] -= x[i] * A[i+j*ldA];
		}
	}
}

/* x = inv(L') * inv(D) * inv(L) * b */
static void ldl_solve(unsigned int n, const double *L, const double *d_inv, const double *b, double *x){
	unsigned int i, j;
	
	/* Forward substitution */
	for(i = 0; i < n; ++i){
		x[i] = b[i];
	}
	for(i = 0; i < n; ++i){
		x[n+i] = b[n+i] - L[i]*x[i];
	}
	for(i = 0; i < 3; ++i){
		x[2*n+i] = b[2*n+i];
		for(j = 0; j < n; ++j){
			x[2*n+i] -= L[n+i+j*3]*x[n+j];
		}
	}
	x[2*n+1] -= x[2*n+0]*L[4*n+0];
	x[2*n+2] -= x[2*n+0]*L[4*n+1] + x[2*n+1]*L[4*n+2];
	/* Diagonal scaling */
	for(i = 0; i < n+n+3; i++){
		x[i] *= d_inv[i];
	}
	/* Backward substitution */
	x[2*n+1] -= x[2*n+2]*L[4*n+2];
	x[2*n+0] -= x[2*n+2]*L[4*n+1] + x[2*n+1]*L[4*n+0];
	for(i = 0; i < n; ++i){
		for(j = 0; j < 3; ++j){
			x[n+i] -= L[n+j+i*3]*x[2*n+j];
		}
	}
	for(i = 0; i < n; ++i){
		x[i] -= L[i]*x[n+i];
	}
}

/*     [    I                             ]
 * L = [ inv(sz)          I               ]
 *     [    0     A*inv(b33-inv(sz))   I  ]
 *
 * D = diag([    sz       b33-inv(sz)    -A*inv(b33-inv(sz))*A'  ])
 */
static void ldl_factor(unsigned int n, const double *s_inv_z, double b33, const double *A, unsigned int ldA, double reg, double *L, double *d_inv){
	unsigned int i;
	
	/* Block wise factorization */
	for(i = 0; i < n; ++i){
		double dii = s_inv_z[i];
		if(fabs(dii) < reg){ dii = reg; }
		d_inv[i] = 1./dii;
		L[i] = d_inv[i];
	}
	for(i = 0; i < n; ++i){
		double dii = b33 - d_inv[i];
		if(fabs(dii) < reg){ dii = reg; }
		d_inv[n+i] = 1./dii;
		L[n+3*i+0] = A[0+i*ldA] * d_inv[n+i];
		L[n+3*i+1] = A[1+i*ldA] * d_inv[n+i];
		L[n+3*i+2] = A[2+i*ldA] * d_inv[n+i];
	}
	/* Compute lower 3x3 diagonal block */
	double d33[9]; /* use d33[3] for dii, d33[6] for L32D2 */
	for(i = 0; i < 3; ++i){
		unsigned j;
		for(j = 0; j <= i; ++j){
			unsigned k;
			d33[i+j*3] = 0;
			for(k = 0; k < n; ++k){
				d33[i+j*3] -= L[n+i+k*3] * A[j+k*ldA];
			}
		}
	}
	/* Perform LDL on 3x3 block */
	d33[3] = d33[0];
	if(fabs(d33[3]) < reg){ d33[3] = reg; }
	d_inv[2*n+0] = 1./d33[3];
	L[4*n+0] = d33[1] * d_inv[2*n+0];
	L[4*n+1] = d33[2] * d_inv[2*n+0];
	
	d33[3] = d33[4] - d33[1] * L[4*n+0];
	if(fabs(d33[3]) < reg){ d33[3] = reg; }
	d_inv[2*n+1] = 1./d33[3];
	d33[6] = d33[5] - d33[1] * L[4*n+1];
	L[4*n+2] = d33[6] * d_inv[2*n+1];
	
	d33[3] = d33[8] - d33[2]*L[4*n+1] - L[4*n+2]*d33[6];
	if(fabs(d33[3]) < reg){ d33[3] = reg; }
	d_inv[2*n+2] = 1./d33[3];
}

/* y = KKT*x.
 *       [ sz  I     ]
 * KKT = [ I  b33 A' ]
 *       [     A  0  ]
 */
static void multKKT(unsigned int n, const double *s_inv_z, const double *A, unsigned int ldA, const double *x, double *y){
	unsigned int i, j;
	for(i = 0; i < n; ++i){
		y[i] = s_inv_z[i]*x[i] + x[n+i];
	}
	for(i = 0; i < n; ++i){
		y[n+i] = x[i];
		for(j = 0; j < 3; ++j){
			y[n+i] += A[j+i*ldA]*x[2*n+j];
		}
	}
	for(i = 0; i < 3; ++i){
		y[2*n+i] = 0;
		for(j = 0; j < n; ++j){
			y[2*n+i] += A[i+j*ldA]*x[n+j];
		}
	}
}

static void refine(unsigned int np, const double *p, const Settings *settings, Workspace *work, const double *b, double *x){
	const unsigned int n23 = 2*np+3;
	unsigned int i, j;

	double *residual = work->buffer;
	double *y = work->rhs;
	for(j = 0; j < settings->refine_steps; j++) {
		multKKT(np, work->s_inv_z, p, 4, x, residual);
		for(i = 0; i < n23; i++) {
			residual[i] -= b[i];
		}
		/* Solve to find new_var = KKT \ (target - A*var) */
		ldl_solve(np, work->L, work->d_inv, residual, y);
		/* Update var += new_var, or var += KKT \ (target - A*var) */
		for(i = 0; i < n23; i++) {
			x[i] -= y[i];
		}
	}
}

static void init_vars2(unsigned int np, const double *p, const double dir[3], const Settings *settings, Workspace *work){
	/* Calculates a better starting point, using a similar approach to CVXOPT */
	unsigned int i;
	double *x, /* *s, */ *z;
	double alpha;

	/* Make sure sinvz is 1 to make hijacked KKT system ok */
	for(i = 0; i < np; i++){
		work->s_inv_z[i] = 1;
	}
	ldl_factor(np, work->s_inv_z, -1., p, 4, settings->kkt_reg, work->L, work->d_inv);
	{ /* Fill rhs with (0, h, c) */
		for(i = 0; i < np; i++){
			work->rhs[i] = 0;
		}
		for(i = 0; i < np; i++){
			work->rhs[np+i] = p[4*i+3];
		}	
		for(i = 0; i < 3; i++){
			work->rhs[2*np+i] = dir[i];
		}
	}
	/* Borrow work.lhs_aff for the solution */
	ldl_solve(np, work->L, work->d_inv, work->rhs, work->lhs_aff);
	/* Don't do any refinement for now. Precision doesn't matter too much */
	/* s = work->lhs_aff; */
	z = work->lhs_aff + np;
	x = work->lhs_aff + 2*np;

	/* Just set x and y as is */
	for(i = 0; i < 3; i++){ work->x[i] = x[i]; }

	/* Now complete the initialization. Start with s
	 * Must have alpha > max(z)
	 */
	alpha = -DBL_MAX;
	for(i = 0; i < np; i++){
		if(alpha < z[i]){
			alpha = z[i];
		}
	}
	if(alpha < 0){
		for(i = 0; i < np; i++){
			work->x[3+i] = -z[i];
		}
	}else{
		alpha += 1;
		for(i = 0; i < np; i++){
			work->x[3+i] = -z[i] + alpha;
		}
	}
	/* Now initialize z
	 * Now must have alpha > max(-z)
	 */
	alpha = -DBL_MAX;
	for(i = 0; i < np; i++){
		if(alpha < -z[i]){
			alpha = -z[i];
		}
	}
	if(alpha < 0){
		for(i = 0; i < np; i++){
			work->x[np+3+i] = z[i];
		}
	}else{
		alpha += 1;
		for(i = 0; i < np; i++){
			work->x[np+3+i] = z[i] + alpha;
		}
	}
}

/* Solves the following linear program:
 *  min -dir' * r
 *  s.t. p' * [r;1] <= 0
 */
int geom_convex_bound3d(unsigned int np, const double *p, const double dir[3], double r[3], double *wksp){
	const unsigned int n23 = 2*np+3;
	Settings settings;
	unsigned int i;
	int iter;

	double *dx, *ds, *dz;
	double minval;
	double alpha;
	
	settings.resid_tol = 1e-6;
	settings.eps = 1e-4;
	settings.max_iters = 25;
	settings.refine_steps = 1;
	settings.s_init = 1;
	settings.z_init = 1;
	settings.kkt_reg = 1e-7;

	double *workalloc = wksp;
	if(NULL == wksp){
		workalloc = (double*)malloc(sizeof(double) * (7*n23+4*np));
	}
	Workspace work;
	work.s_inv = workalloc;
	work.s_inv_z = work.s_inv + np;
	work.rhs = work.s_inv_z + np;
	work.x = work.rhs + n23;
	work.lhs_aff = work.x + n23;
	work.lhs_cc = work.lhs_aff + n23;
	work.buffer = work.lhs_cc + n23;
	work.d_inv = work.buffer + n23;
	work.L = work.d_inv + n23;
	
	/* previously in workspace */
	int converged = 0;
	double gap;
	/* double optval; */
	double ineq_resid_squared;

	/* first 3 elements of work.x are the actual variables
	 * s is &work.x[3], size np
	 * z is &work.x[np+3], size np
	 */
	/*printf("iter     objv        gap       |Gx+s-h|    step\n"); */
	/*
	for(i = 0; i < 3; i++){
		work.x[i] = 0;
	}
	for(i = 0; i < np; i++){
		work.x[3+i] = (p[4*i+3] > 0) ? p[4*i+3] : settings.s_init;
	}
	for(i = 0; i < np; i++){
		work.x[np+3+i] = settings.z_init;
	}
	*/
	init_vars2(np, p, dir, &settings, &work);
	
	for(iter = 0; iter < settings.max_iters; iter++){
		for(i = 0; i < np; i++){
			work.s_inv[i] = 1.0 / work.x[3+i];
			work.s_inv_z[i] = work.s_inv[i]*work.x[np+3+i];
		}

		ldl_factor(np, work.s_inv_z, 0., p, 4, settings.kkt_reg, work.L, work.d_inv);
		{ /* Affine scaling directions */
			/* r1 = -z */
			for (i = 0; i < np; i++){ work.rhs[i] = -work.x[np+3+i]; }
			/* r2 = -Gx - s + h */
			dmmTv(3, np, p, 4, work.x, &work.rhs[np]);
			for(i = 0; i < np; i++){ work.rhs[np+i] += -work.x[3+i] + p[4*i+3]; }
			/* r3 = -A^Ty - G^Tz - Px - q */
			dmmv(3, np, p, 4, &work.x[np+3], &work.rhs[2*np]);
			for(i = 0; i < 3; i++){ work.rhs[2*np+i] += dir[i]; }
		}
		ldl_solve(np, work.L, work.d_inv, work.rhs, work.lhs_aff);
		refine(np, p, &settings, &work, work.rhs, work.lhs_aff);
		{ /* Centering plus corrector directions */
			double *ds_aff = work.lhs_aff, *dz_aff = work.lhs_aff + np;
			double mu = 0;
			double alpha;
			double sigma = 0;
			double smu;
			double minval = 0;

			for(i = 0; i < np; i++){
				mu += work.x[3+i]*work.x[np+3+i];
			}

			/* Find min(min(ds./s), min(dz./z)) */
			for(i = 0; i < np; i++){
				if(ds_aff[i] < minval*work.x[3+i]){
					minval = ds_aff[i]/work.x[3+i];
				}
			}
			for(i = 0; i < np; i++){
				if(dz_aff[i] < minval*work.x[np+3+i]){
					minval = dz_aff[i]/work.x[np+3+i];
				}
			}

			/* Find alpha */
			if(-1 < minval){
				alpha = 1;
			}else{
				alpha = -1/minval;
			}

			sigma = 0;
			for(i = 0; i < np; i++){
				sigma += (work.x[3+i] + alpha*ds_aff[i])*
					(work.x[np+3+i] + alpha*dz_aff[i]);
			}
			sigma /= mu;
			sigma = sigma*sigma*sigma;

			mu *= 0.0909090909090909;
			smu = sigma*mu;

			/* Fill-in the rhs */
			for(i = 0; i < np; i++){
				work.rhs[i] = work.s_inv[i]*(smu - ds_aff[i]*dz_aff[i]);
			}
			for(i = np; i < n23; i++){
				work.rhs[i] = 0;
			}
		}
		ldl_solve(np, work.L, work.d_inv, work.rhs, work.lhs_cc);
		refine(np, p, &settings, &work, work.rhs, work.lhs_cc);

		/* Add the two together and store in aff */
		for(i = 0; i < n23; i++){
			work.lhs_aff[i] += work.lhs_cc[i];
		}

		/* Rename aff to reflect its new meaning */
		ds = work.lhs_aff;
		dz = work.lhs_aff + np;
		dx = work.lhs_aff + 2*np;
		/* Find min(min(ds./s), min(dz./z)) */
		minval = 0;
		for(i = 0; i < np; i++){
			if(ds[i] < minval*work.x[3+i]){
				minval = ds[i]/work.x[3+i];
			}
		}
		for(i = 0; i < np; i++){
			if(dz[i] < minval*work.x[np+3+i]){
				minval = dz[i]/work.x[np+3+i];
			}
		}

		/* Find alpha */
		if(-0.99 < minval){
			alpha = 1;
		}else{
			alpha = -0.99/minval;
		}

		/* Update the primal and dual variables */
		for(i = 0; i < 3; i++){
			work.x[i] += alpha*dx[i];
		}
		for(i = 0; i < np; i++){
			work.x[3+i] += alpha*ds[i];
		}
		for(i = 0; i < np; i++){
			work.x[np+3+i] += alpha*dz[i];
		}
		{
			gap = 0;
			for(i = 0; i < np; i++){
				gap += work.x[np+3+i]*work.x[3+i];
			}
		}
		{ /* Calculate the norm ||-Gx - s + h|| */
			/* Find -Gx */
			dmmTv(3, np, p, 4, work.x, work.buffer);
			/* Add -s + h */
			for(i = 0; i < np; i++){
				work.buffer[i] += -work.x[3+i] + p[4*i+3];
			}
			/* Now find the squared norm */
			ineq_resid_squared = 0;
			for(i = 0; i < np; i++){
				ineq_resid_squared += work.buffer[i]*work.buffer[i];
			}
		}

		/* optval = -dir[0]*work.x[0]-dir[1]*work.x[1]-dir[2]*work.x[2]; */
		/*
		printf("%3d   %10.3e  %9.2e  %9.2e  % 6.4f\n",
			iter+1, optval, gap,
			sqrt(ineq_resid_squared), alpha
		);
		*/

		/* Test termination conditions. Requires optimality, and satisfied constraints */
		if((gap < settings.eps)
		&& (ineq_resid_squared <= settings.resid_tol*settings.resid_tol)
		){
			converged = 1;
			iter++;
			break;
		}
	}
	r[0] = work.x[0];
	r[1] = work.x[1];
	r[2] = work.x[2];
	if(NULL == wksp){
		free(workalloc);
	}
	return !converged;
}

static int triangle_contains(
	const double org[2], /* triangle vertices are {org,org+u,org+v}, in CCW orientation */
	const double u[2],
	const double v[2],
	const double p[2] /* query point */
){
	if(geom_orient2d(org,u,p) >= 0){
		if(geom_orient2d(org,v,p) <= 0){
			/* treat org as origin, we have points u and v, just need p-org */
			double x[2] = {p[0] - org[0], p[1] - org[1]};
			if(geom_orient2d(u,v,x) >= 0){
				return 1;
			}
		}
	}
	return 0;
}
/* An n-sided polygon always has n-2 triangles in its triangulation. */
int geom_polygon_triangulate2d(
	unsigned int n, const double *v, /* the polygon */
	unsigned int *t /* length 3*(n-2), stores the triangle as triples of vertex indices into v */
){
	/* t is used to store a working copy of the currently clipped polygon vertex index, size 3 <= m <= n.
	 * t must also store the resulting triangles, size 3*(n-m).
	 * The total size of these two is 3*n-2*m, which must be <= 3*(n-2). This is true for m >= 3.
	 * Therefore, we keep the working copy at the very end of t, and add triangles to the beginning.
	 */
	unsigned int *V; /* pointer to start of working copy */
	unsigned int nv; /* size of V */
	unsigned int i;
	int count;
	int tc = 0; /* number of triangles currently in t */
	
	if(n < 3){ return -1; }
	if(NULL == v){ return -2; }
	if(NULL == t){ return -3; }

	/* Make a copy of all the vertices */
	V = t + 2*n-6;
	nv = n;
	for(i = 0; i < n; ++i){ V[i] = i; }
	
	count = 2*nv;
	for(i = nv-1; nv > 2; ){
		/*
		fprintf(stderr, "V:");
		for(j = 0; j < nv; ++j){
			fprintf(stderr, " %d", V[j]);
		}fprintf(stderr, "\n");
		fprintf(stderr, "count = %d\n", count);
		*/
		int u, w;
		if(0 >= (count--)){ return 1; } /* bad polygon */
		/* get 3 consecutive vertices */
		u = i; /* if(nv <= u){ u = 0; } */ /* prev */
		i = u+1; if(nv <= i){ i = 0; } /* mid */
		w = i+1; if(nv <= w){ w = 0; } /* next */

		/* Can clip the ear? */
		int can_clip = 1;
		{
			int p;
			double tri_a[2] = {v[2*V[i]+0] - v[2*V[u]+0], v[2*V[i]+1] - v[2*V[u]+1]};
			double tri_b[2] = {v[2*V[w]+0] - v[2*V[u]+0], v[2*V[w]+1] - v[2*V[u]+1]};
			if(geom_orient2d(&v[2*V[u]], &v[2*V[i]], &v[2*V[w]]) < 0){
				can_clip = 0;
			}else{
				/* if the u-i-w triangle contains any other vertex, can't clip. */
				for(p = 0; p < nv; ++p){
					if((p == u) || (p == i) || (p == w)){ continue; }
					if(triangle_contains(&v[2*V[u]],tri_a,tri_b, &v[2*V[p]])){ can_clip = 0; break; }
				}
			}
		}

		/* Clip off the ear */
		if(can_clip){
			unsigned int tri[3] = {V[u], V[i], V[w]};
			/* erase vertex i */
			while(i > 0){
				V[i] = V[i-1];
				--i;
			}
			++V; --nv; count = 2*nv;
			/* Add the new triangle */
			t[3*tc+0] = tri[0];
			t[3*tc+1] = tri[1];
			t[3*tc+2] = tri[2];
			++tc;
		}
	}
	return 0;
}

static int dsign(double x){
	if(0 == x){ return 0; }
	else if(x > 0){ return 1; }
	else{ return -1; }
}
static int segments_intersect(const double a[2], const double b[2], const double c[2], const double d[2], double *x, double *t){
	double c0, c1, c2, c3;
	c0 = geom_orient2d(a,b,c);
	c1 = geom_orient2d(a,b,d);
	if(0 == c0 && 0 == c1){ return 0; } /* collinear -> no intersection */
	if(dsign(c0) != dsign(c1)){
		c2 = geom_orient2d(c, d, a);
		c3 = geom_orient2d(c, d, b);
		if(dsign(c2) != dsign(c3)){
			if(NULL != x){
				double r;
				c1 = c0-c1; c3 = c2-c3;
				if(fabs(c1) > fabs(c3)){
					r = c0/c1;
					x[0] = c[0] + r*(d[0]-c[0]);
					x[1] = c[1] + r*(d[1]-c[1]);
				}else{
					r = c2/c3;
					x[0] = a[0] + r*(b[0]-a[0]);
					x[1] = a[1] + r*(b[1]-a[1]);
				}
				if(NULL != t){ *t = r; }
			}
			return 1;
		}else{ return 0; }
	}else{ return 0; }
}
int geom_convex_polygon_intersection2d(
	unsigned int n, /* n >= 3 */
	const double *P,
	unsigned int m, /* m >= 3 */
	const double *Q,
	unsigned int *ni, /* on input, size of Pi, on output, numer of points in Pi */
	double *Pi /* output intersection polygon */
){
	int i, j;
	if(n < 3){ return -1; }
	if(NULL == P){ return -2; }
	if(m < 3){ return -3; }
	if(NULL == Q){ return -4; }
	if(NULL == ni){ return -5; }
	if(NULL == Pi){ return -6; }

	/* Implementation of:
	 * "A new linear algorithm for intersecting convex polygons"
	 * Joseph O'Rourke, Chi-Bin Chien, Thomas Olson, and David Naddor
	 * Computer Graphics and Image Processing 19, pp. 384-391 (1982)
	 */
	
	const unsigned int nPi = *ni; *ni = 0;
	
	int ip = 1, iq = 1;
	int ipp = 0, iqp = 0; /* prev of ip and iq */
	char inside = ' ';
	/* record first intersection */
	int first_xsected = 0;
	int ipf = n, iqf = m;
	int first_iter = 0;
	
	int Pi_full = 0;
	int iter;

	/* First, a bounding box check */
	{
		int iP_x_min = 0, iP_x_max = 0, iP_y_min = 0, iP_y_max = 0;
		int iQ_x_min = 0, iQ_x_max = 0, iQ_y_min = 0, iQ_y_max = 0;
		for(i = 1; i < n; ++i){
			if(P[2*i+0] < P[2*iP_x_min+0]){ iP_x_min = i; }
			if(P[2*i+1] < P[2*iP_y_min+1]){ iP_y_min = i; }
			if(P[2*i+0] > P[2*iP_x_max+0]){ iP_x_max = i; }
			if(P[2*i+1] > P[2*iP_y_max+1]){ iP_y_max = i; }
		}
		for(i = 1; i < m; ++i){
			if(Q[2*i+0] < Q[2*iQ_x_min+0]){ iQ_x_min = i; }
			if(Q[2*i+1] < Q[2*iQ_y_min+1]){ iQ_y_min = i; }
			if(Q[2*i+0] > Q[2*iQ_x_max+0]){ iQ_x_max = i; }
			if(Q[2*i+1] > Q[2*iQ_y_max+1]){ iQ_y_max = i; }
		}
		if(
			( Q[2*iQ_x_min+0] > P[2*iP_x_max+0] ) ||
			( P[2*iP_x_min+0] > Q[2*iQ_x_max+0] ) ||
			( Q[2*iQ_y_min+1] > P[2*iP_y_max+1] ) ||
			( P[2*iP_y_min+1] > Q[2*iQ_y_max+1] )
		){ return 0; }
	}

	for(iter = 0; iter <= 2*(m+n); ++iter){
/*fprintf(stderr, "iter %d, ip = %d, iq = %d, inside = %c\n", iter, ip, iq, inside); */
		double xp[2];
		if(segments_intersect(&P[2*ipp],&P[2*ip],&Q[2*iqp],&Q[2*iq],xp, NULL)){
/*fprintf(stderr, " xsect! %f,%f %f,%f %f,%f %f,%f\n", P[2*ipp+0],P[2*ipp+1],P[2*ip+0],P[2*ip+1],Q[2*iqp+0],Q[2*iqp+1],Q[2*iq+0],Q[2*iq+1]); */
			if(first_xsected && first_iter != iter-1){ /* if the first intersection was NOT found during the previous iteration */
				if(ip == ipf && iq == iqf){ break; } /* if this intersection is the same as the first intersection */
			}
			if(*ni >= nPi){ Pi_full = 1; }
			if(!Pi_full){ Pi[2*(*ni)+0] = xp[0]; Pi[2*(*ni)+1] = xp[1]; (*ni)++; }
/*fprintf(stderr, "  Adding %f,%f\n", Pi[2*((*ni)-1)+0], Pi[2*((*ni)-1)+1]); */
			if(geom_orient2d(&Q[2*iqp],&Q[2*iq],&P[2*ip]) >= 0){
				inside = 'P';
			}else{ inside = 'Q'; }
			
			if(!first_xsected){
				first_xsected = 1;
				ipf = ip; iqf = iq;
				first_iter = iter;
			}
		}
		xp[0] = P[2*ip+0] + (Q[2*iq+0] - P[2*ipp+0]);
		xp[1] = P[2*ip+1] + (Q[2*iq+1] - P[2*ipp+1]);
		if(geom_orient2d(&Q[2*iqp],&Q[2*iq],xp)/*Cross(Q[2*iq]-Q[2*iqp],P[2*ip]-P[2*ipp])*/ >= 0){
			if(geom_orient2d(&Q[2*iqp],&Q[2*iq],&P[2*ip]) >= 0){ /* advance Q */
				if(inside == 'Q'){
					if(*ni >= nPi){ Pi_full = 1; }
					if(!Pi_full){ Pi[2*(*ni)+0] = Q[2*iq+0]; Pi[2*(*ni)+1] = Q[2*iq+1]; (*ni)++; }
				}
				iqp = iq;
				iq = (iq+1)%m;
			}else{ /* advance P */
				if(inside == 'P'){
					if(*ni >= nPi){ Pi_full = 1; }
					if(!Pi_full){ Pi[2*(*ni)+0] = P[2*ip+0]; Pi[2*(*ni)+1] = P[2*ip+1]; (*ni)++; }
				}
				ipp = ip;
				ip = (ip+1)%n;
			}
		}else{
			if(geom_orient2d(&P[2*ipp],&P[2*ip],&Q[2*iq]) >= 0){ /* advance P */
				if(inside == 'P'){
					if(*ni >= nPi){ Pi_full = 1; }
					if(!Pi_full){ Pi[2*(*ni)+0] = P[2*ip+0]; Pi[2*(*ni)+1] = P[2*ip+1]; (*ni)++; }
				}
				ipp = ip;
				ip = (ip+1)%n;
			}else{ /* advance Q */
				if(inside == 'Q'){
					if(*ni >= nPi){ Pi_full = 1; }
					if(!Pi_full){ Pi[2*(*ni)+0] = Q[2*iq+0]; Pi[2*(*ni)+1] = Q[2*iq+1]; (*ni)++; }
				}
				iqp = iq;
				iq = (iq+1)%m;
			}
		}
	}
	/* At this point, either P in Q, Q in P, or they don't intersect */
	if(*ni == 0){
		int flag = 1;
		for(j = 0; j < n; ++j){ /* really we only need to check j == 0, but due to degeneracy, it is safest to check all */
			for(i = 0; i < m; ++i){
				if(geom_orient2d(&Q[2*i],&Q[2*((i+1)%m)], &P[2*j]) < 0){
					flag = 0; j = n+1; break;
				}
			}
		}
		if(flag){ /* P in Q */
			if(*ni+n >= nPi){ Pi_full = 1; }
			if(!Pi_full){
				for(i = 0; i < n; ++i){
					Pi[2*(*ni)+0] = P[2*i+0]; Pi[2*(*ni)+1] = P[2*i+1]; (*ni)++;
				}
				return 1;
			}
		}else{
			flag = 1;
			for(j = 0; j < m; ++j){ /* really we only need to check j == 0, but due to degeneracy, it is safest to check all */
				for(i = 0; i < n; ++i){
					if(geom_orient2d(&P[2*i],&P[2*((i+1)%n)],&Q[2*j]) < 0){
						flag = 0; j = m+1; break;
					}
				}
			}
			if(flag){ /* Q in P */
				if(*ni+m >= nPi){ Pi_full = 1; }
				if(!Pi_full){
					for(i = 0; i < m; ++i){
						Pi[2*(*ni)+0] = Q[2*i+0]; Pi[2*(*ni)+1] = Q[2*i+1]; (*ni)++;
					}
					return 2;
				}
			}
		}
	}
	if(Pi_full){
		return -10;
	}else{
		return 0;
	}
}
