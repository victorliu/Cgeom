#include <math.h>
#include <geom_la.h>
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

// Returns 1 if inside, 0 otherwise
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

// Returns a normal vector for a point on (near) the polygon's boundary
void geom_polygon_normal2f(unsigned int nv, const float  *v, const float  p[2], float  n[2]){
	unsigned int i, j;
	float mindist = FLT_MAX;
	unsigned int idist = -1;
	n[0] = 0.f; n[1] = 0.f;
	for(j = 0, i = nv-1; j < nv; i = j++){
		// compute distance from r to segment
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
		// compute distance from r to segment
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
	n[0] = 0.f; n[1] = 0.f; n[3] = 0.f;
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
	n[0] = 0.; n[1] = 0.; n[3] = 0.;
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
