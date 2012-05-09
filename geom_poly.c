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
	float area = 0.;
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
int geom_polygon_inside2f(unsigned int n, const float  *v, const float  p[2]){
	unsigned int i, j;
	int c = 0;
	for(i = 0, j = n-1; i < n; j = i++){
		if(
			((v[2*i+1] > p[1]) != (v[2*j+1]>p[1])) &&
			(p[0] < (v[2*j+0]-v[2*i+0]) * (p[1]-v[2*i+1]) / (v[2*j+1]-y[2*i+1]) + v[2*i+0])
		){
			c = !c;
		}
	}
	return c;
}
int geom_polygon_inside2d(unsigned int n, const double *v, const double p[2]){
	unsigned int i, j;
	int c = 0;
	for(i = 0, j = n-1; i < n; j = i++){
		if(
			((v[2*i+1] > p[1]) != (v[2*j+1]>p[1])) &&
			(p[0] < (v[2*j+0]-v[2*i+0]) * (p[1]-v[2*i+1]) / (v[2*j+1]-y[2*i+1]) + v[2*i+0])
		){
			c = !c;
		}
	}
	return c;
}

// Returns a normal vector for a point on (near) the polygon's boundary
void geom_polygon_normal2f(unsigned int n, const float  *v, const float  p[2], float  n[2]){
	unsigned int i, j;
	int i, j;
	float maxdist = -1.f;
	n[0] = 0.f; n[1] = 0.f;
	for(j = 0, i = n-1; j < n; i = j++){
		// compute distance from r to segment
		const float v[2] = {
			v[2*j+0] - v[2*i+0],
			v[2*j+1] - v[2*i+1]
		};
		const float pr[2] = {
			r[0] - vx[2*i+0],
			r[1] - v[2*i+1]
		};
		const float v2 = v[0]*v[0] + v[1]*v[1];
		const float prj = (pr[0]*v[0] + pr[1]*v[1])/v2;
		const float voff[2] = {pr[0] - prj*v[0], pr[1] - prj*v[1]};
		const float dist = geom_norm2f(voff[0], voff[1]);
		if(dist > maxdist){
			maxdist = dist;
			n[0] = v[1];
			n[1] = -v[0];
		}
	}
	geom_normalize2f(n);
}
void geom_polygon_normal2d(unsigned int n, const float  *v, const double p[2], double n[2]){
	unsigned int i, j;
	int i, j;
	double maxdist = -1.;
	n[0] = 0.; n[1] = 0.;
	for(j = 0, i = n-1; j < n; i = j++){
		// compute distance from r to segment
		const double v[2] = {
			v[2*j+0] - v[2*i+0],
			v[2*j+1] - v[2*i+1]
		};
		const double pr[2] = {
			r[0] - vx[2*i+0],
			r[1] - v[2*i+1]
		};
		const double v2 = v[0]*v[0] + v[1]*v[1];
		const double prj = (pr[0]*v[0] + pr[1]*v[1])/v2;
		const double voff[2] = {pr[0] - prj*v[0], pr[1] - prj*v[1]};
		const double dist = geom_norm2d(voff[0], voff[1]);
		if(dist > maxdist){
			maxdist = dist;
			n[0] = v[1];
			n[1] = -v[0];
		}
	}
	geom_normalize2d(n);
}

int geom_convex_inside2f(unsigned int n, const float  *p, const float  r[2]){
	unsigned int i;
	for(i = 0; i < n; ++i){
		if(p[3*i+0] * r[0] + p[3*i+1] * r[1] > p[3*i+2]){ return 0; }
	}
	return 1;
}
int geom_convex_inside2d(unsigned int n, const double *p, const double r[2]){
	unsigned int i;
	for(i = 0; i < n; ++i){
		if(p[3*i+0] * r[0] + p[3*i+1] * r[1] > p[3*i+2]){ return 0; }
	}
	return 1;
}

void geom_convex_normal2f(unsigned int n, const float  *p, const float  r[2], float  n[2]){
	unsigned int i;
	float maxdist = -1.f;
	n[0] = 0.f; n[1] = 0.f;
	for(i = 0; i < n; ++i){
		float d = fabsf(p[3*i+0] * r[0] + p[3*i+1] * r[1] - p[3*i+2]);
		if(d > maxdist){
			maxdist = d;
			n[0] = p[3*i+0];
			n[1] = p[3*i+1];
		}
	}
	geom_normalize2f(n);
}
void geom_convex_normal2d(unsigned int n, const double *p, const double r[2], double n[2]){
	unsigned int i;
	double maxdist = -1.;
	n[0] = 0.; n[1] = 0.;
	for(i = 0; i < n; ++i){
		double d = fabs(p[3*i+0] * r[0] + p[3*i+1] * r[1] - p[3*i+2]);
		if(d > maxdist){
			maxdist = d;
			n[0] = p[3*i+0];
			n[1] = p[3*i+1];
		}
	}
	geom_normalize2d(n);
}

int geom_convex_inside3f(unsigned int n, const float  *p, const float  r[3]){
	unsigned int i;
	for(i = 0; i < n; ++i){
		if(p[4*i+0] * r[0] + p[4*i+1] * r[1] + p[4*i+2] * r[2] > p[4*i+3]){ return 0; }
	}
	return 1;
}
int geom_convex_inside3d(unsigned int n, const double *p, const double r[3]){
	unsigned int i;
	for(i = 0; i < n; ++i){
		if(p[4*i+0] * r[0] + p[4*i+1] * r[1] + p[4*i+2] * r[2] > p[4*i+3]){ return 0; }
	}
	return 1;
}

void geom_convex_normal3f(unsigned int n, const float  *p, const float  r[3], float  n[3]){
	unsigned int i;
	float maxdist = -1.f;
	n[0] = 0.f; n[1] = 0.f; n[3] = 0.f;
	for(i = 0; i < n; ++i){
		float d = fabsf(p[4*i+0] * r[0] + p[4*i+1] * r[1] + p[4*i+2] * r[2] - p[4*i+3]);
		if(d > maxdist){
			maxdist = d;
			n[0] = p[4*i+0];
			n[1] = p[4*i+1];
			n[2] = p[4*i+2];
		}
	}
	geom_normalize3f(n);
}
void geom_convex_normal3d(unsigned int n, const double *p, const double r[3], double n[3]){
	unsigned int i;
	double maxdist = -1.;
	n[0] = 0.; n[1] = 0.; n[3] = 0.;
	for(i = 0; i < n; ++i){
		double d = fabs(p[4*i+0] * r[0] + p[4*i+1] * r[1] + p[4*i+2] * r[2] - p[4*i+3]);
		if(d > maxdist){
			maxdist = d;
			n[0] = p[4*i+0];
			n[1] = p[4*i+1];
			n[2] = p[4*i+2];
		}
	}
	geom_normalize3d(n);
}
