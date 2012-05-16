#include <geom_predicates.h>
#include <geom_la.h>
#include <geom_poly.h>
#include <geom_shapes.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>

#define FFMT "%.14g"

int geom_shape2d_output_postscript(const geom_shape2d *s2d, FILE *fp, int opts){
	if(opts & GEOM_SHAPE2D_OUTPUT_POSTSCRIPT_NOTRANSLATE){
	}else{
		fprintf(fp, "gsave " FFMT " " FFMT " translate\n", s2d->org[0], s2d->org[1]);
	}
	switch(s2d->type){
	case GEOM_SHAPE2D_ELLIPSE:
		{
			double m[4] = {
				s2d->s.ellipse.A[0],
				s2d->s.ellipse.A[1],
				s2d->s.ellipse.A[2],
				s2d->s.ellipse.A[3]
			};
			double u[4], s[2], vt[4];
			geom_matsvd2d(m, u, s, vt);
			fprintf(fp, "matrix currentmatrix [" FFMT " " FFMT " " FFMT " " FFMT " 0 0] concat\n",
				u[0]*s[0], u[1]*s[0],
				u[2]*s[1], u[3]*s[1]
			);
			fprintf(fp, "0 0 1 0 360 arc setmatrix stroke\n");
		}
		break;
	case GEOM_SHAPE2D_POLYGON:
		{
			int i;
			fprintf(fp, FFMT " " FFMT " moveto ", s2d->s.polygon.v[2*0+0], s2d->s.polygon.v[2*0+1]);
			for(i = 1; i < s2d->s.polygon.nv; ++i){
				fprintf(fp, FFMT " " FFMT " lineto ", s2d->s.polygon.v[2*i+0], s2d->s.polygon.v[2*i+1]);
			}
			fprintf(fp, "closepath\n");
		}
		break;
	default:
		break;
	}
	if((GEOM_SHAPE2D_OUTPUT_POSTSCRIPT_FILL & opts) && !(GEOM_SHAPE2D_OUTPUT_POSTSCRIPT_NOSTROKE & opts)){
		fprintf(fp, "gsave fill grestore\n");
	}else if(GEOM_SHAPE2D_OUTPUT_POSTSCRIPT_FILL & opts){
		fprintf(fp, "fill\n");
	}else if(!(GEOM_SHAPE2D_OUTPUT_POSTSCRIPT_NOSTROKE & opts)){
		fprintf(fp, "stroke\n");
	}
	
	if(opts & GEOM_SHAPE2D_OUTPUT_POSTSCRIPT_NOTRANSLATE){
	}else{
		fprintf(fp, "grestore\n");
	}

	return 0;
}

int geom_aabb2d_output_postscript(const geom_aabb2d *s, FILE *fp, int opts){
	if(opts & GEOM_SHAPE2D_OUTPUT_POSTSCRIPT_NOTRANSLATE){
	}else{
		fprintf(fp, "gsave " FFMT " " FFMT " translate\n", s->c[0], s->c[1]);
	}
	fprintf(fp, FFMT " " FFMT " moveto "
		FFMT " " FFMT " lineto "
		FFMT " " FFMT " lineto "
		FFMT " " FFMT " lineto closepath\n",
		-s->h[0], -s->h[1],
		+s->h[0], -s->h[1],
		+s->h[0], +s->h[1],
		-s->h[0], +s->h[1]);
		
	if((GEOM_SHAPE2D_OUTPUT_POSTSCRIPT_FILL & opts) && !(GEOM_SHAPE2D_OUTPUT_POSTSCRIPT_NOSTROKE & opts)){
		fprintf(fp, "gsave fill grestore\n");
	}else if(GEOM_SHAPE2D_OUTPUT_POSTSCRIPT_FILL & opts){
		fprintf(fp, "fill\n");
	}else if(!(GEOM_SHAPE2D_OUTPUT_POSTSCRIPT_NOSTROKE & opts)){
		fprintf(fp, "stroke\n");
	}
	
	if(opts & GEOM_SHAPE2D_OUTPUT_POSTSCRIPT_NOTRANSLATE){
	}else{
		fprintf(fp, "grestore\n");
	}

	return 0;
}


void geom_aabb3d_union(geom_aabb3d *b1, const geom_aabb3d *b2){
	int i;
	for(i = 0; i < 3; ++i){
		/*
		double mn1 = b1->c[i] - b1->h[i];
		double mn2 = b2->c[i] - b2->h[i];
		double mx1 = b1->c[i] + b1->h[i];
		double mx2 = b2->c[i] + b2->h[i];
		if(mn2 < mn1){ mn1 = mn2; }
		if(mx2 > mx1){ mx1 = mx2; }
		b1->c[i] = 0.5*(mx1+mn1);
		b1->h[i] = 0.5*(mx1-mn1);
		*/
		// To prevent under/overflow and rounding errors, we need to do this more carefully.
		// We can factor out the common center:
		//   nc = new center nc, nh = new halfwidth
		//   nc = 0.5*max(c1+h1,c2+h2) + 0.5*min(c1-h1,c2-h2)
		//   nc = 0.5*c1 + 0.5*c2 + 0.5*max(c1-c2+h1,c2-c1+h2) + 0.5*min(c1-c2-h1,c2-c1-h2)
		//   nh = 0.5*max(c1+h1,c2+h2) - 0.5*min(c1-h1,c2-h2)
		//        analogous decomposition
		double cp = 0.5 * b1->c[i] + 0.5 * b2->c[i];
		double cm = 0.5 * b1->c[i] - 0.5 * b2->c[i];
		double mx1 = cm + 0.5 * b1->h[i];
		double mx2 = 0.5 * b2->h[i] - cm;
		double mn1 = cm - 0.5 * b1->h[i];
		double mn2 = -0.5 * b2->h[i] - cm;
		if(mx2 > mx1){ mx1 = mx2; }
		if(mn2 < mn1){ mn1 = mn2; }
		b1->c[i] = cp + mx1 + mn1;
		b1->h[i] = mx1 - mn1;
	}
}
void geom_aabb2d_union(geom_aabb2d *b1, const geom_aabb2d *b2){
	int i;
	for(i = 0; i < 2; ++i){
		double mn1 = b1->c[i] - b1->h[i];
		double mn2 = b2->c[i] - b2->h[i];
		double mx1 = b1->c[i] + b1->h[i];
		double mx2 = b2->c[i] + b2->h[i];
		if(mn2 < mn1){ mn1 = mn2; }
		if(mx2 > mx1){ mx1 = mx2; }
		b1->c[i] = 0.5*(mx1+mn1);
		b1->h[i] = 0.5*(mx1-mn1);
	}
}
void geom_aabb3d_union_pt(geom_aabb3d *b1, const double p[3]){
	int i;
	for(i = 0; i < 3; ++i){
		double mn1 = b1->c[i] - b1->h[i];
		double mx1 = b1->c[i] + b1->h[i];
		if(p[i] < mn1){ mn1 = p[i]; }
		if(p[i] > mx1){ mx1 = p[i]; }
		b1->c[i] = 0.5*(mx1+mn1);
		b1->h[i] = 0.5*(mx1-mn1);
	}
}
void geom_aabb2d_union_pt(geom_aabb2d *b1, const double p[2]){
	int i;
	for(i = 0; i < 2; ++i){
		double mn1 = b1->c[i] - b1->h[i];
		double mx1 = b1->c[i] + b1->h[i];
		if(p[i] < mn1){ mn1 = p[i]; }
		if(p[i] > mx1){ mx1 = p[i]; }
		b1->c[i] = 0.5*(mx1+mn1);
		b1->h[i] = 0.5*(mx1-mn1);
	}
}

int geom_aabb3d_contains(const geom_aabb3d *b, const double p[3]){
	return
		(fabs(p[0]-b->c[0]) < b->h[0]) &&
		(fabs(p[1]-b->c[1]) < b->h[1]) &&
		(fabs(p[2]-b->c[2]) < b->h[2]);
}
int geom_aabb2d_contains(const geom_aabb2d *b, const double p[2]){
	return
		(fabs(p[0]-b->c[0]) < b->h[0]) &&
		(fabs(p[1]-b->c[1]) < b->h[1]);
}

int geom_aabb3d_intersects(const geom_aabb3d *a, const geom_aabb3d *b){
	// separating axis test
	return
		(fabs(b->c[0] - a->c[0]) <= (a->h[0] + b->h[0])) &&
		(fabs(b->c[1] - a->c[1]) <= (a->h[1] + b->h[1])) &&
		(fabs(b->c[2] - a->c[2]) <= (a->h[2] + b->h[2]));
}

int geom_aabb2d_intersects(const geom_aabb2d *a, const geom_aabb2d *b){
	// separating axis test
	return
		(fabs(b->c[0] - a->c[0]) <= (a->h[0] + b->h[0])) &&
		(fabs(b->c[1] - a->c[1]) <= (a->h[1] + b->h[1]));
}

int geom_shape2d_get_aabb(const geom_shape2d *s2d, geom_aabb2d *b){
	switch(s2d->type){
	case GEOM_SHAPE2D_ELLIPSE:
		{
			// The ellipse is defined by the equation
			//   {x,y}.(A.{x,y}) == 1
			// Expanding and setting the descriminant equal to zero gives the following simple solution.
			double At[4] = {
				s2d->s.ellipse.A[0],
				s2d->s.ellipse.A[2],
				s2d->s.ellipse.A[1],
				s2d->s.ellipse.A[3]
			};
			b->c[0] = s2d->org[0]; b->c[1] = s2d->org[1];
			b->h[0] = geom_norm2d(&At[0]);
			b->h[1] = geom_norm2d(&At[2]);
			return 0;
		}
	case GEOM_SHAPE2D_POLYGON:
		{
			int i;
			double mn[2] = {s2d->s.polygon.v[0],s2d->s.polygon.v[1]};
			double mx[2] = {s2d->s.polygon.v[0],s2d->s.polygon.v[1]};
			for(i = 1; i < s2d->s.polygon.nv; ++i){
				if(s2d->s.polygon.v[2*i+0] < mn[0]){
					mn[0] = s2d->s.polygon.v[2*i+0];
				}
				if(s2d->s.polygon.v[2*i+1] < mn[1]){
					mn[1] = s2d->s.polygon.v[2*i+1];
				}
				if(s2d->s.polygon.v[2*i+0] > mx[0]){
					mx[0] = s2d->s.polygon.v[2*i+0];
				}
				if(s2d->s.polygon.v[2*i+1] > mx[1]){
					mx[1] = s2d->s.polygon.v[2*i+1];
				}
			}
			b->c[0] = s2d->org[0] + 0.5*mn[0] + 0.5*mx[0];
			b->c[1] = s2d->org[1] + 0.5*mn[1] + 0.5*mx[1];
			b->h[0] = 0.5*(mx[0] - mn[0]);
			b->h[1] = 0.5*(mx[1] - mn[1]);
			return 0;
		}
	default:
		return -1;
	}
}

int geom_shape3d_get_aabb(const geom_shape3d *s, geom_aabb3d *b){
	switch(s->type){
	case GEOM_SHAPE3D_TET:
		{
			unsigned i, j;
			double mn[3] = { s->s.tet.v[0], s->s.tet.v[1], s->s.tet.v[2] };
			double mx[3] = { s->s.tet.v[0], s->s.tet.v[1], s->s.tet.v[2] };
			for(i = 1; i < 4; ++i){
				for(j = 0; j < 3; ++j){
					if(s->s.tet.v[3*i+j] < mn[j]){ mn[j] = s->s.tet.v[3*i+j]; }
					if(s->s.tet.v[3*i+j] > mx[j]){ mx[j] = s->s.tet.v[3*i+j]; }
				}
			}
			for(j = 0; j < 3; ++j){
				b->c[j] = s->org[j] + 0.5*mn[j] + 0.5*mx[j];
				b->h[j] = 0.5*(mx[j] - mn[j]);
			}
			return 0;
		}
	case GEOM_SHAPE3D_FRUSTUM:
		{
			// compute bounding boxes of end caps and union
			// The end caps are circles with radius R and normal n:
			//   (r-c).n = 0, (r-c).(r-c) = R^2, (r-c) = {x,y,z}
			// Solving for x,y,z and setting discriminant equal to zero gives
			//   x = R/|n| sqrt(ny^2+nz^2)
			//   y = R/|n| sqrt(nx^2+nz^2)
			//   z = R/|n| sqrt(nx^2+ny^2)
			// The centers of the boxes are just the centers of the circles.
			geom_aabb3d btip;
			b->c[0] = s->org[0]; b->c[1] = s->org[1]; b->c[2] = s->org[2];
			btip.c[0] = s->org[0] + s->s.frustum.len * s->s.frustum.Q[6];
			btip.c[1] = s->org[1] + s->s.frustum.len * s->s.frustum.Q[7];
			btip.c[2] = s->org[2] + s->s.frustum.len * s->s.frustum.Q[8];
			b->h[0] = hypot(s->s.frustum.Q[7], s->s.frustum.Q[8]);
			b->h[1] = hypot(s->s.frustum.Q[6], s->s.frustum.Q[8]);
			b->h[2] = hypot(s->s.frustum.Q[6], s->s.frustum.Q[7]);
			btip.h[0] = s->s.frustum.r_tip * b->h[0];
			btip.h[1] = s->s.frustum.r_tip * b->h[1];
			btip.h[2] = s->s.frustum.r_tip * b->h[2];
			b->h[0] *= s->s.frustum.r_base;
			b->h[1] *= s->s.frustum.r_base;
			b->h[2] *= s->s.frustum.r_base;
			geom_aabb3d_union(b, &btip);
			return 0;
		}
	case GEOM_SHAPE3D_ELLIPSOID:
		{
			// The ellipse is defined by the equation
			//   {x,y,z}.(P.{x,y,z}) == 1
			// The bounds are: max = sqrt(diag(inv(P)))
			// Since P = B^T B, and B = inv({e}),
			// inv(P) = inv(B)*inv(B^T) = {e}*{e}^T
			double At[9] = {
				s->s.ellipsoid.A[0],
				s->s.ellipsoid.A[3],
				s->s.ellipsoid.A[6],
				s->s.ellipsoid.A[1],
				s->s.ellipsoid.A[4],
				s->s.ellipsoid.A[7],
				s->s.ellipsoid.A[2],
				s->s.ellipsoid.A[5],
				s->s.ellipsoid.A[8]
			};
			b->c[0] = s->org[0]; b->c[1] = s->org[1]; b->c[2] = s->org[2];
			b->h[0] = geom_norm3d(&At[0]);
			b->h[1] = geom_norm3d(&At[3]);
			b->h[2] = geom_norm3d(&At[6]);
		}
		return 0;
	case GEOM_SHAPE3D_BLOCK:
		{
			static const double sgn[8] = {
				1.,1.,
				1.,-1.,
				-1.,1.,
				-1.,-1.
			};
			unsigned i;
			b->c[0] = s->org[0]; b->c[1] = s->org[1]; b->c[2] = s->org[2];
			b->h[2] = 0; b->h[1] = 0; b->h[0] = 0;
			
			for(i = 0; i < 4; ++i){
				const double u[3] = {
					fabs(s->s.block.A[3*0+0] + sgn[2*i+0]*s->s.block.A[3*1+0] + sgn[2*i+1]*s->s.block.A[3*2+0]),
					fabs(s->s.block.A[3*0+1] + sgn[2*i+0]*s->s.block.A[3*1+1] + sgn[2*i+1]*s->s.block.A[3*2+1]),
					fabs(s->s.block.A[3*0+2] + sgn[2*i+0]*s->s.block.A[3*1+2] + sgn[2*i+1]*s->s.block.A[3*2+2])
				};
				if(u[0] > b->h[0]){ b->h[0] = u[0]; }
				if(u[1] > b->h[1]){ b->h[1] = u[1]; }
				if(u[2] > b->h[2]){ b->h[2] = u[2]; }
			}
		}
		return 0;
	case GEOM_SHAPE3D_POLY:
		b->c[0] = 0;
		b->c[1] = 0;
		b->c[2] = 0;
		b->h[0] = DBL_MAX;
		b->h[1] = DBL_MAX;
		b->h[2] = DBL_MAX;
		{
			int ret = 0;
			unsigned a;
			double *wksp = (double*)malloc(sizeof(double) * (18*s->s.poly.np+21));
			for(a = 0; a < 3; ++a){
				double rmx[3], rmn[3];
				double dir[3] = {0,0,0};
				dir[a] = 1;
				if(0 == geom_convex_bound3d(s->s.poly.np, s->s.poly.p, dir, rmx, wksp)){
					dir[a] = -1;
					if(0 == geom_convex_bound3d(s->s.poly.np, s->s.poly.p, dir, rmn, wksp)){
						b->c[a] = 0.5*rmn[a] + 0.5*rmx[a];
						b->h[a] = 0.5*(rmx[a] - rmn[a]);
					}else{
						ret = 1;
					}
				}else{
					ret = 1;
				}
			}
			free(wksp);
			return ret;
		}
	default:
		return -1;
	}
}


int geom_shape2d_init(geom_shape2d *s){
	switch(s->type){
	case GEOM_SHAPE2D_ELLIPSE:
		{
			s->s.ellipse.B[0] = s->s.ellipse.A[0];
			s->s.ellipse.B[1] = s->s.ellipse.A[1];
			s->s.ellipse.B[2] = s->s.ellipse.A[2];
			s->s.ellipse.B[3] = s->s.ellipse.A[3];
			geom_matinv2d(s->s.ellipse.B);
			return 0;
		}
	default:
		return 0;
	}
}
int geom_shape3d_init(geom_shape3d *s){
	switch(s->type){
	case GEOM_SHAPE3D_TET:
		{
			double a[3] = {
				s->s.tet.v[3*1+0] - s->s.tet.v[3*0+0],
				s->s.tet.v[3*1+1] - s->s.tet.v[3*0+1],
				s->s.tet.v[3*1+2] - s->s.tet.v[3*0+2]
			};
			double b[3] = {
				s->s.tet.v[3*2+0] - s->s.tet.v[3*0+0],
				s->s.tet.v[3*2+1] - s->s.tet.v[3*0+1],
				s->s.tet.v[3*2+2] - s->s.tet.v[3*0+2]
			};
			double c[3] = {
				s->s.tet.v[3*3+0] - s->s.tet.v[3*0+0],
				s->s.tet.v[3*3+1] - s->s.tet.v[3*0+1],
				s->s.tet.v[3*3+2] - s->s.tet.v[3*0+2]
			};
			double d[3]; geom_cross3d(a, b, d);
			if(geom_dot3d(c,d) < 0){
				double t;
				t = s->s.tet.v[3*2+0]; s->s.tet.v[3*2+0] = s->s.tet.v[3*3+0]; s->s.tet.v[3*3+0] = t;
				t = s->s.tet.v[3*2+1]; s->s.tet.v[3*2+1] = s->s.tet.v[3*3+1]; s->s.tet.v[3*3+1] = t;
				t = s->s.tet.v[3*2+2]; s->s.tet.v[3*2+2] = s->s.tet.v[3*3+2]; s->s.tet.v[3*3+2] = t;
			}
			return 0;
		}
	case GEOM_SHAPE3D_BLOCK:
		{
			s->s.block.B[0] = s->s.block.A[0];
			s->s.block.B[1] = s->s.block.A[1];
			s->s.block.B[2] = s->s.block.A[2];
			s->s.block.B[3] = s->s.block.A[3];
			s->s.block.B[4] = s->s.block.A[4];
			s->s.block.B[5] = s->s.block.A[5];
			s->s.block.B[6] = s->s.block.A[6];
			s->s.block.B[7] = s->s.block.A[7];
			s->s.block.B[8] = s->s.block.A[8];
			geom_matinv3d(s->s.block.B);
			return 0;
		}
	case GEOM_SHAPE3D_ELLIPSOID:
		{
			s->s.ellipsoid.B[0] = s->s.ellipsoid.A[0];
			s->s.ellipsoid.B[1] = s->s.ellipsoid.A[1];
			s->s.ellipsoid.B[2] = s->s.ellipsoid.A[2];
			s->s.ellipsoid.B[3] = s->s.ellipsoid.A[3];
			s->s.ellipsoid.B[4] = s->s.ellipsoid.A[4];
			s->s.ellipsoid.B[5] = s->s.ellipsoid.A[5];
			s->s.ellipsoid.B[6] = s->s.ellipsoid.A[6];
			s->s.ellipsoid.B[7] = s->s.ellipsoid.A[7];
			s->s.ellipsoid.B[8] = s->s.ellipsoid.A[8];
			geom_matinv3d(s->s.ellipsoid.B);
			return 0;
		}
	case GEOM_SHAPE3D_FRUSTUM:
		{
			geom_maketriad3d(&s->s.frustum.Q[6], &s->s.frustum.Q[0], &s->s.frustum.Q[3]);
			double len = geom_norm3d(&s->s.frustum.Q[6]);
			s->s.frustum.Q[6] /= len;
			s->s.frustum.Q[7] /= len;
			s->s.frustum.Q[8] /= len;
			s->s.frustum.len *= len;
			return 0;
		}
	case GEOM_SHAPE3D_EXTRUSION:
		{
			geom_maketriad3d(&s->s.extrusion.Q[6], &s->s.extrusion.Q[0], &s->s.extrusion.Q[3]);
			double len = geom_norm3d(&s->s.extrusion.Q[6]);
			s->s.extrusion.Q[6] /= len;
			s->s.extrusion.Q[7] /= len;
			s->s.extrusion.Q[8] /= len;
			s->s.extrusion.len *= len;
			geom_shape2d_init(&s->s.extrusion.s2);
			return 0;
		}
	default:
		return 0;
	}
}

static int geom_shape2d_contains_org(const geom_shape2d *s, const double p[2]){
	switch(s->type){
	case GEOM_SHAPE2D_ELLIPSE:
		{
			double v[2];
			geom_matvec2d(s->s.ellipse.B, p, v);
			return geom_norm2d(v) <= 1.;
		}
	case GEOM_SHAPE2D_POLYGON:
		return geom_polygon_inside2d(s->s.polygon.nv, s->s.polygon.v, p);
	default:
		return 0;
	}
}

int geom_shape2d_contains(const geom_shape2d *s, const double p[2]){
	// offset vector from local origin
	const double po[2] = {p[0]-s->org[0], p[1]-s->org[1]};
	return geom_shape2d_contains_org(s, po);
}

static int geom_shape3d_contains_org(const geom_shape3d *s, const double p[3]){
	switch(s->type){
	case GEOM_SHAPE3D_TET:
		{
			return
			   (geom_orient3d(&s->s.tet.v[0], &s->s.tet.v[3], &s->s.tet.v[6], p) < 0)
			&& (geom_orient3d(&s->s.tet.v[0], &s->s.tet.v[9], &s->s.tet.v[3], p) < 0)
			&& (geom_orient3d(&s->s.tet.v[0], &s->s.tet.v[6], &s->s.tet.v[9], p) < 0)
			&& (geom_orient3d(&s->s.tet.v[3], &s->s.tet.v[9], &s->s.tet.v[6], p) < 0);
		}
	case GEOM_SHAPE3D_BLOCK:
		{
			unsigned i, j;
			for(i = 0; i < 3; ++i){
				double v = 0;
				for(j = 0; j < 3; ++j){
					v += s->s.block.B[i+j*3]*p[j];
				}
				if(fabs(v) > 1){ return 0; }
			}
			return 1;
		}
	case GEOM_SHAPE3D_POLY:
		return geom_convex_inside3d(s->s.poly.np, s->s.poly.p, p);
	case GEOM_SHAPE3D_ELLIPSOID:
		{
			double v[3];
			geom_matvec3d(s->s.ellipsoid.B, p, v);
			return geom_norm3d(v) <= 1.;
		}
	case GEOM_SHAPE3D_FRUSTUM:
		{
			double xyz[3];
			geom_matTvec3d(s->s.frustum.Q, p, xyz);
			xyz[0] /= s->s.frustum.r_base;
			xyz[1] /= s->s.frustum.r_base;
			xyz[2] /= s->s.frustum.len;
			if(xyz[2] < 0. || xyz[2] > 1.){ return 0; }
			return geom_norm2d(xyz) <= (1. + (s->s.frustum.r_tip/s->s.frustum.r_base-1.) * xyz[2]);
		}
	case GEOM_SHAPE3D_EXTRUSION:
		{
			double xyz[3];
			geom_matTvec3d(s->s.frustum.Q, p, xyz);
			if(xyz[2] < 0. || xyz[2] > s->s.frustum.len){ return 0; }
			return geom_shape2d_contains(&s->s.extrusion.s2, xyz);
		}
	default:
		return 0;
	}
}

int geom_shape3d_contains(const geom_shape3d *s, const double p[3]){
	// offset vector from local origin
	const double po[3] = {p[0]-s->org[0], p[1]-s->org[1], p[2]-s->org[2]};
	return geom_shape3d_contains_org(s, po);
}

int geom_shape2d_normal(const geom_shape2d *s, const double p[2], double n[2]){
	const double po[3] = {p[0]-s->org[0], p[1]-s->org[1]};
	switch(s->type){
	case GEOM_SHAPE2D_ELLIPSE:
		{
			// The ellipsoid is defined by
			//   (B.po)^2 = 1
			// The gradient is
			//   B^T.B.po
			double Bpo[2];
			geom_matvec2d(s->s.ellipse.B, po, Bpo);
			geom_matTvec2d(s->s.ellipse.B, Bpo, n);
		}
		break;
	case GEOM_SHAPE2D_POLYGON:
		geom_polygon_normal2d(s->s.polygon.nv, s->s.polygon.v, p, n);
		return 0;
	default:
		n[0] = 0;
		n[1] = 0;
		return 0;
	}
	geom_normalize2d(n);
	return 0;
}

int geom_shape3d_normal(const geom_shape3d *s, const double p[3], double n[3]){
	const double po[3] = {p[0]-s->org[0], p[1]-s->org[1], p[2]-s->org[2]};
	switch(s->type){
	case GEOM_SHAPE3D_TET:
		// Decompose into barycentric coords
		// Face opposite minimum vertex controls normal
		// Barycentric coords are ratios of volumes of subdividing tets to total volume
		{
			double u[4] = {
				geom_orient3d(&s->s.tet.v[3], &s->s.tet.v[9], po, &s->s.tet.v[6]),
				geom_orient3d(&s->s.tet.v[0], &s->s.tet.v[6], po, &s->s.tet.v[9]),
				geom_orient3d(&s->s.tet.v[0], &s->s.tet.v[9], po, &s->s.tet.v[3]),
				geom_orient3d(&s->s.tet.v[0], &s->s.tet.v[3], po, &s->s.tet.v[6])
			};
			unsigned i = geom_imin4d(u);
			unsigned i1 = (i+1)%4;
			unsigned i2 = (i+2)%4;
			unsigned i3 = (i+3)%4;
			double a[3] = {
				s->s.tet.v[3*i2+0] - s->s.tet.v[3*i1+0],
				s->s.tet.v[3*i2+1] - s->s.tet.v[3*i1+1],
				s->s.tet.v[3*i2+2] - s->s.tet.v[3*i1+2]
			};
			double b[3] = {
				s->s.tet.v[3*i3+0] - s->s.tet.v[3*i1+0],
				s->s.tet.v[3*i3+1] - s->s.tet.v[3*i1+1],
				s->s.tet.v[3*i3+2] - s->s.tet.v[3*i1+2]
			};
			geom_cross3d(a, b, n);
			if(0 != i%2){
				n[0] = -n[0];
				n[1] = -n[1];
				n[2] = -n[2];
			}
		}
		break;
	case GEOM_SHAPE3D_POLY:
		geom_convex_normal3d(s->s.poly.np, s->s.poly.p, p, n);
		return 0;
	case GEOM_SHAPE3D_BLOCK:
		{
			// The cube is defined by
			//   norm(B.po, infty) = 1
			// where the columns of B are s->block.e. Let Bpo = B.po.
			//   norm(Bpo,infty) = max(abs(Bpo[0]),abs(Bpo[1]),abs(Bpo[2]))
			// We compute the gradient:
			//   Let I = argmax(abs(Bpo[i]))
			// Then
			//   grad = sgn(Bpo[I]) {I-th row of B}
			double Bpo[3];
			geom_matvec3d(s->s.block.B, po, Bpo);
			const double ab[3] = {fabs(Bpo[0]),fabs(Bpo[1]),fabs(Bpo[2])};
			int I = 0;
			if(ab[0] > ab[1]){
				if(ab[2] > ab[0]){
					I = 2;
				}else{
					// I = 0;
				}
			}else{
				if(ab[2] > ab[1]){
					I = 2;
				}else{
					I = 1;
				}
			}
			double sgn = (Bpo[I] > 0 ? 1. : -1);
			n[0] = sgn * s->s.block.B[I+0];
			n[1] = sgn * s->s.block.B[I+3];
			n[2] = sgn * s->s.block.B[I+6];
		}
		break;
	case GEOM_SHAPE3D_ELLIPSOID:
		{
			// The ellipsoid is defined by
			//   (B.po)^2 = 1
			// The gradient is
			//   B^T.B.po
			double Bpo[3];
			geom_matvec3d(s->s.ellipsoid.B, po, Bpo);
			geom_matTvec3d(s->s.ellipsoid.B, Bpo, n);
		}
		break;
	case GEOM_SHAPE3D_FRUSTUM:
		// Let Q = {v1, v2, axis} be an orthogonal set.
		// Let the un-transformed point r = inv(Q).po;
		// If Q = U.diag{S} where U is orthogonal, then inv(Q) = inv(diag(S)).U^T
		// The frustum is defined by
		//   max(
		//     2*norm_inf(r.z - 0.5),
		//     norm_2(r.xy) / (1 + (r_tip-1)*r.z)
		//   ) <= 1
		// Let denom = (1 + (r_tip-1)*r.z)
		// To compute the gradient, determine which of the two arguments
		// to max is active, then,
		//   grad = Q.{
		//     sgn(r.z - 0.5) * <0,0,1>
		//     r.xy / (|r.xy| denom) + (1-r_tip) norm_2(r.xy) / denom * <0,0,1>
		//   }
		{
			double r[3]; geom_matTvec3d(s->s.frustum.Q, po, r);
			r[0] /= s->s.frustum.r_base;
			r[1] /= s->s.frustum.r_base;
			r[2] /= s->s.frustum.len;
			const double nrxy = geom_norm2d(r);
			if(2.*fabs(r[2]-0.5) > 1.){
				double sgn = (r[2] > 0.5 ? 1. : -1.);
				n[0] = sgn*s->s.frustum.Q[6];
				n[1] = sgn*s->s.frustum.Q[7];
				n[2] = sgn*s->s.frustum.Q[8];
			}else{
				// multiply through by denom:
				//   r.xy / |r.xy| + (1-r_tip) |r.xy| * <0,0,1>
				r[0] /= nrxy;
				r[1] /= nrxy;
				r[2] = (1.-s->s.frustum.r_tip) * nrxy;
				geom_matvec3d(s->s.frustum.Q, r, n);
			}
		}
		break;
	case GEOM_SHAPE3D_EXTRUSION:
		{
			double xyz[3];
			geom_matTvec3d(s->s.frustum.Q, po, xyz);
			if(xyz[2] < 0.){
				n[0] = -s->s.frustum.Q[6];
				n[1] = -s->s.frustum.Q[7];
				n[2] = -s->s.frustum.Q[8];
			}else if(xyz[2] > s->s.frustum.len){
				n[0] = s->s.frustum.Q[6];
				n[1] = s->s.frustum.Q[7];
				n[2] = s->s.frustum.Q[8];
			}else{
				double n2[2];
				geom_shape2d_normal(&s->s.extrusion.s2, xyz, n2);
				n[0] = s->s.frustum.Q[0] * n2[0] + s->s.frustum.Q[3] * n2[1];
				n[1] = s->s.frustum.Q[1] * n2[0] + s->s.frustum.Q[4] * n2[1];
				n[2] = s->s.frustum.Q[2] * n2[0] + s->s.frustum.Q[5] * n2[1];
			}
		}
		break;
	default:
		n[0] = 0;
		n[1] = 0;
		n[2] = 0;
		return 0;
	}
	geom_normalize3d(n);
	return 0;
}

double geom_shape2d_approx_simplex_overlap(const geom_shape2d *s, const double torg[2], const double t[6], unsigned int n){
	const double org[2] = { torg[0]-s->org[0], torg[1]-s->org[1] };
	unsigned int i, j;
	unsigned int count = 0;
	double in = 1./n;
	for(i = 0; i < n; ++i){
		for(j = 0; j < n; ++j){
			double a = in*(i + (1./3.));
			double b = in*(j + (1./3.));
			double c = 1-a-b;
			double p[2] = {
				org[0] + a*t[0] + b*t[2] + c*t[4],
				org[1] + a*t[1] + b*t[3] + c*t[5]
			};
			if(geom_shape2d_contains_org(s, p)){
				count++;
			}
		}
	}
	return (double)count*2. / (double)(n*(n+1));
}

double geom_shape3d_approx_simplex_overlap(const geom_shape3d *s,const double torg[3],  const double t[12], unsigned int n){
	const double org[3] = { torg[0]-s->org[0], torg[1]-s->org[1], torg[2]-s->org[2] };
	unsigned int i, j, k;
	unsigned int count = 0;
	double in = 1./n;
	for(i = 0; i < n; ++i){
		for(j = 0; j < n-i; ++j){
			for(k = 0; k < n-i-j; ++k){
				double a = in*(i + 0.25);
				double b = in*(j + 0.25);
				double c = in*(k + 0.25);
				double d = 1-a-b-c;
				double p[3] = {
					org[0] + a*t[0] + b*t[3] + c*t[6] + d*t[ 9],
					org[1] + a*t[1] + b*t[4] + c*t[7] + d*t[10],
					org[2] + a*t[2] + b*t[5] + c*t[8] + d*t[11]
				};
				if(geom_shape3d_contains(s, p)){
					count++;
				}
			}
		}
	}
	return (double)count*6. / (double)(n*(n+1)*(n+2));
}

int geom_shape2d_intersects_simplex(const geom_shape2d *s, const double torg[2], const double t[6]){
	const double org[2] = { torg[0]-s->org[0], torg[1]-s->org[1] };
	double to[6] = {
		t[0] + org[0], t[1] + org[1],
		t[2] + org[0], t[3] + org[1],
		t[4] + org[0], t[5] + org[1]
	};
	switch(s->type){
	case GEOM_SHAPE2D_ELLIPSE:
		{
			// transform to circle
			double Bto[6];
			geom_matvec2d(s->s.ellipse.B, &to[0], &Bto[0]);
			geom_matvec2d(s->s.ellipse.B, &to[2], &Bto[2]);
			geom_matvec2d(s->s.ellipse.B, &to[4], &Bto[4]);
			unsigned i, j;
			int nout = 0;
			for(i = 0; i < 3; ++i){
				if(geom_norm2d(&Bto[2*i]) > 1.){
					nout++;
					break;
				}
			}
			// If all 3 vertices are inside, then the simplex is inside
			if(0 == nout){ return 1; }
			else if(nout < 3){ return 3; }
			// An arc could intersect an edge of the triangle twice
			// Suppose an edge has endpoints a and b, vector u=b-a.
			// Then |a + ru|^2 = 1 would have a solution with 0 < r < 1
			//    = a.a + 2ra.u + r^2 u.u = 1
			// We simply need to check the discrimiant:
			//    disc = [a.u]^2 + (u.u)(1-a.a)
			for(i = 2, j = 0; j < 3; i = j++){
				double r[2];
				//// computing vector u could have catastrophic roundoff
				//const double u[2] = { Bto[2*j+0] - Bto[2*i+0], Bto[2*j+1] - Bto[2*i+1] };
				double dt[2] = { t[2*j+0] - t[2*i+0], t[2*j+1] - t[2*i+1] };
				double u[2]; geom_matvec2d(s->s.ellipse.B, dt, u);
				double aa = Bto[2*i+0]*Bto[2*i+0] + Bto[2*i+1]*Bto[2*i+1];
				double uu = u[0]*u[0] + u[1]*u[1];
				double au = Bto[2*i+0]*u[0] + Bto[2*i+1]*u[1];
				int nsol = geom_quadraticd(uu, au, 1.-aa, r);
				if(nsol < 1){ continue; } // no intersection
				else if(nsol == 2){ // we won't worry about n == 1 tangent case
					if((0 < r[0] && r[0] < 1) || (0 < r[1] && r[1] < 1)){
						return 3;
					}
				}
			}
			// simplex could completely surround the ellipse
			const double zz[2] = {0,0};
			if(geom_polygon_inside2d(3, Bto, zz)){ return 2; }
			else{ return 0; }
		}
		break;
	case GEOM_SHAPE2D_POLYGON:
		{
			// The simplex is inside the polygon iff all 3 vertices
			// are inside and no vertex of the polygon is in the simplex.
			// The polygon is inside the simplex iff all the polygon
			// vertices are inside the simplex.
			int nv_in_tri  = 0;
			int nv_in_poly = 0;
			unsigned int i;
			for(i = 0; i < s->s.polygon.nv; ++i){
				if(geom_polygon_inside2d(3, to, &(s->s.polygon.v[2*i]))){
					nv_in_tri++;
				}
			}
			if(s->s.polygon.nv == nv_in_tri){ return 2; }
			for(i = 0; i < 3; ++i){
				if(geom_polygon_inside2d(s->s.polygon.nv, s->s.polygon.v, &to[2*i])){
					nv_in_poly++;
				}
			}
			if(3 == nv_in_poly){
				if(0 == nv_in_tri){
					return 1;
				}else{
					return 3;
				}
			}else if(0 == nv_in_poly){
				// If no vertices of the triangle are in the polygon, then either they
				// are disjoint or the polygon is inside the triangle. We already checked
				// the latter case above.
				return 0;
			}else{
				return 3;
			}
		}
	default:
		return 0;
	}
	return 0;
}
int geom_shape3d_intersects_simplex(const geom_shape3d *s, const double torg[3], const double t[12]){
	const double org[3] = { torg[0]-s->org[0], torg[1]-s->org[1], torg[2]-s->org[2] };
	double to[12] = {
		t[0] + org[0], t[1] + org[1], t[2] + org[2],
		t[3] + org[0], t[4] + org[1], t[5] + org[2],
		t[6] + org[0], t[7] + org[1], t[8] + org[2],
		t[9] + org[0], t[10] + org[1], t[11] + org[2]
	};
	switch(s->type){
	case GEOM_SHAPE3D_TET:
		{
		}
	case GEOM_SHAPE3D_POLY:
		{
		}
	case GEOM_SHAPE3D_BLOCK:
		{
		}
	case GEOM_SHAPE3D_ELLIPSOID:
		{
			// transform to sphere
			double Bto[12];
			geom_matvec3d(s->s.ellipsoid.B, &to[0], &Bto[0]);
			geom_matvec3d(s->s.ellipsoid.B, &to[3], &Bto[3]);
			geom_matvec3d(s->s.ellipsoid.B, &to[6], &Bto[6]);
			geom_matvec3d(s->s.ellipsoid.B, &to[9], &Bto[9]);
			unsigned i, j;
			int nout = 0;
			for(i = 0; i < 4; ++i){
				if(geom_norm3d(&Bto[3*i]) > 1.){
					nout++;
					break;
				}
			}
			// If all 4 vertices are inside, then the tet is inside
			if(0 == nout){ return 1; }
			else if(nout < 3){ return 3; }
			
			// Assuming the simplex is small relative to the ellipsoid, ...
			return 0;
			
			// The ellipsoid could intersect a face of the tet.
			// First check if the plane of the face goes near the origin
			// If yes, then reduce to a 2D circle triangle problem in
			// the plane.
			static const unsigned fi[4][3] = {
				{1,2,3},
				{0,3,2},
				{0,1,3},
				{0,2,1}
			};
			for(i = 0; i < 4; ++i){ // i is the face index
			}
		}
	case GEOM_SHAPE3D_FRUSTUM:
		{
		}
	case GEOM_SHAPE3D_EXTRUSION:
		{
		}
	default:
		return 0;
	}
	return 0;
}
