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
			{
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
			}
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

int geom_shape2d_contains(const geom_shape2d *s, const double p[2]){
	// offset vector from local origin
	const double po[2] = {p[0]-s->org[0], p[1]-s->org[1]};
	switch(s->type){
	case GEOM_SHAPE2D_ELLIPSE:
		{
			double v[2];
			geom_matvec2d(s->s.ellipse.B, po, v);
			return geom_norm2d(v) <= 1.;
		}
	case GEOM_SHAPE2D_POLYGON:
		return geom_polygon_inside2d(s->s.polygon.nv, s->s.polygon.v, po);
	default:
		return 0;
	}
}

int geom_shape3d_contains(const geom_shape3d *s, const double p[3]){
	const double po[3] = {p[0]-s->org[0], p[1]-s->org[1], p[2]-s->org[2]};
	switch(s->type){
	case GEOM_SHAPE3D_BLOCK:
		{
			unsigned i, j;
			for(i = 0; i < 3; ++i){
				double v = 0;
				for(j = 0; j < 3; ++j){
					v += s->s.block.B[i+j*3]*po[j];
				}
				if(fabs(v) > 1){ return 0; }
			}
			return 1;
		}
	case GEOM_SHAPE3D_POLY:
		return geom_convex_inside3d(s->s.poly.np, s->s.poly.p, po);
	case GEOM_SHAPE3D_ELLIPSOID:
		{
			double v[3];
			geom_matvec3d(s->s.ellipsoid.B, po, v);
			return geom_norm3d(v) <= 1.;
		}
	case GEOM_SHAPE3D_FRUSTUM:
		{
			double xyz[3];
			geom_matTvec3d(s->s.frustum.Q, po, xyz);
			xyz[0] /= s->s.frustum.r_base;
			xyz[1] /= s->s.frustum.r_base;
			xyz[2] /= s->s.frustum.len;
			if(xyz[2] < 0. || xyz[2] > 1.){ return 0; }
			return geom_norm2d(xyz) <= (1. + (s->s.frustum.r_tip/s->s.frustum.r_base-1.) * xyz[2]);
		}
	case GEOM_SHAPE3D_EXTRUSION:
		{
			double xyz[3];
			geom_matTvec3d(s->s.frustum.Q, po, xyz);
			if(xyz[2] < 0. || xyz[2] > s->s.frustum.len){ return 0; }
			return geom_shape2d_contains(&s->s.extrusion.s2, xyz);
		}
	default:
		return 0;
	}
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
	case GEOM_SHAPE3D_POLY:
		geom_convex_normal3d(s->s.poly.np, s->s.poly.p, p, n);
		break;
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

