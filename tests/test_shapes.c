#include <geom_la.h>
#include <geom_shapes.h>
#include <test_common.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void print_postscript_header(FILE *fp);
void print_postscript_footer(FILE *fp);
void draw_aabb3(FILE *fp, const geom_aabb3d *b);
void draw_samples(FILE *fp, const geom_shape3d *s3d, unsigned int n);

int main(int argc, char *argv[]){
	if(argc > 1){
		srand(atoi(argv[1]));
	}
	unsigned int i, j, k;
	unsigned int n = 256, m;
	double d1[4], d2[4];
	geom_shape2d *s2d;
	geom_shape3d *s3d;
	geom_aabb2d  b2d1, b2d2;
	geom_aabb3d  b3d1, b3d2;
	
	if(0){
		printf("Testing ellipse and polygon\n");
		FILE *fp = fopen("test_shapes2d.ps", "wb");
		print_postscript_header(fp);
		
		// Test ellipse
		s2d = (geom_shape2d*)malloc(sizeof(geom_shape2d));
		s2d->type = GEOM_SHAPE2D_ELLIPSE;
		geom_randd(2, s2d->org);
		geom_randd(4, s2d->s.ellipse.A);
		s2d->s.ellipse.B[0] = s2d->s.ellipse.A[0];
		s2d->s.ellipse.B[1] = s2d->s.ellipse.A[1];
		s2d->s.ellipse.B[2] = s2d->s.ellipse.A[2];
		s2d->s.ellipse.B[3] = s2d->s.ellipse.A[3];
		geom_matinv2d(s2d->s.ellipse.B);
		
		// Draw axes
		fprintf(fp, "%g %g moveto %g %g rlineto stroke\n",
			s2d->org[0],
			s2d->org[1],
			s2d->s.ellipse.A[0],
			s2d->s.ellipse.A[1]
		);
		fprintf(fp, "%g %g moveto %g %g rlineto stroke\n",
			s2d->org[0],
			s2d->org[1],
			s2d->s.ellipse.A[2],
			s2d->s.ellipse.A[3]
		);
		
		geom_shape2d_output_postscript(s2d, fp, 0);
		geom_shape2d_get_aabb(s2d, &b2d1);
		geom_aabb2d_output_postscript(&b2d1, fp, 0);
		fprintf(fp, "0 1 0 setrgbcolor\n");
		for(i = 0; i < n; ++i){
			geom_randd(2, d1);
			if(geom_aabb2d_contains(&b2d1, d1)){
				fprintf(fp, "%g %g 0.01 0 360 arc fill\n", d1[0], d1[1]);
			}
		}
		fprintf(fp, "0 0 1 setrgbcolor\n");
		for(i = 0; i < n; ++i){
			geom_randd(2, d1);
			if(geom_shape2d_contains(s2d, d1)){
				fprintf(fp, "%g %g 0.01 0 360 arc fill\n", d1[0], d1[1]);
			}
		}
		free(s2d);
		
		// Test polygon
		m = 7;
		s2d = (geom_shape2d*)malloc(sizeof(geom_shape2d) + 2*(m-1)*sizeof(double));
		s2d->type = GEOM_SHAPE2D_POLYGON;
		s2d->s.polygon.nv = m;
		geom_randd(2, s2d->org);
		geom_randd(2*m, s2d->s.polygon.v);
		geom_shape2d_output_postscript(s2d, fp, 0);
		geom_shape2d_get_aabb(s2d, &b2d2);
		geom_aabb2d_output_postscript(&b2d2, fp, 0);

		fprintf(fp, "1 0 0 setrgbcolor\n");
		for(i = 0; i < n; ++i){
			geom_randd(2, d1);
			if(geom_shape2d_contains(s2d, d1)){
				fprintf(fp, "%g %g 0.01 0 360 arc fill\n", d1[0], d1[1]);
			}
		}
		free(s2d);
		
		geom_aabb2d_union(&b2d1, &b2d2);
		geom_aabb2d_output_postscript(&b2d1, fp, 0);
		
		print_postscript_footer(fp);
		fclose(fp);
		printf("Check file test_shapes2d.ps\n");
	}
	n = 16;
	if(0){
		printf("Testing ellipsoid\n");
		FILE *fp = fopen("test_shapes3d_ellipsoid.pvf", "wb");
		fprintf(fp, "PVF(1)\n");
		fprintf(fp, "bs = %g\n", 1./n);
		
		s3d = (geom_shape3d*)malloc(sizeof(geom_shape3d));
		s3d->type = GEOM_SHAPE3D_ELLIPSOID;
		geom_randd(3, s3d->org);
		s3d->org[0] = 0; s3d->org[1] = 0; s3d->org[2] = 0;
		geom_randd(9, s3d->s.ellipsoid.A);
		for(i = 0; i < 9; ++i){
			s3d->s.ellipsoid.B[i] = s3d->s.ellipsoid.A[i];
		}
		geom_matinv3d(s3d->s.ellipsoid.B);
		geom_shape3d_get_aabb(s3d, &b3d1);
		draw_aabb3(fp, &b3d1);
		
		draw_samples(fp, s3d, n);
		free(s3d);
	}
	if(0){
		printf("Testing block\n");
		FILE *fp = fopen("test_shapes3d_block.pvf", "wb");
		fprintf(fp, "PVF(1)\n");
		fprintf(fp, "bs = %g\n", 1./n);
		
		s3d = (geom_shape3d*)malloc(sizeof(geom_shape3d));
		s3d->type = GEOM_SHAPE3D_BLOCK;
		geom_randd(3, s3d->org);
		s3d->org[0] = 0; s3d->org[1] = 0; s3d->org[2] = 0;
		geom_randd(9, s3d->s.block.A);
		for(i = 0; i < 9; ++i){
			s3d->s.block.B[i] = s3d->s.block.A[i];
		}
		geom_matinv3d(s3d->s.block.B);
		geom_shape3d_get_aabb(s3d, &b3d1);
		draw_aabb3(fp, &b3d1);
		
		draw_samples(fp, s3d, n);
		free(s3d);
	}
	if(0){
		printf("Testing frustum\n");
		FILE *fp = fopen("test_shapes3d_frustum.pvf", "wb");
		fprintf(fp, "PVF(1)\n");
		fprintf(fp, "bs = %g\n", 1./n);
		
		s3d = (geom_shape3d*)malloc(sizeof(geom_shape3d));
		s3d->type = GEOM_SHAPE3D_FRUSTUM;
		geom_randd(3, s3d->org);
		//s3d->org[0] = 0; s3d->org[1] = 0; s3d->org[2] = 0;
		geom_randd(3, s3d->s.frustum.Q);
		geom_normalize3d(s3d->s.frustum.Q);
		geom_maketriad3d(&s3d->s.frustum.Q[0], &s3d->s.frustum.Q[3], &s3d->s.frustum.Q[6]);
		geom_randd(1, &s3d->s.frustum.len);
		geom_randd(1, &s3d->s.frustum.r_base);
		geom_randd(1, &s3d->s.frustum.r_tip);
		s3d->s.frustum.len = 0.5+fabs(s3d->s.frustum.len);
		s3d->s.frustum.r_base = 0.5+fabs(s3d->s.frustum.r_base);
		s3d->s.frustum.r_tip = 0.5+fabs(s3d->s.frustum.r_tip);
		/*
		s3d->org[0] = 0; s3d->org[1] = 0; s3d->org[2] = 0;
		s3d->s.frustum.Q[0] = 1;
		s3d->s.frustum.Q[1] = 0;
		s3d->s.frustum.Q[2] = 0;
		s3d->s.frustum.Q[3] = 0;
		s3d->s.frustum.Q[4] = 1;
		s3d->s.frustum.Q[5] = 0;
		s3d->s.frustum.Q[6] = 0;
		s3d->s.frustum.Q[7] = 0;
		s3d->s.frustum.Q[8] = 1;
		s3d->s.frustum.len = 1;
		s3d->s.frustum.r_base = 1;
		s3d->s.frustum.r_tip = 0.5;
		*/
		
		geom_shape3d_get_aabb(s3d, &b3d1);
		draw_aabb3(fp, &b3d1);
		
		draw_samples(fp, s3d, n);
		free(s3d);
	}
	if(1){
		printf("Testing polyhedron\n");
		m = 11;
		FILE *fp = fopen("test_shapes3d_poly.pvf", "wb");
		fprintf(fp, "PVF(1)\n");
		fprintf(fp, "bs = %g\n", 1./n);
		
		s3d = (geom_shape3d*)malloc(sizeof(geom_shape3d) + sizeof(double) * 4*(m-1));
		s3d->type = GEOM_SHAPE3D_POLY;
		//geom_randd(3, s3d->org);
		s3d->org[0] = 0; s3d->org[1] = 0; s3d->org[2] = 0;
		s3d->s.poly.np = m;
		geom_randd(4*m, s3d->s.poly.p);
		for(i = 0; i < m; ++i){
			geom_normalize3d(&s3d->s.poly.p[4*i+0]);
			s3d->s.poly.p[4*i+3] = fabs(s3d->s.poly.p[4*i+3]);
		}
		
		geom_shape3d_get_aabb(s3d, &b3d1);
		draw_aabb3(fp, &b3d1);
		
		draw_samples(fp, s3d, n);
		free(s3d);
	}
	return 0;
}

void print_postscript_header(FILE *fp){
	fprintf(fp, "72 72 scale\n 4 7 translate\n0.002 setlinewidth\n");
	fprintf(fp,
		"/ellipse {\n"
		" /yrad exch def\n"
		" /xrad exch def\n"
		" /savematrix matrix currentmatrix def\n"
		" xrad yrad scale\n"
		" 0 0 1 0 360 arc\n"
		" savematrix setmatrix\n"
		"} def\n"
	);
}
void print_postscript_footer(FILE *fp){
	fprintf(fp, "showpage\n");
}

void draw_aabb3(FILE *fp, const geom_aabb3d *B){
	unsigned int a, b, c;
	int s, i, j;
	double shrink = 0.5;
	double r[3];
	for(a = 0; a < 3; ++a){
		b = (a+1)%3; c = (a+2)%3;
		for(s = -1; s <= 1; s += 2){
			r[a] = B->c[a] + s*B->h[a];
			for(i = -1; i <= 1; i += 2){
				for(j = -1; j <= 1; j += 2){
					r[b] = B->c[b] + i*shrink*B->h[b];
					r[c] = B->c[c] + j*shrink*B->h[c];
					fprintf(fp, "p{%g,%g,%g}\n", r[0], r[1], r[2]);
				}
			}
			fprintf(fp, "q{-4,-3,-2,-1}\n");
		}
	}
}

void draw_samples(FILE *fp, const geom_shape3d *s3d, unsigned int n){
	unsigned int i, j, k;
	double scale = 2;
	double r[3], v[3];
	for(k = 0; k < n; ++k){
		r[2] = ((double)k+0.50)/(double)n * scale - 0.5*scale;
		for(j = 0; j < n; ++j){
			r[1] = ((double)j+0.50)/(double)n * scale - 0.5*scale;
			for(i = 0; i < n; ++i){
				r[0] = ((double)i+0.50)/(double)n * scale - 0.5*scale;
				if(geom_shape3d_contains(s3d, r)){
					fprintf(fp, "p{%g,%g,%g}\nb{-1,bs}\n", r[0], r[1], r[2]);
				}
				
				// draw normals
				geom_shape3d_normal(s3d, r, v);
				fprintf(fp, "p{%g,%g,%g}\np{%g,%g,%g}\nv{-2,-1}\n",
					r[0], r[1], r[2],
					r[0]+v[0]/n, r[1]+v[1]/n, r[2]+v[2]/n
				);
			}
		}
	}
}
