#include <geom_la.h>
#include <geom_poly.h>
#include <test_common.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(){
	unsigned int i, j, k;
	unsigned int n = 256, m;
	float  fv[256], f1[16], f2[16], f3[16], fa, fb, fc, fd;
	double dv[256], d1[16], d2[16], d3[16], da, db, dc, dd;
	
	printf("Checking polygon area\n");
	{
		// concave example
		fv[ 0] = 0.f; fv[ 1] = 0.f;
		fv[ 2] = 2.f; fv[ 3] = 0.f;
		fv[ 4] = 2.f; fv[ 5] = 1.f;
		fv[ 6] = 1.f; fv[ 7] = 1.f;
		fv[ 8] = 1.f; fv[ 9] = 2.f;
		fv[10] = 2.f; fv[11] = 2.f;
		fv[12] = 2.f; fv[13] = 3.f;
		fv[14] = 0.f; fv[15] = 3.f;
		dv[ 0] = 0.0; dv[ 1] = 0.0;
		dv[ 2] = 2.0; dv[ 3] = 0.0;
		dv[ 4] = 2.0; dv[ 5] = 1.0;
		dv[ 6] = 1.0; dv[ 7] = 1.0;
		dv[ 8] = 1.0; dv[ 9] = 2.0;
		dv[10] = 2.0; dv[11] = 2.0;
		dv[12] = 2.0; dv[13] = 3.0;
		dv[14] = 0.0; dv[15] = 3.0;
		fa = geom_polygon_area2f(8, fv);
		da = geom_polygon_area2d(8, dv);
		checkf_small("geom_polygon_area2f concave check", fa - 5.f);
		checkd_small("geom_polygon_area2d concave check", da - 5.0);
		
		// regular polygons
		for(i = 3; i <= 60; ++i){
			for(j = 0; j < i; ++j){
				double angle = 2.*M_PI*(double)j/(double)i;
				fv[2*j+0] = cosf(angle);
				fv[2*j+1] = sinf(angle);
				dv[2*j+0] = cos(angle);
				dv[2*j+1] = sin(angle);
			}
			double area = 0.5*(double)i * sin(2.*M_PI/(double)i);
			fa = geom_polygon_area2f(i, fv);
			da = geom_polygon_area2d(i, dv);
			checkf_small("geom_polygon_area2f regular polygon check", 0.5f*(fa - (float)area));
			checkd_small("geom_polygon_area2d regular polygon check", 0.5*(da - area));
		}
	}
	printf("Checking polygon inside\n");
	{
		// We still have the 60 sided polygon from above
		int fcount = 0;
		int dcount = 0;
		j = 999999;
		for(i = 0; i < j; ++i){
			geom_randf(2, f1); f1[0] *= 2.f; f1[1] *= 2.f;
			geom_randd(2, d1); d1[0] *= 2.0; d1[1] *= 2.0;
			if(geom_polygon_inside2f(60, fv, f1)){ ++fcount; }
			if(geom_polygon_inside2d(60, dv, d1)){ ++dcount; }
		}
		double area = 0.5*(double)60. * sin(2.*M_PI/(double)60.);
		printf(" This should be small: %.14g\n", 4.*(float )fcount / (float )j - area);
		printf(" This should be small: %.14g\n", 4.*(double)dcount / (double)j - area);
		
		// Generate random star
		m = 17;
		for(i = 0; i < m; ++i){
			double angle = 2.*M_PI*(double)i/(double)m;
			geom_randf(1, f1);
			geom_randd(1, d1);
			fv[2*i+0] = fabsf(f1[0])*cosf(angle);
			fv[2*i+1] = fabsf(f1[0])*sinf(angle);
			dv[2*i+0] = fabs(d1[0])*cos(angle);
			dv[2*i+1] = fabs(d1[0])*sin(angle);
		}
		// Dump file
		FILE *fpf = fopen("test_poly_insidef.dat", "wb");
		FILE *fpd = fopen("test_poly_insided.dat", "wb");
		n = 256;
		for(j = 0; j < n; ++j){
			f1[1] = ((float )j+0.5f)/(float )n - 0.5f;
			d1[1] = ((double)j+0.50)/(double)n - 0.50;
			for(i = 0; i < n; ++i){
				f1[0] = ((float )i+0.5f)/(float )n - 0.5f;
				d1[0] = ((double)i+0.50)/(double)n - 0.50;
				fprintf(fpf, "%d\t%d\t%d\n", i, j, geom_polygon_inside2f(m, fv, f1));
				fprintf(fpd, "%d\t%d\t%d\n", i, j, geom_polygon_inside2d(m, dv, d1));
			}
			fprintf(fpf, "\n");
			fprintf(fpd, "\n");
		}
		fclose(fpf);
		fclose(fpd);
		printf("Check files test_poly_insidef.dat and test_poly_insided.dat\n");
	}
	printf("Checking polygon normal\n");
	{
		FILE *fpf = fopen("test_poly_normalf.dat", "wb");
		FILE *fpd = fopen("test_poly_normald.dat", "wb");
		n = 32;
		for(j = 0; j < n; ++j){
			f1[1] = ((float )j+0.5f)/(float )n - 0.5f;
			d1[1] = ((double)j+0.50)/(double)n - 0.50;
			for(i = 0; i < n; ++i){
				f1[0] = ((float )i+0.5f)/(float )n - 0.5f;
				d1[0] = ((double)i+0.50)/(double)n - 0.50;
				geom_polygon_normal2f(m, fv, f1, f2);
				geom_polygon_normal2d(m, dv, d1, d2);
				fprintf(fpf, "%d\t%d\t%.14g\t%.14g\n", i, j, f2[0], f2[1]);
				fprintf(fpd, "%d\t%d\t%.14g\t%.14g\n", i, j, d2[0], d2[1]);
			}
			fprintf(fpf, "\n");
			fprintf(fpd, "\n");
		}
		fclose(fpf);
		fclose(fpd);
		printf("Check files test_poly_normalf.dat and test_poly_normald.dat\n");
	}
	printf("Checking polygon inside\n");
	{
		// Generate random inequalities
		m = 7;
		geom_randf(3*m, fv);
		geom_randd(3*m, dv);
		// Normalize and make feasible
		for(i = 0; i < m; ++i){
			geom_normalize2f(&fv[3*i+0]);
			geom_normalize2d(&dv[3*i+0]);
			fv[3*i+2] = fabsf(fv[3*i+2]);
			dv[3*i+2] = fabs (dv[3*i+2]);
		}
		// Dump file
		FILE *fpf = fopen("test_cvx2_insidef.dat", "wb");
		FILE *fpd = fopen("test_cvx2_insided.dat", "wb");
		n = 256;
		for(j = 0; j < n; ++j){
			f1[1] = ((float )j+0.5f)/(float )n - 0.5f;
			d1[1] = ((double)j+0.50)/(double)n - 0.50;
			for(i = 0; i < n; ++i){
				f1[0] = ((float )i+0.5f)/(float )n - 0.5f;
				d1[0] = ((double)i+0.50)/(double)n - 0.50;
				fprintf(fpf, "%d\t%d\t%d\n", i, j, geom_convex_inside2f(m, fv, f1));
				fprintf(fpd, "%d\t%d\t%d\n", i, j, geom_convex_inside2d(m, dv, d1));
			}
			fprintf(fpf, "\n");
			fprintf(fpd, "\n");
		}
		fclose(fpf);
		fclose(fpd);
		printf("Check files test_cvx2_insidef.dat and test_cvx2_insided.dat\n");
	}
	printf("Checking polygon normal\n");
	{
		FILE *fpf = fopen("test_cvx2_normalf.dat", "wb");
		FILE *fpd = fopen("test_cvx2_normald.dat", "wb");
		n = 32;
		for(j = 0; j < n; ++j){
			f1[1] = ((float )j+0.5f)/(float )n - 0.5f;
			d1[1] = ((double)j+0.50)/(double)n - 0.50;
			for(i = 0; i < n; ++i){
				f1[0] = ((float )i+0.5f)/(float )n - 0.5f;
				d1[0] = ((double)i+0.50)/(double)n - 0.50;
				geom_convex_normal2f(m, fv, f1, f2);
				geom_convex_normal2d(m, dv, d1, d2);
				fprintf(fpf, "%d\t%d\t%.14g\t%.14g\n", i, j, f2[0], f2[1]);
				fprintf(fpd, "%d\t%d\t%.14g\t%.14g\n", i, j, d2[0], d2[1]);
			}
			fprintf(fpf, "\n");
			fprintf(fpd, "\n");
		}
		fclose(fpf);
		fclose(fpd);
		printf("Check files test_cvx2_normalf.dat and test_cvx2_normald.dat\n");
	}
	printf("Checking polyhedron inside\n");
	{
		// Generate random inequalities
		m = 17;
		n = 64;
		geom_randf(4*m, fv);
		geom_randd(4*m, dv);
		// Dump file
		FILE *fpf = fopen("test_cvx3_insidef.pvf", "wb");
		FILE *fpd = fopen("test_cvx3_insided.pvf", "wb");
		fprintf(fpf, "PVF(1)\n");
		fprintf(fpd, "PVF(1)\n");
		fprintf(fpf, "ps = 0.1\nbs = %g\n", 0.5/(double)n); // plane size, ball size
		fprintf(fpd, "ps = 0.1\nbs = %g\n", 0.5/(double)n);
		// Normalize and make feasible
		for(i = 0; i < m; ++i){
			geom_normalize3f(&fv[4*i+0]);
			geom_normalize3d(&dv[4*i+0]);
			fv[4*i+3] = fabsf(fv[4*i+3]);
			dv[4*i+3] = fabs (dv[4*i+3]);
			// Output the plane
			geom_maketriad3f(&fv[4*i+0], f1, f2);
			fprintf(fpf, "p{%g,%g,%g}\n",
				fv[4*i+3]*fv[4*i+0],
				fv[4*i+3]*fv[4*i+1],
				fv[4*i+3]*fv[4*i+2]
				);
			fprintf(fpf, "p{%g+ps*%g,%g+ps*%g,%g+ps*%g}\n",
				fv[4*i+3]*fv[4*i+0], f1[0],
				fv[4*i+3]*fv[4*i+1], f1[1],
				fv[4*i+3]*fv[4*i+2], f1[2]
				);
			fprintf(fpf, "p{%g+ps*%g,%g+ps*%g,%g+ps*%g}\n",
				fv[4*i+3]*fv[4*i+0], f2[0],
				fv[4*i+3]*fv[4*i+1], f2[1],
				fv[4*i+3]*fv[4*i+2], f2[2]
				);
			fprintf(fpf, "t{-3,-2,-1}\n");
			
			geom_maketriad3d(&dv[4*i+0], d1, d2);
			fprintf(fpd, "p{%g,%g,%g}\n",
				dv[4*i+3]*dv[4*i+0],
				dv[4*i+3]*dv[4*i+1],
				dv[4*i+3]*dv[4*i+2]
				);
			fprintf(fpd, "p{%g+ps*%g,%g+ps*%g,%g+ps*%g}\n",
				dv[4*i+3]*dv[4*i+0], d1[0],
				dv[4*i+3]*dv[4*i+1], d1[1],
				dv[4*i+3]*dv[4*i+2], d1[2]
				);
			fprintf(fpd, "p{%g+ps*%g,%g+ps*%g,%g+ps*%g}\n",
				dv[4*i+3]*dv[4*i+0], d2[0],
				dv[4*i+3]*dv[4*i+1], d2[1],
				dv[4*i+3]*dv[4*i+2], d2[2]
				);
			fprintf(fpd, "t{-3,-2,-1}\n");
		}
		for(k = 0; k < n; ++k){
			f1[2] = ((float )k+0.5f)/(float )n - 0.5f;
			d1[2] = ((double)k+0.50)/(double)n - 0.50;
			for(j = 0; j < n; ++j){
				f1[1] = ((float )j+0.5f)/(float )n - 0.5f;
				d1[1] = ((double)j+0.50)/(double)n - 0.50;
				for(i = 0; i < n; ++i){
					f1[0] = ((float )i+0.5f)/(float )n - 0.5f;
					d1[0] = ((double)i+0.50)/(double)n - 0.50;
					if(geom_convex_inside3f(m, fv, f1)){
						fprintf(fpf, "p{%g,%g,%g}\nb{-1,bs}\n", f1[0], f1[1], f1[2]);
					}
					if(geom_convex_inside3d(m, dv, d1)){
						fprintf(fpd, "p{%g,%g,%g}\nb{-1,bs}\n", f1[0], f1[1], f1[2]);
					}
				}
			}
		}
		fclose(fpf);
		fclose(fpd);
		printf("Check files test_cvx3_insidef.pvf and test_cvx3_insided.pvf\n");
	}
	return 0;
}

