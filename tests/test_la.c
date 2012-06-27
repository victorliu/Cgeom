#include <geom_la.h>
#include <test_common.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(){
	unsigned int i, j;
	const unsigned int n = 256;
	float  f1[16], f2[16], f3[16], f4[16], fa, fb, fc, fd;
	double d1[16], d2[16], d3[16], d4[16], da, db, dc, dd;
	const float  fz = sinf(0.f);
	const double dz = sin (0. );
	const float  fNaN = sqrtf(-1.f);
	const double dNaN = sqrt (-1.0);
	const float fInf = 1.f / fz;
	const float dInf = 1.0 / dz;
	
	printf("Checking norms\n");
	for(i = 0; i < n; ++i){
		// Check zero
		geom_setzerof(4, f1);
		geom_setzerod(4, d1);
		fa = geom_norm2f(f1);
		da = geom_norm2d(d1);
		checkf_equal("geom_norm2f zero check", fa, 0.f);
		checkd_equal("geom_norm2d zero check", da, 0.0);
		fb = geom_norm3f(f1);
		db = geom_norm3d(d1);
		checkf_equal("geom_norm3f zero check", fb, 0.f);
		checkd_equal("geom_norm3d zero check", db, 0.0);
		fc = geom_norm4f(f1);
		dc = geom_norm4d(d1);
		checkf_equal("geom_norm4f zero check", fc, 0.f);
		checkd_equal("geom_norm4d zero check", dc, 0.0);
		
		// Check singles
		geom_setzerof(2, f1);
		geom_setzerod(2, d1);
		geom_randf(1, &fd);
		geom_randd(1, &dd);
		f1[rand()%2] = fd;
		d1[rand()%2] = dd;
		fa = geom_norm2f(f1);
		da = geom_norm2d(d1);
		checkf_equal("geom_norm2f singleton check", fabsf(fd), fa);
		checkd_equal("geom_norm2d singleton check", fabs (dd), da);
		geom_setzerof(3, f1);
		geom_setzerod(3, d1);
		geom_randf(1, &fd);
		geom_randd(1, &dd);
		f1[rand()%3] = fd;
		d1[rand()%3] = dd;
		fb = geom_norm3f(f1);
		db = geom_norm3d(d1);
		checkf_equal("geom_norm3f singleton check", fabsf(fd), fb);
		checkd_equal("geom_norm3d singleton check", fabs (dd), db);
		geom_setzerof(4, f1);
		geom_setzerod(4, d1);
		geom_randf(1, &fd);
		geom_randd(1, &dd);
		f1[rand()%4] = fd;
		d1[rand()%4] = dd;
		fc = geom_norm4f(f1);
		dc = geom_norm4d(d1);
		checkf_equal("geom_norm4f singleton check", fabsf(fd), fc);
		checkd_equal("geom_norm4d singleton check", fabs (dd), dc);
	
		// Check general
		geom_randf(4, f1);
		geom_randd(4, d1);
		fa = geom_norm2f(f1);
		da = geom_norm2d(d1);
		checkf_small("geom_norm2f general check", fa*fa - f1[0]*f1[0] - f1[1]*f1[1]);
		checkd_small("geom_norm2d general check", da*da - d1[0]*d1[0] - d1[1]*d1[1]);
		fb = geom_norm3f(f1);
		db = geom_norm3d(d1);
		checkf_small("geom_norm3f general check", fb*fb - fa*fa - f1[2]*f1[2]);
		checkd_small("geom_norm3d general check", db*db - da*da - d1[2]*d1[2]);
		fc = geom_norm4f(f1);
		dc = geom_norm4d(d1);
		checkf_small("geom_norm4f general check", fc*fc - fb*fb - f1[3]*f1[3]);
		checkd_small("geom_norm4d general check", dc*dc - db*db - d1[3]*d1[3]);
		
		// Check NaN conditions
		geom_randf(2, f1);
		geom_randd(2, d1);
		f1[rand()%2] = fNaN;
		d1[rand()%2] = dNaN;
		fa = geom_norm2f(f1);
		da = geom_norm2d(d1);
		checkf_NaN("geom_norm2f NaN check", fa);
		checkd_NaN("geom_norm2d NaN check", da);
		
		geom_randf(3, f1);
		geom_randd(3, d1);
		f1[rand()%3] = fNaN;
		d1[rand()%3] = dNaN;
		fb = geom_norm3f(f1);
		db = geom_norm3d(d1);
		checkf_NaN("geom_norm3f NaN check", fb);
		checkd_NaN("geom_norm3d NaN check", db);
		
		geom_randf(4, f1);
		geom_randd(4, d1);
		f1[rand()%4] = fNaN;
		d1[rand()%4] = dNaN;
		fc = geom_norm4f(f1);
		dc = geom_norm4d(d1);
		checkf_NaN("geom_norm4f NaN check", fc);
		checkd_NaN("geom_norm4d NaN check", dc);
	}
	printf("Checking normalization\n");
	for(i = 0; i < n; ++i){
		// Check NaN
		geom_setzerof(2, f1);
		geom_setzerod(2, d1);
		geom_normalize2f(f1);
		geom_normalize2d(d1);
		checkf_NaN("geom_normalize2f NaN check", f1[rand()%2]);
		checkd_NaN("geom_normalize2d NaN check", d1[rand()%2]);
		geom_setzerof(3, f1);
		geom_setzerod(3, d1);
		geom_normalize3f(f1);
		geom_normalize3d(d1);
		checkf_NaN("geom_normalize3f NaN check", f1[rand()%3]);
		checkd_NaN("geom_normalize3d NaN check", d1[rand()%3]);
		geom_setzerof(4, f1);
		geom_setzerod(4, d1);
		geom_normalize4f(f1);
		geom_normalize4d(d1);
		checkf_NaN("geom_normalize3f NaN check", f1[rand()%4]);
		checkd_NaN("geom_normalize3d NaN check", d1[rand()%4]);
		
		geom_randf(2, f1);
		geom_randd(2, d1);
		geom_normalize2f(f1);
		geom_normalize2d(d1);
		checkf_small("geom_normalize2f general check", geom_norm2f(f1) - 1.f);
		checkd_small("geom_normalize2d general check", geom_norm2d(d1) - 1.0);
		geom_randf(3, f1);
		geom_randd(3, d1);
		geom_normalize3f(f1);
		geom_normalize3d(d1);
		checkf_small("geom_normalize3f general check", geom_norm3f(f1) - 1.f);
		checkd_small("geom_normalize3d general check", geom_norm3d(d1) - 1.0);
		geom_randf(4, f1);
		geom_randd(4, d1);
		geom_normalize4f(f1);
		geom_normalize4d(d1);
		checkf_small("geom_normalize4f general check", geom_norm4f(f1) - 1.f);
		checkd_small("geom_normalize4d general check", geom_norm4d(d1) - 1.0);
	}
	printf("Checking cross products\n");
	for(i = 0; i < n; ++i){
		geom_randf(3, f1);
		geom_randd(3, d1);
		geom_randf(3, f2);
		geom_randd(3, d2);
		geom_cross3f(f1, f2, f3);
		geom_cross3d(d1, d2, d3);
		checkf_small("geom_cross3f orth check1", f1[0]*f3[0] + f1[1]*f3[1] + f1[2]*f3[2]);
		checkf_small("geom_cross3f orth check2", f2[0]*f3[0] + f2[1]*f3[1] + f2[2]*f3[2]);
		checkd_small("geom_cross3d orth check1", d1[0]*d3[0] + d1[1]*d3[1] + d1[2]*d3[2]);
		checkd_small("geom_cross3d orth check2", d2[0]*d3[0] + d2[1]*d3[1] + d2[2]*d3[2]);
	}
	printf("Checking maketriad\n");
	for(i = 0; i < n; ++i){
		geom_randf(3, f1);
		geom_randd(3, d1);
		geom_maketriad3f(f1, f2, f3);
		geom_maketriad3d(d1, d2, d3);
		checkf_small("geom_maketriad3f orth check1", f1[0]*f2[0] + f1[1]*f2[1] + f1[2]*f2[2]);
		checkf_small("geom_maketriad3f orth check2", f2[0]*f3[0] + f2[1]*f3[1] + f2[2]*f3[2]);
		checkf_small("geom_maketriad3f orth check3", f3[0]*f1[0] + f3[1]*f1[1] + f3[2]*f1[2]);
		checkf_small("geom_maketriad3f norm check1", geom_norm3f(f2) - 1.f);
		checkf_small("geom_maketriad3f norm check2", geom_norm3f(f3) - 1.f);
		geom_cross3f(f1, f2, f4);
		checkf_pos  ("geom_maketriad3f orient check", f3[0]*f4[0] + f3[1]*f4[1] + f3[2]*f4[2]); 
		
		checkf_small("geom_maketriad3d orth check1", d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2]);
		checkd_small("geom_maketriad3d orth check2", d2[0]*d3[0] + d2[1]*d3[1] + d2[2]*d3[2]);
		checkd_small("geom_maketriad3d orth check3", d3[0]*d1[0] + d3[1]*d1[1] + d3[2]*d1[2]);
		checkd_small("geom_maketriad3d norm check1", geom_norm3d(d2) - 1.d);
		checkd_small("geom_maketriad3d norm check2", geom_norm3d(d3) - 1.d);
		geom_cross3d(d1, d2, d4);
		checkf_pos  ("geom_maketriad3d orient check", d3[0]*d4[0] + d3[1]*d4[1] + d3[2]*d4[2]); 
	}
	printf("Checking matvec\n");
	{
		for(i = 0; i < 16; ++i){
			f1[i] = (float)(i+1);
			d1[i] = (double)(i+1);
		}
		for(i = 0; i < 4; ++i){
			f2[i] = (float)(i+17);
			d2[i] = (double)(i+17);
		}
		geom_matvec2f(f1, f2, f3);
		geom_matvec2d(d1, d2, d3);
		checkf_small("geom_matvec2f check1", f3[0] - 71.f);
		checkf_small("geom_matvec2f check2", f3[1] - 106.f);
		checkf_small("geom_matvec2d check1", d3[0] - 71.0);
		checkf_small("geom_matvec2d check2", d3[1] - 106.0);
		for(i = 0; i < 16; ++i){
			f1[i] = (float)(i+1);
			d1[i] = (double)(i+1);
		}
		for(i = 0; i < 4; ++i){
			f2[i] = (float)(i+17);
			d2[i] = (double)(i+17);
		}
		geom_matvec3f(f1, f2, f3);
		geom_matvec3d(d1, d2, d3);
		checkf_small("geom_matvec3f check1", f3[0] - 222.f);
		checkf_small("geom_matvec3f check2", f3[1] - 276.f);
		checkf_small("geom_matvec3f check3", f3[2] - 330.f);
		checkf_small("geom_matvec3d check1", d3[0] - 222.0);
		checkf_small("geom_matvec3d check2", d3[1] - 276.0);
		checkf_small("geom_matvec3d check3", d3[2] - 330.0);
		for(i = 0; i < 16; ++i){
			f1[i] = (float)(i+1);
			d1[i] = (double)(i+1);
		}
		for(i = 0; i < 4; ++i){
			f2[i] = (float)(i+17);
			d2[i] = (double)(i+17);
		}
		geom_matvec4f(f1, f2, f3);
		geom_matvec4d(d1, d2, d3);
		checkf_small("geom_matvec4f check1", f3[0] - 538.f);
		checkf_small("geom_matvec4f check2", f3[1] - 612.f);
		checkf_small("geom_matvec4f check3", f3[2] - 686.f);
		checkf_small("geom_matvec4f check4", f3[3] - 760.f);
		checkf_small("geom_matvec4d check1", d3[0] - 538.0);
		checkf_small("geom_matvec4d check2", d3[1] - 612.0);
		checkf_small("geom_matvec4d check3", d3[2] - 686.0);
		checkf_small("geom_matvec4d check4", d3[3] - 760.0);
	}
	printf("Checking matmat and matinv\n");
	for(i = 0; i < 1; ++i){
		geom_randf(16, f1);
		geom_randd(16, d1);
		for(j = 0; j < 16; ++j){
			f2[j] = f1[j];
			d2[j] = d1[j];
		}
		geom_matinv2f(f1);
		geom_matinv2d(d1);
		geom_matmat2f(f1, f2, f3);
		geom_matmat2d(d1, d2, d3);
		checkf_small("geom_matinv2f and geom_matmat2f check1", 0.125f*(f3[0] - 1.f));
		checkf_small("geom_matinv2f and geom_matmat2f check2", 0.125f*(f3[1] - 0.f));
		checkf_small("geom_matinv2f and geom_matmat2f check3", 0.125f*(f3[2] - 0.f));
		checkf_small("geom_matinv2f and geom_matmat2f check4", 0.125f*(f3[3] - 1.f));
		checkd_small("geom_matinv2d and geom_matmat2d check1", 0.125*(d3[0] - 1.));
		checkd_small("geom_matinv2d and geom_matmat2d check2", 0.125*(d3[1] - 0.));
		checkd_small("geom_matinv2d and geom_matmat2d check3", 0.125*(d3[2] - 0.));
		checkd_small("geom_matinv2d and geom_matmat2d check4", 0.125*(d3[3] - 1.));
		
		geom_randf(16, f1);
		geom_randd(16, d1);
		for(j = 0; j < 16; ++j){
			f2[j] = f1[j];
			d2[j] = d1[j];
		}
		geom_matinv3f(f1);
		geom_matinv3d(d1);
		geom_matmat3f(f1, f2, f3);
		geom_matmat3d(d1, d2, d3);
		fd = 0.1f;
		dd = 0.1;
		checkf_small("geom_matinv3f and geom_matmat3f check1", fd*(f3[0] - 1.f));
		checkf_small("geom_matinv3f and geom_matmat3f check2", fd*(f3[1] - 0.f));
		checkf_small("geom_matinv3f and geom_matmat3f check3", fd*(f3[2] - 0.f));
		checkf_small("geom_matinv3f and geom_matmat3f check4", fd*(f3[3] - 0.f));
		checkf_small("geom_matinv3f and geom_matmat3f check5", fd*(f3[4] - 1.f));
		checkf_small("geom_matinv3f and geom_matmat3f check6", fd*(f3[5] - 0.f));
		checkf_small("geom_matinv3f and geom_matmat3f check7", fd*(f3[6] - 0.f));
		checkf_small("geom_matinv3f and geom_matmat3f check8", fd*(f3[7] - 0.f));
		checkf_small("geom_matinv3f and geom_matmat3f check9", fd*(f3[8] - 1.f));
		checkd_small("geom_matinv3d and geom_matmat3d check1", dd*(d3[0] - 1.));
		checkd_small("geom_matinv3d and geom_matmat3d check2", dd*(d3[1] - 0.));
		checkd_small("geom_matinv3d and geom_matmat3d check3", dd*(d3[2] - 0.));
		checkd_small("geom_matinv3d and geom_matmat3d check4", dd*(d3[3] - 0.));
		checkd_small("geom_matinv3d and geom_matmat3d check5", dd*(d3[4] - 1.));
		checkd_small("geom_matinv3d and geom_matmat3d check6", dd*(d3[5] - 0.));
		checkd_small("geom_matinv3d and geom_matmat3d check7", dd*(d3[6] - 0.));
		checkd_small("geom_matinv3d and geom_matmat3d check8", dd*(d3[7] - 0.));
		checkd_small("geom_matinv3d and geom_matmat3d check9", dd*(d3[8] - 1.));
		
		geom_randf(16, f1);
		geom_randd(16, d1);
		for(j = 0; j < 16; ++j){
			f2[j] = f1[j];
			d2[j] = d1[j];
		}
		geom_matinv4f(f1);
		geom_matinv4d(d1);
		geom_matmat4f(f1, f2, f3);
		geom_matmat4d(d1, d2, d3);
		fd = 0.05f;
		dd = 0.05;
		checkf_small("geom_matinv4f and geom_matmat4f check1" , fd*(f3[ 0] - 1.f));
		checkf_small("geom_matinv4f and geom_matmat4f check2" , fd*(f3[ 1] - 0.f));
		checkf_small("geom_matinv4f and geom_matmat4f check3" , fd*(f3[ 2] - 0.f));
		checkf_small("geom_matinv4f and geom_matmat4f check4" , fd*(f3[ 3] - 0.f));
		checkf_small("geom_matinv4f and geom_matmat4f check5" , fd*(f3[ 4] - 0.f));
		checkf_small("geom_matinv4f and geom_matmat4f check6" , fd*(f3[ 5] - 1.f));
		checkf_small("geom_matinv4f and geom_matmat4f check7" , fd*(f3[ 6] - 0.f));
		checkf_small("geom_matinv4f and geom_matmat4f check8" , fd*(f3[ 7] - 0.f));
		checkf_small("geom_matinv4f and geom_matmat4f check9" , fd*(f3[ 8] - 0.f));
		checkf_small("geom_matinv4f and geom_matmat4f check10", fd*(f3[ 9] - 0.f));
		checkf_small("geom_matinv4f and geom_matmat4f check11", fd*(f3[10] - 1.f));
		checkf_small("geom_matinv4f and geom_matmat4f check12", fd*(f3[11] - 0.f));
		checkf_small("geom_matinv4f and geom_matmat4f check13", fd*(f3[12] - 0.f));
		checkf_small("geom_matinv4f and geom_matmat4f check14", fd*(f3[13] - 0.f));
		checkf_small("geom_matinv4f and geom_matmat4f check15", fd*(f3[14] - 0.f));
		checkf_small("geom_matinv4f and geom_matmat4f check16", fd*(f3[15] - 1.f));
		checkd_small("geom_matinv4d and geom_matmat4d check1" , dd*(d3[ 0] - 1.));
		checkd_small("geom_matinv4d and geom_matmat4d check2" , dd*(d3[ 1] - 0.));
		checkd_small("geom_matinv4d and geom_matmat4d check3" , dd*(d3[ 2] - 0.));
		checkd_small("geom_matinv4d and geom_matmat4d check4" , dd*(d3[ 3] - 0.));
		checkd_small("geom_matinv4d and geom_matmat4d check5" , dd*(d3[ 4] - 0.));
		checkd_small("geom_matinv4d and geom_matmat4d check6" , dd*(d3[ 5] - 1.));
		checkd_small("geom_matinv4d and geom_matmat4d check7" , dd*(d3[ 6] - 0.));
		checkd_small("geom_matinv4d and geom_matmat4d check8" , dd*(d3[ 7] - 0.));
		checkd_small("geom_matinv4d and geom_matmat4d check9" , dd*(d3[ 8] - 0.));
		checkd_small("geom_matinv4d and geom_matmat4d check10", dd*(d3[ 9] - 0.));
		checkd_small("geom_matinv4d and geom_matmat4d check11", dd*(d3[10] - 1.));
		checkd_small("geom_matinv4d and geom_matmat4d check12", dd*(d3[11] - 0.));
		checkd_small("geom_matinv4d and geom_matmat4d check13", dd*(d3[12] - 0.));
		checkd_small("geom_matinv4d and geom_matmat4d check14", dd*(d3[13] - 0.));
		checkd_small("geom_matinv4d and geom_matmat4d check15", dd*(d3[14] - 0.));
		checkd_small("geom_matinv4d and geom_matmat4d check16", dd*(d3[15] - 1.));
	}
	printf("Failures in matrix inversion routines may be ignored if they are rare.\n");
	printf("(Randomly generated matrices may be nearly singular.)\n");
	
	printf("Checking matsvd\n");
	for(i = 0; i < 2; ++i){
		if(0 == i){
			d1[0] = 1;
			d1[1] = 0;
			d1[2] = 0;
			d1[3] = 1;
		}else{
			geom_randd(4, d1);
		}
		geom_matsvd2d(d1, &d2[0], d3, &d2[4]);
		// Check for orthogonality
		checkd_small("geom_matsvd2d orth u check1", dd*(d2[0]*d2[0] + d2[1]*d2[1] - 1.));
		checkd_small("geom_matsvd2d orth u check2", dd*(d2[0]*d2[1] + d2[2]*d2[3] - 0.));
		checkd_small("geom_matsvd2d orth u check3", dd*(d2[0]*d2[2] + d2[1]*d2[3] - 0.));
		checkd_small("geom_matsvd2d orth u check4", dd*(d2[2]*d2[2] + d2[3]*d2[3] - 1.));
		checkd_small("geom_matsvd2d orth v check1", dd*(d2[4]*d2[4] + d2[5]*d2[5] - 1.));
		checkd_small("geom_matsvd2d orth v check2", dd*(d2[4]*d2[5] + d2[6]*d2[7] - 0.));
		checkd_small("geom_matsvd2d orth v check3", dd*(d2[4]*d2[6] + d2[5]*d2[7] - 0.));
		checkd_small("geom_matsvd2d orth v check4", dd*(d2[6]*d2[6] + d2[7]*d2[7] - 1.));
		d2[0] *= d3[0];
		d2[1] *= d3[0];
		d2[2] *= d3[1];
		d2[3] *= d3[1];
		geom_matmat2d(&d2[0], &d2[4], d4);
		// Check that the decomposition is correct
		checkd_small("geom_matsvd2d check1", dd*(d4[0] - d1[0]));
		checkd_small("geom_matsvd2d check2", dd*(d4[1] - d1[1]));
		checkd_small("geom_matsvd2d check3", dd*(d4[2] - d1[2]));
		checkd_small("geom_matsvd2d check4", dd*(d4[3] - d1[3]));
	}
	return 0;
}

