#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <float.h>

void geom_randf(unsigned int n, float  *v){
	while(n --> 0){
		*v = ((float)rand() / (float)RAND_MAX) - 0.5f;
		++v;
	}
}
void geom_randd(unsigned int n, double *v){
	while(n --> 0){
		*v = ((double)rand() / (double)RAND_MAX) - 0.5f;
		++v;
	}
}

float  geom_maxf(unsigned int n, float  *v){
	float m = 0;
	while(n --> 0){
		float a = fabsf(*v);
		if(a > m){ m = a; }
		++v;
	}
	return m;
}
double geom_maxd(unsigned int n, double *v){
	double m = 0;
	while(n --> 0){
		double a = fabs(*v);
		if(a > m){ m = a; }
		++v;
	}
	return m;
}

void geom_setzerof(unsigned int n, float  *v){
	while(n --> 0){
		*v = 0.f;
		++v;
	}
}
void geom_setzerod(unsigned int n, double *v){
	while(n --> 0){
		*v = 0.;
		++v;
	}
}


int checkf_small(const char *func, float v){
	v = fabsf(v);
	if(v > FLT_EPSILON){
		printf(" Failure in %s: %.14g > %.14g\n", func, v, FLT_EPSILON);
		return 1;
	}
	return 0;
}
int checkd_small(const char *func, double v){
	v = fabs(v);
	if(v > DBL_EPSILON){
		printf(" Failure in %s: %.14g > %.14g\n", func, v, DBL_EPSILON);
		return 1;
	}
	return 0;
}
int checkf_pos(const char *func, float v){
	if(v <= 0.f){
		printf(" Failure in %s: %.14g not positive\n", func, v);
		return 1;
	}
	return 0;
}
int checkd_pos(const char *func, double v){
	if(v <= 0.){
		printf(" Failure in %s: %.14g not positive\n", func, v);
		return 1;
	}
	return 0;
}
int checkf_equal(const char *func, float v, float w){
	if(v != w){
		printf(" Failure in %s: %.14g != %.14g\n", func, v, w);
		return 1;
	}
	return 0;
}
int checkd_equal(const char *func, double v, double w){
	if(v != w){
		printf(" Failure in %s: %.14g != %.14g\n", func, v, w);
		return 1;
	}
	return 0;
}
int checkf_NaN(const char *func, float v){
	if(v == v){
		printf(" Failure in %s: %.14g != NaN, %d\n", func, v);
		return 1;
	}
	return 0;
}
int checkd_NaN(const char *func, double v){
	if(v == v){
		printf(" Failure in %s: %.14g != NaN\n", func, v);
		return 1;
	}
	return 0;
}
