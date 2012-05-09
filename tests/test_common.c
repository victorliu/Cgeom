#include <stdlib.h>
#include <math.h>

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
