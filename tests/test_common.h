#ifndef GEOM_TEST_COMMON_H_INCLUDED
#define GEOM_TEST_COMMON_H_INCLUDED

void geom_randf(unsigned int n, float  *v);
void geom_randd(unsigned int n, double *v);
float  geom_maxf(unsigned int n, float  *v);
double geom_maxd(unsigned int n, double *v);
void geom_setzerof(unsigned int n, float  *v);
void geom_setzerod(unsigned int n, double *v);

#endif // GEOM_TEST_COMMON_H_INCLUDED
