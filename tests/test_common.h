#ifndef GEOM_TEST_COMMON_H_INCLUDED
#define GEOM_TEST_COMMON_H_INCLUDED

void geom_randf(unsigned int n, float  *v);
void geom_randd(unsigned int n, double *v);
float  geom_maxf(unsigned int n, float  *v);
double geom_maxd(unsigned int n, double *v);
void geom_setzerof(unsigned int n, float  *v);
void geom_setzerod(unsigned int n, double *v);

int checkf_small(const char *func, float v);
int checkd_small(const char *func, double v);
int checkf_pos(const char *func, float v);
int checkd_pos(const char *func, double v);
int checkf_equal(const char *func, float v, float w);
int checkd_equal(const char *func, double v, double w);
int checkf_NaN(const char *func, float v);
int checkd_NaN(const char *func, double v);

#endif // GEOM_TEST_COMMON_H_INCLUDED
