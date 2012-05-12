#include <geom_poly.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double frand(){
	return (double)rand() / (double)RAND_MAX - 0.5;
}

int main(int argc, char **argv) {
	int num_iters;
	srand(23);

	unsigned int i, n = 11;
	double c[3];
	double x[3];
	double pn[48];
	c[0] = frand();
	c[1] = frand();
	c[2] = frand();
	for(i = 0; i < n; ++i){
		pn[4*i+0] = frand();
		pn[4*i+1] = frand();
		pn[4*i+2] = frand();
		pn[4*i+3] = fabs(frand());
	}

	if(geom_convex_bound3d(n, pn, c, x, NULL)){
		printf("unbounded\n");
	}
	printf("%g %g %g\n", x[0], x[1], x[2]);

	return 0;
}
