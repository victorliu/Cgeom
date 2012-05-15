#include <geom_bvh.h>
#include <test_common.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int gen_shape(double c[2], double h[2], int *tag, void *data){
	int *t = (int*)data;
	*tag = *t;
	(*t)++;
	
	geom_randd(2, c); c[0] *= 8; c[1] *= 8;
	geom_randd(2, h);
	h[0] = 2*fabs(h[0]);
	h[1] = 2*fabs(h[1]);

	return 0;
}

int main(){
	int t;
	for(t = 0; t < 99999; ++t){
		const int n = 100;
		int tag = 0;
		geom_bvh2d bvh = geom_bvh2d_new(n, &gen_shape, (void*)&tag);
		geom_bvh2d_destroy(bvh);
	}
	return 0;
}
