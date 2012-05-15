#include <geom_shapeset.h>
#include <stdlib.h>
#include <stdio.h>

int main(){
	double x, y;
	FILE *fp = fopen("test_shapeset2d.dat", "wb");
	const double L[4] = {
		1,0,
		0,1
	};
	int i, j;
	
	geom_shapeset2d ss = geom_shapeset2d_new(L);
	
	int n = 10;
	for(i = 0; i < n; ++i){
		for(j = 0; j < n; ++j){
			const double r = 0.2;
			geom_shape2d *s = (geom_shape2d*)malloc(sizeof(geom_shape2d));
			s->type = GEOM_SHAPE2D_ELLIPSE;
			s->s.ellipse.A[0] = r;
			s->s.ellipse.A[1] = 0;
			s->s.ellipse.A[2] = 0;
			s->s.ellipse.A[3] = r;
			s->s.ellipse.B[0] = 1./r;
			s->s.ellipse.B[1] = 0;
			s->s.ellipse.B[2] = 0;
			s->s.ellipse.B[3] = 1./r;
			s->org[0] = 0.1*i;
			s->org[1] = 0.1*j;
			geom_shapeset2d_add(ss, s);
		}
	}
	geom_shapeset2d_finalize(ss);
	
	for(x = -0.5; x <= 1.5; x += 0.005){
		for(y = -0.5; y <= 1.5; y += 0.005){
			double p[2] = {x, y};
			double v = -1;
			int index = geom_shapeset2d_query_pt(ss, p);
			if(index >= 0){
				v = (double)index/(double)(n*n);
			}
			fprintf(fp, "%f %f %f %d\n", x, y, v, index);
		}
		fprintf(fp, "\n");
	}
	
	n = geom_shapeset2d_size(ss);
	fprintf(stderr, "n = %d\n", n);
	for(i = 0; i < n; ++i){
		geom_shape2d *s = geom_shapeset2d_index(ss, i);
		free(s);
	}
	geom_shapeset2d_destroy(ss);
	
	fclose(fp);
	return 0;
}
