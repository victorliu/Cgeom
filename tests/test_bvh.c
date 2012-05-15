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

int traverse(int tag, const double c[2], const double h[2], int leaf, void *data){
	FILE *fp = (FILE*)data;
	if(leaf){
		fprintf(fp, "0 setgray\n");
	}else{
		fprintf(fp, "0.5 setgray\n");
	}
	fprintf(fp, "newpath %f %f moveto %f %f lineto %f %f lineto %f %f lineto closepath stroke\n",
		c[0]-h[0], c[1]-h[1],
		c[0]+h[0], c[1]-h[1],
		c[0]+h[0], c[1]+h[1],
		c[0]-h[0], c[1]+h[1]
		);
	fprintf(fp, "newpath %f %f moveto (%d) show\n", c[0], c[1], tag);
	return 1;
}

int query(int tag, const double c[2], const double h[2], void *data){
	FILE *fp = (FILE*)data;
	fprintf(fp, "newpath %f %f moveto %f %f lineto %f %f lineto %f %f lineto closepath stroke\n",
		c[0]-h[0], c[1]-h[1],
		c[0]+h[0], c[1]-h[1],
		c[0]+h[0], c[1]+h[1],
		c[0]-h[0], c[1]+h[1]
		);
	return 1;
}

int main(){
	const int n = 100;
	int tag = 0;
	FILE *fp = fopen("test_bvh2d.ps", "wb");
	geom_bvh2d bvh = geom_bvh2d_new(n, &gen_shape, (void*)&tag);
	
	srand(0);
	fprintf(fp, "32 32 scale\n5 5 translate\n0.01 setlinewidth\n/Helvetica findfont 0.1 scalefont setfont\n");
	geom_bvh2d_traverse(bvh, &traverse, (void*)fp);
	
	fprintf(fp, "\n1 0 0 setrgbcolor\n");
	{
		double p[2];
		geom_randd(2, p); p[0] *= 3;  p[1] *= 3;
		fprintf(fp, "newpath %f %f 0.1 0 360 arc fill\n", p[0], p[1]);
		geom_bvh2d_query_pt(bvh, p, &query, (void*)fp);
	}
	
	fprintf(fp, "showpage\n");
	
	geom_bvh2d_destroy(bvh);
	fclose(fp);
	return 0;
}
