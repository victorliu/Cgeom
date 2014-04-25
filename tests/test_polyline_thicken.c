#include <geom_arclinegraph.h>
#include <stdio.h>

int main(){
	unsigned int n = 9;
	unsigned int i, j;
	double v[18] = {
		0,0,
		2,0,
		2,4,
		4,4,
		4,-4,
		6,-4,
		6,8,
		8,8,
		8,-8
	};
	double vo[36];
	
	printf("16 16 scale\n8 20 translate\n0.01 setlinewidth\n");
	
	printf("1 0 0 setrgbcolor\n");
	printf("newpath %g %g moveto\n", v[0], v[1]);
	for(i = 1; i < n; ++i){
		printf("%g %g lineto\n", v[2*i+0], v[2*i+1]);
	}
	printf("stroke\n");
	polyline_thicken(n, v, 1, vo);
	
	printf("0 0 0 setrgbcolor\n");
	printf("newpath %g %g moveto\n", vo[0], vo[1]);
	for(i = 1; i < 2*n; ++i){
		printf("%g %g lineto\n", vo[2*i+0], vo[2*i+1]);
	}
	printf("closepath stroke\n");
	
	// reflect over 0,1
	for(i = 0; i < n; ++i){
		v[2*i+0] = 0 - v[2*i+0];
		v[2*i+1] = 2 - v[2*i+1];
	}
	printf("1 0 0 setrgbcolor\n");
	printf("newpath %g %g moveto\n", v[0], v[1]);
	for(i = 1; i < n; ++i){
		printf("%g %g lineto\n", v[2*i+0], v[2*i+1]);
	}
	printf("stroke\n");
	polyline_thicken(n, v, 1, vo);
	
	printf("0 0 0 setrgbcolor\n");
	printf("newpath %g %g moveto\n", vo[0], vo[1]);
	for(i = 1; i < 2*n; ++i){
		printf("%g %g lineto\n", vo[2*i+0], vo[2*i+1]);
	}
	printf("closepath stroke\n");
	
	// draw orientation arrows
	printf("0.02 setlinewidth\n");
	for(j = 2*n-1, i = 0; i < 2*n; j=i++){
		double c[2] = { 0.5*(vo[2*j+0]+vo[2*i+0]), 0.5*(vo[2*j+1]+vo[2*i+1]) };
		printf("newpath %g %g moveto\n", vo[2*j+0], vo[2*j+1]);
		printf("%g %g lineto stroke\n", c[0], c[1]);
	}
	
	printf("showpage\n");
	
	return 0;
}
