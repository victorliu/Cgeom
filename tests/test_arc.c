#include <geom_arc.h>
#include <test_common.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double frand(){
	static const double nrm = 1. / (double)RAND_MAX;
	double f1 = (double)rand() * nrm;
	double f2 = (double)rand() * nrm;
	return f1 + f2*nrm;
}

int main(int argc, char *argv[]){
	unsigned int i, j, k;
	
	if(argc < 2){
		srand(0);
	}else{
		srand(atoi(argv[1]));
	}
	
	printf(
		"72 72 scale\n"
		"8 8 scale\n"
		"0.0625 0.375 translate\n"
		"0.001 setlinewidth\n"
		"/Courier findfont 0.03 scalefont setfont\n"
	);
	
	for(i = 0; i < 1; ++i){
		double a[2] = { frand(), frand() };
		double b[2] = { frand(), frand() };
		double g = exp(8.*(frand() - 0.5));
		if(1 & rand()){
			g = -g;
		}
		
		printf("0 0.06 moveto (a = %g, %g) show\n", a[0], a[1]);
		printf("0 0.03 moveto (b = %g, %g) show\n", b[0], b[1]);
		printf("0 0 moveto (g = %g) show\n", g);
		
		printf("%g %g moveto %g %g lineto stroke\n",
			a[0], a[1], b[0], b[1]
		);
		
		{
			double xb[2], yb[2];
			geom_arc_bound_rect(a, b, g, xb, yb);
			printf("%g %g moveto %g %g lineto %g %g lineto %g %g lineto closepath stroke\n",
				xb[0], yb[0],
				xb[1], yb[0],
				xb[1], yb[1],
				xb[0], yb[1]
			);
		}
		
		{
			double c[2], r;
			geom_arc_bound_circle(a, b, g, c, &r);
			printf("gsave 1 0 0 setrgbcolor\n");
			printf("%g %g moveto %g %g %g 0 360 arc closepath stroke\n", c[0]+r, c[1], c[0], c[1], r);
			printf("grestore\n");
		}
		
		k = 100;
		for(i = 0; i <= k; ++i){
			double s = (double)i/(double)k;
			double p[2];
			geom_arc_param(a, b, g, s, p, NULL);
			if(0 == i){
				printf("%g %g moveto ", p[0], p[1]);
			}else{
				printf("%g %g lineto ", p[0], p[1]);
			}
		}
		printf("stroke\n");
		
		
	}
	
	
	printf("showpage\n");
	
	return 0;
}

