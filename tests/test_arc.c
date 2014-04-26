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

void draw_arc(const double a[2], const double b[2], double g){
	unsigned int i;
	// Draw discretization
	const unsigned int k = 100;
	for(i = 0; i <= k; ++i){
		double s = (double)i/(double)k;
		double p[2];
		geom_arc_param(a, b, g, s, p, NULL);
		if(0 == i){
			printf("newpath %g %g moveto ", p[0], p[1]);
		}else{
			printf("%g %g lineto ", p[0], p[1]);
		}
	}
	printf("stroke\n");
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
		"4 4 scale\n"
		"0.25 0.75 translate\n"
		"0.001 setlinewidth\n"
		"/Courier findfont 0.03 scalefont setfont\n"
	);
	
	for(i = 0; i < 1; ++i){
		double a[2], b[2], c[2], g;
		for(j = 0; j <= 0; ++j){
			a[0] = frand();
			a[1] = frand();
			b[0] = frand();
			b[1] = frand();
			c[0] = frand();
			c[1] = frand();
		}
		
		g = exp(8.*(frand() - 0.5));
		if(1 & rand()){
			g = -g;
		}
		
		
		printf("%g %g moveto %g %g 0.002 0 360 arc fill\n",
			c[0], c[1], c[0], c[1]
		);
		
		//g = geom_arc_g_from_pt(a, b, c);
		
		printf("newpath %g %g moveto (a) show\n", a[0], a[1]);
		printf("newpath %g %g moveto (b) show\n", b[0], b[1]);
		printf("newpath %g %g moveto (c) show\n", c[0], c[1]);
		
		printf("gsave 0 -0.5 translate\n");
		printf("0 0.09 moveto (a = %g, %g) show\n", a[0], a[1]);
		printf("0 0.06 moveto (b = %g, %g) show\n", b[0], b[1]);
		printf("0 0.03 moveto (c = %g, %g) show\n", c[0], c[1]);
		printf("0 0 moveto (g = %g) show grestore\n", g);
		
		if(0){
			printf("%g %g moveto %g %g lineto stroke\n",
				a[0], a[1], b[0], b[1]
			);
		}
		
		if(0){
			double xb[2], yb[2];
			geom_arc_bound_rect(a, b, g, xb, yb);
			printf("%g %g moveto %g %g lineto %g %g lineto %g %g lineto closepath stroke\n",
				xb[0], yb[0],
				xb[1], yb[0],
				xb[1], yb[1],
				xb[0], yb[1]
			);
		}
		
		if(0){
			double c[2], r;
			geom_arc_bound_circle(a, b, g, c, &r);
			printf("gsave 1 0 0 setrgbcolor\n");
			printf("%g %g moveto %g %g %g 0 360 arc closepath stroke\n", c[0]+r, c[1], c[0], c[1], r);
			printf("grestore\n");
		}
		
		// Draw discretization
		printf("0.01 setlinewidth\n");
		printf("1 0 0 setrgbcolor\n");
		draw_arc(a, b, g);
		printf("0.001 setlinewidth\n");
		printf("0 setgray\n");
		
		if(0){
			// Draw tangents
			for(i = 0; i <= k; ++i){
				double s = (double)i/(double)k;
				double p[2], t[2];
				geom_arc_param(a, b, g, s, p, t);
				printf("newpath %g %g 0.002 0 360 arc fill\n", p[0], p[1]);
				printf("newpath %g %g moveto %g %g rlineto stroke\n", p[0], p[1], t[0], t[1]);
			}
		}
		
		
		if(0){
			double c[2], theta[2];
			double r;
			geom_arc_circle(a, b, g, c, &r, theta);
			printf("newpath %g %g %g %g %g arc%c stroke\n",
				c[0], c[1], r * 0.9, theta[0]*180/M_PI, theta[1]*180/M_PI,
				(g < 0 ? 'n' : ' ')
			);
		}
		
		
		if(0){
			// Draw the monotone splitting
			double pg[17];
			int n = geom_arc_split_monotone(a, b, g, pg);
			for(i = 0; i < n; ++i){
				double c[2];
				double r;
				double theta[2];
				geom_arc_circle(&pg[3*i+0], &pg[3*(i+1)+0], pg[3*i+2], c, &r, theta);
				printf(" newpath %g %g %g %g %g arc%c stroke\n",
					c[0], c[1], r * (0.7 + 0.05*i),
					theta[0]*180/M_PI, theta[1]*180/M_PI,
					(g < 0 ? 'n' : ' ')
				);
			}
			for(i = 0; i <= n; ++i){
				printf("newpath %g %g 0.005 0 360 arc fill\n", pg[3*i+0], pg[3*i+1]);
			}
		}
		
		if(1){
			double ao[2], bo[2], go = g;
			geom_arc_offset(a, b, g, 0.1, ao, bo);
			draw_arc(ao, bo, go);
			geom_arc_offset(a, b, g, -0.1, ao, bo);
			draw_arc(ao, bo, go);
		}
		
		if(1){
			double ao[2], bo[2], go = g;
			double d[2] = { 0.1, 0.2 };
			geom_arc_extend(a, b, g, d, ao, bo, &go);
			draw_arc(ao, bo, go);
		}
	}
	
	
	printf("showpage\n");
	
	return 0;
}

