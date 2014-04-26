#include <stdlib.h>
#include <stdio.h>
#include <Cgeom/geom_la.h>
#include <Cgeom/geom_arc.h>

/* Given a simply connected polyline (not closed) defined by a
 * sequence of n points (and so there are n-1 segments), produce
 * a closed, positively (CCW) oriented polygon that is the offset
 * curve of the polyline with a specified total thickness.
 * This function assumes the input does not have adjacent segments
 * that double back on themselves.
 */
int polyline_thicken(unsigned int n, const double *vin, double thickness, double *vout){
	const double h = 0.5 * thickness;
	double tc[2]; /* current directed tangent */
	double tp[2]; /* previous directed tangent */
	unsigned int i;
	if(0 == n){ return 0; }
	if(1 == n){
		vout[0] = vin[0] - h;
		vout[1] = vin[1] - h;
		vout[2] = vin[2] + h;
		vout[3] = vin[3] - h;
		vout[4] = vin[4] + h;
		vout[5] = vin[5] + h;
		vout[6] = vin[6] - h;
		vout[7] = vin[7] + h;
		return 4;
	}
	
	/* right side normal: tc[1], -tc[0] */
	/* left side normal: -tc[1], tc[0] */
	
	/* for each point, we add the right side offset, then at index j add the left side offset */
	
	/* process initial point */
	{
		const unsigned int j = 2*n-1;
		tc[0] = vin[2]-vin[0]; tc[1] = vin[3]-vin[1]; geom_normalize2d(tc);
		vout[0] = vin[0] - h*(tc[0]-tc[1]);
		vout[1] = vin[1] - h*(tc[1]+tc[0]);
		vout[2*j+0] = vin[0] - h*(tc[0]+tc[1]);
		vout[2*j+1] = vin[1] - h*(tc[1]-tc[0]);
	}
	/* process middle points */
	for(i = 1; i+1 < n; ++i){
		const unsigned int j = 2*n-i-1;
		tp[0] = tc[0]; tp[1] = tc[1];
		tc[0] = vin[2*i+2]-vin[2*i+0]; tc[1] = vin[2*i+3]-vin[2*i+1]; geom_normalize2d(tc);
		/* Compute intersection of two line offsets
		 * Let the right side normals before and after point p be u and v.
		 * The desired point q = p + alpha*(u+v) where alpha = 1/(1+u.v)
		 */
		const double dot = tc[0]*tp[0]+tc[1]*tp[1];
		double alpha = h/(1.+dot);
		vout[2*i+0] = vin[2*i+0] + alpha*(tc[1]+tp[1]);
		vout[2*i+1] = vin[2*i+1] - alpha*(tc[0]+tp[0]);
		/* For the left side normal, we want q = p - alpha*(u+v) */
		vout[2*j+0] = vin[2*i+0] - alpha*(tc[1]+tp[1]);
		vout[2*j+1] = vin[2*i+1] + alpha*(tc[0]+tp[0]);
	}
	/* process final point */
	{
		i = n-1;
		vout[2*i+0] = vin[2*i+0] + h*(tc[0]+tc[1]);
		vout[2*i+1] = vin[2*i+1] + h*(tc[1]-tc[0]);
		vout[2*i+2] = vin[2*i+0] + h*(tc[0]-tc[1]);
		vout[2*i+3] = vin[2*i+1] + h*(tc[1]+tc[0]);
	}
	return 2*n;
}

int arclinegraph_thicken(
	int join_type, const double *v,
	/* input arcs: */
	unsigned int n, const unsigned int *iab, const double *g,
	double thickness,
	/* outputs: */
	double *w, unsigned int *icd, double *h,
	int *nw, int *nseg
){
	unsigned int i, j;
	const int n2 = 2*n;
	double dum[2];
	/* First generate the halfedge structures.
	 * half: { from, next, flip, edge }
	 * We don't need flip, since flip(i) == i^1
	 * We don't need edge, since edge(i) == i/2
	 */
	int *half = (int*)malloc(sizeof(int) * 2*n2);
	for(i = 0; i < n; ++i){
		half[2*(2*i+0)+0] = iab[2*i+0];
		half[2*(2*i+0)+1] = 2*i+1;
		half[2*(2*i+1)+0] = iab[2*i+1];
		half[2*(2*i+1)+1] = 2*i+0;
	}
	/* Generate tangents to endpoints of arcs */
	double *tv = (double*)malloc(sizeof(double) * 2*n2);
	for(i = 0; i < n; ++i){
		geom_arc_param(&v[2*iab[2*i+0]], &v[2*iab[2*i+1]], g[i], 0., dum, &tv[2*(2*i+0)]);
		geom_normalize2d(&tv[2*(2*i+0)]);
		geom_arc_param(&v[2*iab[2*i+0]], &v[2*iab[2*i+1]], g[i], 1., dum, &tv[2*(2*i+1)]);
		geom_normalize2d(&tv[2*(2*i+1)]);
		tv[2*(2*i+1)+0] = -tv[2*(2*i+1)+0];
		tv[2*(2*i+1)+1] = -tv[2*(2*i+1)+1];
	}
	/* tangents point towards middle of arc (inward instead of outward) */
	
	/* Join the halfedges by finding the next element */
	for(i = 0; i < n2; ++i){
		/* skip if already set */
		if(half[2*i+1] >= 0){ continue; }
		const int ito = half[2*(i^1)+0]; /* vertex index of dst */
		/* Get the incident tangent (pointing toward dst) */
		double bsin = 0;
		double bcos = 1;
		const unsigned int iopp = i^1;
		for(j = 0; j < n2; ++j){
			if((ito != half[2*j+0]) || (j == iopp)){ continue; }
			/* Get the outgoing tangent (pointing away from dst(i)) */
			double csin = tv[2*(2*i+1)+0]*tv[2*(2*j+0)+1] - tv[2*(2*i+1)+1]*tv[2*(2*j+0)+0];
			double ccos = tv[2*(2*i+1)+0]*tv[2*(2*j+0)+0] + tv[2*(2*i+1)+1]*tv[2*(2*j+0)+1];
			/* Determine if jtan is pointing more CW than current best */
			int is_better = 0;
			if(iopp == half[2*i+1]){
				is_better = 1;
			}else{
				if(bsin > 0){
					is_better = (csin > 0 && ccos > bcos);
				}else{
					is_better = (csin > 0 || ccos < bcos);
				}
			}
			if(is_better){
				half[2*i+1] = j;
				bsin = csin; bcos = ccos;
			}
		}
	}
	/* At this point the arrangement has been computed */
	for(i = 0; i < n2; ++i){
		if(half[2*i+1] < 0){
			fprintf(stderr, "BAD: i = %d\n", i);
		}
	}
	
	/* Now walk the halfedges and produce offsets */
	unsigned int ks = 0; /* Running index of end of icd and h */
	unsigned int kw = 0; /* RUnning index of end of w */
	for(i = n2-1, j = 0; j < n2; i=j++){
		const unsigned int ie = i/2;
		const unsigned int je = j/2;
		const unsigned io = i&1;
		const unsigned jo = j&1;
		double dot = tv[2*i+0]*tv[2*j+0] + tv[2*i+1]*tv[2*j+1];
		if(half[2*i+1] == j || 1 == dot){ /* endpoint cap */
			
		}else{
		}
	}
	
	free(tv);
	free(half);
	return 0;
}
