#include <geom_bvh.h>
#include <stdlib.h>

//#define BVH_DEBUG

#ifdef BVH_DEBUG
# include <stdio.h>
# define BVHDBG(...) fprintf(stderr, __VA_ARGS__)
#else
# define BVHDBG(...)
#endif

// The internal representation of boxes is min/max instead of center/halfwidth

struct geom_bvh2d_struct{
	geom_bvh2d child[4];
	double b[4]; // x-min, x-max, y-min, y-max
	int tag; // for internal nodes, set to max of all subnodes
};
struct geom_bvh3d_struct{
	geom_bvh3d child[8];
	double b[6]; // x-min, x-max, y-min, y-max, z-min, z-max
	int tag;
};

static int isqrt_ceil(int u){
	if(u <= 0){ return 0; }
	int a = 2;
	int b = u/a;
	while(((a-b) > 1) || ((b-a) > 1)){
		a = (a+b)/2;
		b = u/a;
	}
	a = (a<b)?a:b;
	return (a*a == u) ? a : a+1;
}
static int icbrt_ceil(int x){
	const int xsave = x;
	int s = 30;
	int y = 0, b;
	// the 3 add's of 1 can be replaced by or's of 1 since the arguments are always even
	while(s >= 0){ // Do 11 times.
		y *= 2;
		//b = (3*y*(y + 1) + 1) << s;
		b = ((3*y*(y | 1)) | 1) << s;
		s -= 3;
		if(x >= b){
			x -= b;
			y |= 1; //y++;
		}
	}
	return (y*y*y == xsave) ? y : y+1;
}

static void sort_bvh2(geom_bvh2d *B, int beg, int end, int d){
	if(end > beg+1){
		geom_bvh2d t, piv = B[beg];
		int l = beg+1, r = end;
		while(l < r){
			if(B[l]->b[2*d+0]+B[l]->b[2*d+1] <= piv->b[2*d+0]+piv->b[2*d+1]){
				++l;
			}else{
				--r;
				t = B[l];
				B[l] = B[r];
				B[r] = t;
			}
		}
		--l;
		t = B[l];
		B[l] = B[beg];
		B[beg] = t;
		sort_bvh2(B, beg, l, d);
		sort_bvh2(B, r, end, d);
	}
}
static void sort_bvh3(geom_bvh3d *B, int beg, int end, int d){
	if(end > beg+1){
		geom_bvh3d t, piv = B[beg];
		int l = beg+1, r = end;
		while(l < r){
			if(B[l]->b[2*d+0]+B[l]->b[2*d+1] <= piv->b[2*d+0]+piv->b[2*d+1]){
				++l;
			}else{
				--r;
				t = B[l];
				B[l] = B[r];
				B[r] = t;
			}
		}
		--l;
		t = B[l];
		B[l] = B[beg];
		B[beg] = t;
		sort_bvh3(B, beg, l, d);
		sort_bvh3(B, r, end, d);
	}
}

static geom_bvh2d geom_bvh2d_make_blank(){
	geom_bvh2d b = (geom_bvh2d)malloc(sizeof(struct geom_bvh2d_struct));
	int j;
	for(j = 0; j < 4; ++j){
		b->child[j] = NULL;
	}
	b->tag = 0;
	b->b[0] = 0;
	b->b[1] = 0;
	b->b[2] = 0;
	b->b[3] = 0;
	return b;
}
static geom_bvh3d geom_bvh3d_make_blank(){
	geom_bvh3d b = (geom_bvh3d)malloc(sizeof(struct geom_bvh3d_struct));
	int j;
	for(j = 0; j < 8; ++j){
		b->child[j] = NULL;
	}
	b->tag = 0;
	b->b[0] = 0;
	b->b[1] = 0;
	b->b[2] = 0;
	b->b[3] = 0;
	b->b[4] = 0;
	b->b[5] = 0;
	return b;
}

static void geom_bvh2d_STR(unsigned int *n_, geom_bvh2d *B){
	int i, j, k;
	geom_bvh2d *T = B;
	const int n = *n_;
	const int P = (n+3)/4;
	const int S = isqrt_ceil(P);
	*n_ = 0;
	
	// Sort rectangles by first coordinate
	sort_bvh2(B, 0, n, 0);
	
	// Partition into S slices and sort each slice by second coordinate
	i = n;
	while(i > 0){
		int slice_size = 4*S;
		if(slice_size > i){ slice_size = i; }
		BVHDBG(" Slice size = %d\n", slice_size);
		
		sort_bvh2(T, 0, slice_size, 1);
		
		// Now pack all the nodes into runs of length 4
		j = slice_size;
		while(j > 0){
			int pack_size = 4;
			if(pack_size > j){ pack_size = j; }
			
			// Make a new node
			geom_bvh2d b = geom_bvh2d_make_blank();
			b->child[0] = T[0];
			b->tag = T[0]->tag;
			b->b[0] = T[0]->b[0];
			b->b[1] = T[0]->b[1];
			b->b[2] = T[0]->b[2];
			b->b[3] = T[0]->b[3];
			for(k = 1; k < pack_size; ++k){
				b->child[k] = T[k];
				if(T[k]->tag > b->tag){ b->tag = T[k]->tag; }
				if(T[k]->b[0] < b->b[0]){ b->b[0] = T[k]->b[0]; }
				if(T[k]->b[1] > b->b[1]){ b->b[1] = T[k]->b[1]; }
				if(T[k]->b[2] < b->b[2]){ b->b[2] = T[k]->b[2]; }
				if(T[k]->b[3] > b->b[3]){ b->b[3] = T[k]->b[3]; }
			}
			BVHDBG("  Making internal node %p, tag=%d, b=%f,%f,%f,%f\n", b, b->tag, b->b[0], b->b[1], b->b[2], b->b[3]);
			for(k = 0; k < pack_size; ++k){
				BVHDBG("   Child node %p, tag=%d, b=%f,%f,%f,%f\n", T[k], T[k]->tag, T[k]->b[0], T[k]->b[1], T[k]->b[2], T[k]->b[3]);
			}
			*B = b;
			T += pack_size;
			B++;
			(*n_)++;
			
			j -= pack_size;
		}
		
		i -= slice_size;
	}
	BVHDBG("Iteration of STR done; n=%d -> %d\n", n, *n_);
}
geom_bvh2d geom_bvh2d_new(unsigned int n, int (*shape_iterator)(double c[2], double h[2], int *tag, void *data), void *data){
	unsigned int i;
	
	geom_bvh2d *B = (geom_bvh2d*)malloc(sizeof(geom_bvh2d) * n);
	geom_bvh2d ret;
	for(i = 0; i < n; ++i){
		double c[2], h[2];
		B[i] = geom_bvh2d_make_blank();
		shape_iterator(c, h, &(B[i]->tag), data);
		B[i]->b[0] = c[0]-h[0];
		B[i]->b[1] = c[0]+h[0];
		B[i]->b[2] = c[1]-h[1];
		B[i]->b[3] = c[1]+h[1];
		BVHDBG("Making leaf %p, tag=%d, b=%f,%f,%f,%f\n", B[i], B[i]->tag, B[i]->b[0], B[i]->b[1], B[i]->b[2], B[i]->b[3]);
	}
	
	while(n > 1){
		geom_bvh2d_STR(&n, B);
	}
	ret = B[0];
	free(B);
	return ret;
}

static void geom_bvh3d_STR(unsigned int *n_, geom_bvh3d *B){
	int i, j, k, l;
	geom_bvh3d *T = B;
	const unsigned int n = *n_;
	const unsigned int P = (n+7)/8;
	const unsigned int S = icbrt_ceil(P);
	*n_ = 0;
	
	// Sort rectangles by first coordinate
	sort_bvh3(B, 0, n, 0);
	
	// Partition into S slices and sort each slice by second coordinate
	i = n;
	while(i > 0){
		int slice_size = 8*S;
		if(slice_size > i){ slice_size = i; }
		BVHDBG(" Slice size = %d\n", slice_size);
		
		sort_bvh3(T, 0, slice_size, 1);
		
		// Now slice by second coordinate and sort by third coord
		j = slice_size;
		while(j > 0){
			int slice_size2 = 8*S;
			if(slice_size2 > j){ slice_size2 = j; }
			BVHDBG("  Slice2 size = %d\n", slice_size2);
			
			sort_bvh3(T, 0, slice_size2, 2);
			
			k = slice_size2;
			while(k > 0){
				int pack_size = 8;
				if(pack_size > k){ pack_size = k; }
				
				// Make a new node
				geom_bvh3d b = geom_bvh3d_make_blank();
				b->child[0] = T[0];
				b->tag = T[0]->tag;
				b->b[0] = T[0]->b[0];
				b->b[1] = T[0]->b[1];
				b->b[2] = T[0]->b[2];
				b->b[3] = T[0]->b[3];
				b->b[4] = T[0]->b[4];
				b->b[5] = T[0]->b[5];
				for(l = 1; l < pack_size; ++l){
					b->child[l] = T[l];
					if(T[l]->tag > b->tag){ b->tag = T[l]->tag; }
					if(T[l]->b[0] < b->b[0]){ b->b[0] = T[l]->b[0]; }
					if(T[l]->b[1] > b->b[1]){ b->b[1] = T[l]->b[1]; }
					if(T[l]->b[2] < b->b[2]){ b->b[2] = T[l]->b[2]; }
					if(T[l]->b[3] > b->b[3]){ b->b[3] = T[l]->b[3]; }
					if(T[l]->b[4] < b->b[4]){ b->b[4] = T[l]->b[4]; }
					if(T[l]->b[5] > b->b[5]){ b->b[5] = T[l]->b[5]; }
				}
				BVHDBG("  Making internal node %p, tag=%d, b=%f,%f,%f,%f,%f,%f\n", b, b->tag, b->b[0], b->b[1], b->b[2], b->b[3], b->b[4], b->b[5]);
				for(l = 0; l < pack_size; ++l){
					BVHDBG("   Child node %p, tag=%d, b=%f,%f,%f,%f,%f,%f\n", T[l], T[l]->tag, T[l]->b[0], T[l]->b[1], T[l]->b[2], T[l]->b[3], T[l]->b[4], T[l]->b[5]);
				}
				*B = b;
				T += pack_size;
				B++;
				(*n_)++;
				k -= pack_size;
			}
			
			j -= slice_size2;
		}
		
		i -= slice_size;
	}
	BVHDBG("Iteration of STR done; n=%d -> %d\n", n, *n_);
}
geom_bvh3d geom_bvh3d_new(unsigned int n, int (*shape_iterator)(double c[3], double h[3], int *tag, void *data), void *data){
	unsigned int i;
	
	geom_bvh3d *B = (geom_bvh3d*)malloc(sizeof(geom_bvh3d) * n);
	geom_bvh3d ret;
	for(i = 0; i < n; ++i){
		double c[3], h[3];
		B[i] = geom_bvh3d_make_blank();
		shape_iterator(c, h, &(B[i]->tag), data);
		B[i]->b[0] = c[0]-h[0];
		B[i]->b[1] = c[0]+h[0];
		B[i]->b[2] = c[1]-h[1];
		B[i]->b[3] = c[1]+h[1];
		B[i]->b[4] = c[2]-h[2];
		B[i]->b[5] = c[2]+h[2];
		BVHDBG("Making leaf %p, tag=%d, b=%f,%f,%f,%f,%f,%f\n", B[i], B[i]->tag, B[i]->b[0], B[i]->b[1], B[i]->b[2], B[i]->b[3], B[i]->b[4], B[i]->b[5]);
	}
	
	while(n > 1){
		geom_bvh3d_STR(&n, B);
	}
	ret = B[0];
	free(B);
	return ret;
}

void geom_bvh2d_destroy(geom_bvh2d bvh){
	int i;
	if(NULL == bvh){ return; }
	for(i = 0; i < 4; ++i){
		geom_bvh2d_destroy(bvh->child[i]);
	}
	free(bvh);
}
void geom_bvh3d_destroy(geom_bvh3d bvh){
	int i;
	if(NULL == bvh){ return; }
	for(i = 0; i < 8; ++i){
		geom_bvh3d_destroy(bvh->child[i]);
	}
	free(bvh);
}

int geom_bvh2d_query_pt(geom_bvh2d bvh, const double p[2], int (*query_func)(int tag, const double c[2], const double h[2], void *data), void *data){
	if(NULL == bvh){
		BVHDBG("bvh == NULL in BVH2_query_pt; this should never happen\n");
		return 1;
	}
	BVHDBG("Visiting node %p, tag=%d, b=%f,%f,%f,%f, int=%d\n", bvh, bvh->tag, bvh->b[0], bvh->b[1], bvh->b[2], bvh->b[3], (NULL != bvh->child[0]));
	if(!(bvh->b[0] <= p[0] && p[0] <= bvh->b[1] && bvh->b[2] <= p[1] && p[1] <= bvh->b[3])){
		// not in box
		return 1;
	}
	if(NULL == bvh->child[0]){
		double c[2], h[2];
		c[0] = 0.5*bvh->b[0] + 0.5*bvh->b[1];
		c[1] = 0.5*bvh->b[2] + 0.5*bvh->b[3];
		h[0] = 0.5*bvh->b[1] - 0.5*bvh->b[0];
		h[1] = 0.5*bvh->b[3] - 0.5*bvh->b[2];
		return query_func(bvh->tag, c, h, data);
	}else{
		int i;
		for(i = 0; i < 4; ++i){
			if(NULL == bvh->child[i]){ break; }
			if(0 == geom_bvh2d_query_pt(bvh->child[i], p, query_func, data)){ return 0; }
		}
		return 1;
	}
}
int geom_bvh3d_query_pt(geom_bvh3d bvh, const double p[3], int (*query_func)(int tag, const double c[3], const double h[3], void *data), void *data){
	if(NULL == bvh){
		BVHDBG("bvh == NULL in BVH3_query_pt; this should never happen\n");
		return 1;
	}
	BVHDBG("Visiting node %p, tag=%d, b=%f,%f,%f,%f,%f,%f, int=%d\n", bvh, bvh->tag, bvh->b[0], bvh->b[1], bvh->b[2], bvh->b[3], bvh->b[4], bvh->b[5], (NULL != bvh->child[0]));
	if(!(bvh->b[0] <= p[0] && p[0] <= bvh->b[1] && bvh->b[2] <= p[1] && p[1] <= bvh->b[3] && bvh->b[4] <= p[2] && p[2] <= bvh->b[5])){
		// not in box
		return 1;
	}
	if(NULL == bvh->child[0]){
		double c[3], h[3];
		c[0] = 0.5*bvh->b[0] + 0.5*bvh->b[1];
		c[1] = 0.5*bvh->b[2] + 0.5*bvh->b[3];
		c[2] = 0.5*bvh->b[4] + 0.5*bvh->b[5];
		h[0] = 0.5*bvh->b[1] - 0.5*bvh->b[0];
		h[1] = 0.5*bvh->b[3] - 0.5*bvh->b[2];
		h[2] = 0.5*bvh->b[5] - 0.5*bvh->b[4];
		return query_func(bvh->tag, c, h, data);
	}else{
		int i;
		for(i = 0; i < 8; ++i){
			if(NULL == bvh->child[i]){ break; }
			if(0 == geom_bvh3d_query_pt(bvh->child[i], p, query_func, data)){ return 0; }
		}
		return 1;
	}
}

int geom_bvh2d_query_box(geom_bvh2d bvh, const double c[2], const double h[2], int (*query_func)(int tag, const double c[2], const double h[2], void *data), void *data){
	if(NULL == bvh){
		BVHDBG("bvh == NULL in BVH2_query_box; this should never happen\n");
		return 1;
	}
	BVHDBG("Visiting node %p, tag=%d, b=%f,%f,%f,%f, int=%d\n", bvh, bvh->tag, bvh->b[0], bvh->b[1], bvh->b[2], bvh->b[3], (NULL != bvh->child[0]));
	if((bvh->b[1] < c[0]-h[0] || c[0]+h[0] < bvh->b[0]) && (bvh->b[3] < c[1]-h[1] || c[1]+h[1] < bvh->b[2])){
		// not in box
		return 1;
	}
	if(NULL == bvh->child[0]){
		double c[2], h[2];
		c[0] = 0.5*bvh->b[0] + 0.5*bvh->b[1];
		c[1] = 0.5*bvh->b[2] + 0.5*bvh->b[3];
		h[0] = 0.5*bvh->b[1] - 0.5*bvh->b[0];
		h[1] = 0.5*bvh->b[3] - 0.5*bvh->b[2];
		return query_func(bvh->tag, c, h, data);
	}else{
		int i;
		for(i = 0; i < 4; ++i){
			if(NULL == bvh->child[i]){ break; }
			if(0 == geom_bvh2d_query_box(bvh->child[i], c, h, query_func, data)){ return 0; }
		}
		return 1;
	}
}

int geom_bvh3d_query_box(geom_bvh3d bvh, const double c[3], const double h[3], int (*query_func)(int tag, const double c[3], const double h[3], void *data), void *data){
	if(NULL == bvh){
		BVHDBG("bvh == NULL in BVH3_query_box; this should never happen\n");
		return 1;
	}
	BVHDBG("Visiting node %p, tag=%d, b=%f,%f,%f,%f,%f,%f, int=%d\n", bvh, bvh->tag, bvh->b[0], bvh->b[1], bvh->b[2], bvh->b[3], bvh->b[4], bvh->b[5], (NULL != bvh->child[0]));
	if((bvh->b[1] < c[0]-h[0] || c[0]+h[0] < bvh->b[0]) && (bvh->b[3] < c[1]-h[1] || c[1]+h[1] < bvh->b[2]) && (bvh->b[5] < c[2]-h[2] || c[2]+h[2] < bvh->b[4])){
		// not in box
		return 1;
	}
	if(NULL == bvh->child[0]){
		double c[3], h[3];
		c[0] = 0.5*bvh->b[0] + 0.5*bvh->b[1];
		c[1] = 0.5*bvh->b[2] + 0.5*bvh->b[3];
		c[2] = 0.5*bvh->b[4] + 0.5*bvh->b[5];
		h[0] = 0.5*bvh->b[1] - 0.5*bvh->b[0];
		h[1] = 0.5*bvh->b[3] - 0.5*bvh->b[2];
		h[2] = 0.5*bvh->b[5] - 0.5*bvh->b[4];
		return query_func(bvh->tag, c, h, data);
	}else{
		int i;
		for(i = 0; i < 8; ++i){
			if(NULL == bvh->child[i]){ break; }
			if(0 == geom_bvh3d_query_box(bvh->child[i], c, h, query_func, data)){ return 0; }
		}
		return 1;
	}
}

int geom_bvh2d_traverse(geom_bvh2d bvh, int (*func)(int tag, const double c[2], const double h[2], int leaf, void *data), void *data){
	if(NULL == bvh){
		BVHDBG("bvh == NULL in BVH2_traverse; this should never happen\n");
		return 1;
	}
	BVHDBG("Visiting node %p, tag=%d, b=%f,%f,%f,%f, int=%d\n", bvh, bvh->tag, bvh->b[0], bvh->b[1], bvh->b[2], bvh->b[3], (NULL != bvh->child[0]));
	int i;
	double c[2], h[2];
	c[0] = 0.5*bvh->b[0] + 0.5*bvh->b[1];
	c[1] = 0.5*bvh->b[2] + 0.5*bvh->b[3];
	h[0] = 0.5*bvh->b[1] - 0.5*bvh->b[0];
	h[1] = 0.5*bvh->b[3] - 0.5*bvh->b[2];
	if(0 == func(bvh->tag, c, h, (NULL == bvh->child[0]), data)){ return 0; }
	for(i = 0; i < 4; ++i){
		if(NULL == bvh->child[i]){ break; }
		if(0 == geom_bvh2d_traverse(bvh->child[i], func, data)){ return 0; }
	}
	return 1;
}

int geom_bvh3d_traverse(geom_bvh3d bvh, int (*func)(int tag, const double c[3], const double h[3], int leaf, void *data), void *data){
	if(NULL == bvh){
		BVHDBG("bvh == NULL in BVH3_traverse; this should never happen\n");
		return 1;
	}
	BVHDBG("Visiting node %p, tag=%d, b=%f,%f,%f,%f,%f,%f, int=%d\n", bvh, bvh->tag, bvh->b[0], bvh->b[1], bvh->b[2], bvh->b[3], bvh->b[4], bvh->b[5], (NULL != bvh->child[0]));
	int i;
	double c[3], h[3];
	c[0] = 0.5*bvh->b[0] + 0.5*bvh->b[1];
	c[1] = 0.5*bvh->b[2] + 0.5*bvh->b[3];
	c[2] = 0.5*bvh->b[4] + 0.5*bvh->b[5];
	h[0] = 0.5*bvh->b[1] - 0.5*bvh->b[0];
	h[1] = 0.5*bvh->b[3] - 0.5*bvh->b[2];
	h[2] = 0.5*bvh->b[5] - 0.5*bvh->b[4];
	if(0 == func(bvh->tag, c, h, (NULL == bvh->child[0]), data)){ return 0; }
	for(i = 0; i < 8; ++i){
		if(NULL == bvh->child[i]){ break; }
		if(0 == geom_bvh3d_traverse(bvh->child[i], func, data)){ return 0; }
	}
	return 1;
}
