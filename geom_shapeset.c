#include <Cgeom/geom_shapeset.h>
#include <Cgeom/geom_bvh.h>
#include <stdlib.h>
#include <string.h>

#define geom_shapeset2d_threshold 32
#define geom_shapeset3d_threshold 32

//#define SS_DEBUG

#ifdef SS_DEBUG
# include <stdio.h>
# define SSDBG(...) fprintf(stderr, __VA_ARGS__)
#else
# define SSDBG(...)
#endif

typedef struct geom_shape2d_info_struct{
	geom_shape2d *s;
	geom_aabb2d box;
	unsigned int flags;
} geom_shape2d_info;

struct geom_shapeset2d_struct{
	unsigned int n;
	unsigned int n_alloc;
	geom_shape2d_info *info;
	
	geom_bvh2d bvh; // may not be used
	int use_bvh;
	
	int periodic;
	double lattice[4];
};

geom_shapeset2d geom_shapeset2d_new(){
	geom_shapeset2d ss = (geom_shapeset2d)malloc(sizeof(struct geom_shapeset2d_struct));
	ss->n = 0;
	ss->n_alloc = 0;
	ss->info = NULL;
	ss->bvh = NULL;
	ss->use_bvh = 0;
	ss->periodic = 0;
	return ss;
}

void geom_shapeset2d_destroy(geom_shapeset2d ss){
	if(NULL == ss){ return; }
	if(ss->use_bvh){
		geom_bvh2d_destroy(ss->bvh);
	}
	free(ss->info);
	free(ss);
}

int geom_shapeset2d_set_lattice(geom_shapeset2d ss, const double lattice[4]){
	if(NULL == ss){ return -1; }
	if(NULL == lattice){
		ss->periodic = 0;
	}else{
		ss->periodic = 1;
		ss->lattice[0] = lattice[0];
		ss->lattice[1] = lattice[1];
		ss->lattice[2] = lattice[2];
		ss->lattice[3] = lattice[3];
	}
	return 0;
}

struct shape2d_iter_data{
	unsigned int index;
	geom_shape2d_info *info;
};
static int shape2d_iter(double c[2], double h[2], int *tag, void *data){
	struct shape2d_iter_data *i = (struct shape2d_iter_data*)data;
	*tag = i->index; // the BVH's tag is actually the index into the vector info
	c[0] = i->info->box.c[0];
	c[1] = i->info->box.c[1];
	h[0] = i->info->box.h[0];
	h[1] = i->info->box.h[1];
	i->index++;
	i->info++;
	return 1;
}

int geom_shapeset2d_add(geom_shapeset2d ss, geom_shape2d *s){
	unsigned int i;
	if(NULL == ss){ return -1; }
	if(NULL == s){ return -2; }
	
	if(ss->n >= ss->n_alloc){
		if(0 == ss->n_alloc){
			ss->n_alloc = geom_shapeset2d_threshold;
		}else{
			ss->n_alloc *= 2;
		}
		ss->info = (geom_shape2d_info*)realloc(ss->info, sizeof(geom_shape2d_info) * ss->n_alloc);
	}
	i = ss->n;
	ss->info[i].s = s;
	ss->info[i].flags = 0;
	ss->info[i].flags |= geom_shape2d_get_aabb(s, &(ss->info[i].box)) ? GEOM_SHAPESET2D_FLAG_UNBOUNDED : 0;
	ss->n++;
	
	if(ss->use_bvh){
		ss->use_bvh = 0;
		geom_bvh2d_destroy(ss->bvh);
	}
	return i;
}

void geom_shapeset2d_finalize(geom_shapeset2d ss){
	if(NULL == ss || ss->use_bvh){ return; }
	struct shape2d_iter_data d;
	d.index = 0;
	d.info = ss->info;
	ss->bvh = geom_bvh2d_new(ss->n, &shape2d_iter, (void*)&d);
	ss->use_bvh = (NULL != ss->bvh);
}

unsigned int geom_shapeset2d_size(geom_shapeset2d ss){
	if(NULL == ss){ return 0; }
	return ss->n;
}

geom_shape2d* geom_shapeset2d_index(geom_shapeset2d ss, int index){
	if(NULL == ss){ return NULL; }
	if(index < 0 || index >= ss->n){ return NULL; }
	return ss->info[index].s;
}
int geom_shapeset2d_index_aabb(geom_shapeset2d ss, int index, geom_aabb2d *box){
	if(NULL == ss){ return -1; }
	if(index < 0 || index >= ss->n){ return -2; }
	if(NULL == box){ return -3; }
	memcpy(box, &ss->info[index].box, sizeof(geom_aabb2d));
	return ss->info[index].flags;
}

struct query_pt2d_data{
	geom_shape2d_info *info;
	double pc[2];
	int ibest;
};
static int query_pt2d(int tag, const double c[2], const double h[2], void *data){
	struct query_pt2d_data *d = (struct query_pt2d_data*)data;
	if(tag > d->ibest){
		if(geom_shape2d_contains(d->info[tag].s, d->pc)){	
			d->ibest = tag;
		}
	}
	return 1;
}

int geom_shapeset2d_query_pt(geom_shapeset2d ss, const double p[2]){
	unsigned int c, clim = 1;
	static const int off[] = {
		 0,  0, // must be first
		 1,  0,
		-1,  0,
		 0,  1,
		 0, -1,
		 1,  1,
		 1, -1,
		-1,  1,
		-1, -1
	};
	if(NULL == ss){ return -1; }
	if(NULL == p){ return -2; }
	
	if(ss->periodic){
		clim = 9;
	}
	struct query_pt2d_data d;
	d.info = ss->info;
	d.ibest = -1;
	for(c = 0; c < clim; ++c){
		d.pc[0] = p[0];
		d.pc[1] = p[1];
		if(0 != off[2*c+0]){
			d.pc[0] += (double)off[2*c+0] * ss->lattice[0];
			d.pc[1] += (double)off[2*c+0] * ss->lattice[1];
		}
		if(0 != off[2*c+1]){
			d.pc[0] += (double)off[2*c+1] * ss->lattice[2];
			d.pc[1] += (double)off[2*c+1] * ss->lattice[3];
		}
		if(ss->use_bvh){
			geom_bvh2d_query_pt(ss->bvh, d.pc, &query_pt2d, &d);
		}else{
			int i;
			for(i = 0; i < ss->n; ++i){
				if(i <= d.ibest){ continue; } // skip anything less the current best
				if(GEOM_SHAPESET2D_FLAG_UNBOUNDED & ss->info[i].flags){
					if(geom_shape2d_contains(ss->info[i].s, d.pc)){
						d.ibest = i;
					}
				}else{
					if(geom_aabb2d_contains(&(ss->info[i].box), d.pc)){
						if(geom_shape2d_contains(ss->info[i].s, d.pc)){
							d.ibest = i;
						}
					}
				}
			}
		}
	}
	return d.ibest;
}
int geom_shapeset2d_foreach(
	geom_shapeset2d ss,
	int (*func)(geom_shape2d *s, const geom_aabb2d *box, unsigned int flags, void *data),
	void *data
){
	if(NULL == ss){ return -1; }
	if(NULL == func){ return -2; }
	unsigned int i;
	for(i = 0; i < ss->n; ++i){
		if(!func(ss->info[i].s, &ss->info[i].box, ss->info[i].flags, data)){ break; }
	}
	return i;
}

























typedef struct geom_shape3d_info_struct{
	geom_shape3d *s;
	geom_aabb3d box;
	unsigned int flags;
} geom_shape3d_info;

struct geom_shapeset3d_struct{
	unsigned int n;
	unsigned int n_alloc;
	geom_shape3d_info *info;
	
	geom_bvh3d bvh; // may not be used
	int use_bvh;
	
	int periodic;
	double lattice[9];
};

geom_shapeset3d geom_shapeset3d_new(){
	geom_shapeset3d ss = (geom_shapeset3d)malloc(sizeof(struct geom_shapeset3d_struct));
	ss->n = 0;
	ss->n_alloc = 0;
	ss->info = NULL;
	ss->bvh = NULL;
	ss->use_bvh = 0;
	ss->periodic = 0;
	return ss;
}

void geom_shapeset3d_destroy(geom_shapeset3d ss){
	if(NULL == ss){ return; }
	if(ss->use_bvh){
		geom_bvh3d_destroy(ss->bvh);
	}
	free(ss->info);
	free(ss);
}

int geom_shapeset3d_set_lattice(geom_shapeset3d ss, const double lattice[9]){
	if(NULL == ss){ return -1; }
	if(NULL == lattice){
		ss->periodic = 0;
	}else{
		ss->periodic = 1;
		ss->lattice[0] = lattice[0];
		ss->lattice[1] = lattice[1];
		ss->lattice[2] = lattice[2];
		ss->lattice[3] = lattice[3];
		ss->lattice[4] = lattice[4];
		ss->lattice[5] = lattice[5];
		ss->lattice[6] = lattice[6];
		ss->lattice[7] = lattice[7];
		ss->lattice[8] = lattice[8];
	}
	return 0;
}

struct shape3d_iter_data{
	unsigned int index;
	geom_shape3d_info *info;
};
static int shape3d_iter(double c[2], double h[2], int *tag, void *data){
	struct shape3d_iter_data *i = (struct shape3d_iter_data*)data;
	*tag = i->index; // the BVH's tag is actually the index into the vector info
	c[0] = i->info->box.c[0];
	c[1] = i->info->box.c[1];
	c[2] = i->info->box.c[2];
	h[0] = i->info->box.h[0];
	h[1] = i->info->box.h[1];
	h[2] = i->info->box.h[2];
	i->index++;
	i->info++;
	return 1;
}

int geom_shapeset3d_add(geom_shapeset3d ss, geom_shape3d *s){
	unsigned int i;
	if(NULL == ss){ return -1; }
	if(NULL == s){ return -2; }
	
	if(ss->n >= ss->n_alloc){
		if(0 == ss->n_alloc){
			ss->n_alloc = geom_shapeset3d_threshold;
		}else{
			ss->n_alloc *= 2;
		}
		ss->info = (geom_shape3d_info*)realloc(ss->info, sizeof(geom_shape3d_info) * ss->n_alloc);
	}
	i = ss->n;
	ss->info[i].s = s;
	ss->info[i].flags = 0;
	ss->info[i].flags |= geom_shape3d_get_aabb(s, &(ss->info[i].box)) ? GEOM_SHAPESET3D_FLAG_UNBOUNDED : 0;
	ss->n++;
	
	if(ss->use_bvh){
		ss->use_bvh = 0;
		geom_bvh3d_destroy(ss->bvh);
	}
	return i;
}

void geom_shapeset3d_finalize(geom_shapeset3d ss){
	if(NULL == ss || ss->use_bvh){ return; }
	struct shape3d_iter_data d;
	d.index = 0;
	d.info = ss->info;
	ss->bvh = geom_bvh3d_new(ss->n, &shape3d_iter, (void*)&d);
	ss->use_bvh = (NULL != ss->bvh);
}

unsigned int geom_shapeset3d_size(geom_shapeset3d ss){
	if(NULL == ss){ return 0; }
	return ss->n;
}

geom_shape3d* geom_shapeset3d_index(geom_shapeset3d ss, int index){
	if(NULL == ss){ return NULL; }
	if(index < 0 || index >= ss->n){ return NULL; }
	return ss->info[index].s;
}
int geom_shapeset3d_index_aabb(geom_shapeset3d ss, int index, geom_aabb3d *box){
	if(NULL == ss){ return -1; }
	if(index < 0 || index >= ss->n){ return -2; }
	if(NULL == box){ return -3; }
	memcpy(box, &ss->info[index].box, sizeof(geom_aabb3d));
	return ss->info[index].flags;
}



static int query_pt3d(int tag, const double c[3], const double h[3], void *data){
	int *bbest = (int*)data;
	if(tag > *bbest){ *bbest = tag; }
	return 1;
}

int geom_shapeset3d_query_pt(geom_shapeset3d ss, const double p[3]){
	unsigned int c, clim = 1;
	static const int off[] = {
		 0,  0,  0, // must be first
		 0,  0, -1,
		 1,  0, -1,
		-1,  0, -1,
		 0,  1, -1,
		 0, -1, -1,
		 1,  1, -1,
		 1, -1, -1,
		-1,  1, -1,
		-1, -1, -1,
		 1,  0,  0,
		-1,  0,  0,
		 0,  1,  0,
		 0, -1,  0,
		 1,  1,  0,
		 1, -1,  0,
		-1,  1,  0,
		-1, -1,  0,
		 0,  0,  1,
		 1,  0,  1,
		-1,  0,  1,
		 0,  1,  1,
		 0, -1,  1,
		 1,  1,  1,
		 1, -1,  1,
		-1,  1,  1,
		-1, -1,  1
	};
	if(NULL == ss){ return -1; }
	if(NULL == p){ return -2; }

	if(ss->periodic){
		clim = 27;
	}
	int ibest = -1;
	for(c = 0; c < clim; ++c){
		const double pc[3] = {
			p[0] + (double)off[3*c+0] * ss->lattice[0] + (double)off[3*c+1] * ss->lattice[3] + (double)off[3*c+2] * ss->lattice[6],
			p[1] + (double)off[3*c+0] * ss->lattice[1] + (double)off[3*c+1] * ss->lattice[4] + (double)off[3*c+2] * ss->lattice[7],
			p[2] + (double)off[3*c+0] * ss->lattice[2] + (double)off[3*c+1] * ss->lattice[5] + (double)off[3*c+2] * ss->lattice[8]
		};
		if(ss->use_bvh){
			int bbest = -1;
			geom_bvh3d_query_pt(ss->bvh, pc, &query_pt3d, &bbest);
			if(bbest > ibest){
				if(geom_shape3d_contains(ss->info[bbest].s, pc)){
					ibest = bbest;
				}
			}
		}else{
			int i;
			for(i = 0; i < ss->n; ++i){
				if(i <= ibest){ continue; } // skip anything less the current best
				if(GEOM_SHAPESET3D_FLAG_UNBOUNDED & ss->info[i].flags){
					if(geom_shape3d_contains(ss->info[i].s, pc)){
						ibest = i;
					}
				}else{
					if(geom_aabb3d_contains(&(ss->info[i].box), pc)){
						if(geom_shape3d_contains(ss->info[i].s, pc)){
							ibest = i;
						}
					}
				}
			}
		}
	}
	return 0;
}
int geom_shapeset3d_foreach(
	geom_shapeset3d ss,
	int (*func)(geom_shape3d *s, const geom_aabb3d *box, unsigned int flags, void *data),
	void *data
){
	if(NULL == ss){ return -1; }
	if(NULL == func){ return -2; }
	unsigned int i;
	for(i = 0; i < ss->n; ++i){
		if(!func(ss->info[i].s, &ss->info[i].box, ss->info[i].flags, data)){ break; }
	}
	return i;
}
