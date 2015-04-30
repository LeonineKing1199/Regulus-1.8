/*
 * structures.h
 *
 *  Created on: Apr 9, 2015
 *      Author: christian
 */

#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include <cassert>
#include <iostream>

#include <thrust/device_ptr.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/find.h>
#include <thrust/distance.h>

const int bpg = 512,
		  tpb = 256;

// float or double
#define REAL_TYPE_FP32

#ifdef REAL_TYPE_FP32
typedef float RealType;
typedef RealType real;
#else
typedef double RealType;
typedef RealType real;
#endif


__forceinline__
__device__
int currThreadID(void)
{
	const int tid = blockIdx.x * blockDim.x + threadIdx.x;
	return tid;
}

__forceinline__
__device__
int grid_size(void)
{
	const int gs = blockDim.x * gridDim.x;
	return gs;
}

template<typename T>
T* host_malloc(const int num_elements)
{
	T *allocation = 0;
	const size_t num_bytes = num_elements * sizeof(T);
	cudaMallocHost(&allocation, num_bytes);
	return allocation;
}

template<typename T>
void host_free(T*& allocation)
{
	cudaFreeHost(allocation);
	allocation = 0;
}

template<typename T>
T* device_malloc(const int num_elements)
{
	T *allocation = 0;
	const size_t num_bytes = num_elements * sizeof(T);
	cudaMalloc(&allocation, num_bytes);
	return allocation;
}

template<typename T>
void device_free(T*& allocation)
{
	cudaFree(allocation);
	allocation = 0;
}

__global__
void find_boundaries
(
const int  num_keys,
const int  num_buckets,
const int * __restrict__ which_bucket,
int * __restrict__ bucket_starts
);

/*
 * Basic tetrahedral structure
 */

struct tetra
{
	int v[4];
};

/*
 * Simple point structure
 */

struct point
{
	real c[3];
};



/*
 * Declared for partial vectorized loads of
 * point structure (taken from gDel3d()
 * See license in KerShewchuk.h
 */

#ifdef REAL_TYPE_FP32

struct __align__(8) real2
{
    float a, b;
};

#else

struct __align__(16) real2
{
    double a, b;
};

#endif

/*
 * Credit for these loading routines goes to the gDel3D authors
 * See the license in KerShewchuk.h
 */
__forceinline__
__device__
point load_point(const point * __restrict__ point_set, const int idx)
{
    if (idx % 2 == 0)
    {
        const real2 ptxy = ((real2* ) point_set)[idx * 3 / 2];
        const point pt = { ptxy.a,
        		           ptxy.b,
        		           point_set[idx].c[2] };

        return pt;
    }
    else
    {
        const real2 ptyz = ((real2* ) point_set)[idx * 3 / 2 + 1];
        const point pt = { point_set[ idx ].c[0],
        				   ptyz.a,
        				   ptyz.b };

        return pt;
    }
}

__forceinline__
__device__
tetra load_tetra(const tetra * __restrict__ mesh, const int idx)
{
	const int4 tmp = ((int4* ) mesh)[idx];
	const tetra t = { tmp.x, tmp.y, tmp.z, tmp.w };
	return t;
}

__forceinline__
__device__
void store_tetra(tetra *mesh, const int idx, const tetra t)
{
	int4 tmp =  { t.v[0], t.v[1], t.v[2], t.v[3] };
	((int4* ) mesh)[idx] = tmp;
}

/*
 * User-oriented mesh class
 */

struct mesh
{
private :
	point *user_points;
	tetra *user_mesh;

	point *device_points;
	tetra *device_mesh;

	int grid_length;
	int root_edge_length;
	int num_grid_points;
	int num_root_points;
	int num_tetra;

	void build_domain(void);

public :
	mesh(const int gl);
	~mesh(void);

	point* get_user_points(void) const;
	tetra* get_user_mesh(void) const;

	void print_user_points(void) const;

	void triangulate(void);
};

#endif /* STRUCTURES_H_ */
