/*
 * triangulate.cu
 *
 *  Created on: Apr 9, 2015
 *      Author: christian
 */




#include "structures.h"
#include "KerShewchuk.h"

/*
 * We must also make sure that there are no duplicates as well!
 */

void assert_no_duplicates
(

)
{
	//int *face_hashes = device_malloc<int>
}

/*
 * Make a routine to ensure there are no flat tetrahedra
 */

__global__
void assert_no_flat_tetra
(
const tetra * __restrict__ mesh,
const point * __restrict__ point_set,
const float * __restrict__ predConsts,
const int num_tetra
)
{
	const int thread_id = currThreadID();

	for (int tid = thread_id; tid < num_tetra; tid += grid_size())
	{
		const tetra t = load_tetra(mesh, tid);

		const point a = load_point(point_set, t.v[0]);
		const point b = load_point(point_set, t.v[1]);
		const point c = load_point(point_set, t.v[2]);
		const point d = load_point(point_set, t.v[3]);

		const real ort = orient3dFast(predConsts, (float* ) &a, (float* ) &b, (float* ) &c, (float* ) &d);

		assert(ort >= 0);
	}
}

/*
 * Type of debug routine to visually examine data
 */

__global__
void array_dump
(
const int array_capacity,
const int * __restrict__ pa,
const int * __restrict__ ta,
const int * __restrict__ la,
const int * __restrict__ fs
)
{
	const int thread_id = currThreadID();

	for (int tid = thread_id; tid < array_capacity; tid += grid_size())
	{
		printf("Tuple[%d] : { pa = %d, ta = %d, la = %d, fs = %d }\n", tid, pa[tid], ta[tid], la[tid], fs[tid]);
	}
}

/*
 * Light kernel to clear all inserted points as well
 * from the associated arrays
 */

template<typename T>
struct tuple_comp
{
    __host__ __device__
    bool operator()(const thrust::tuple<T, T, T, T, T> t,
                    const thrust::tuple<T, T, T, T, T> v)
    {
        return ((unsigned& ) thrust::get<0>(t)) < ((unsigned& ) thrust::get<0>(v));
    }
};

__global__
void clear_inserted
(
const int array_size,
int * __restrict__ pa,
int * __restrict__ ta,
int * __restrict__ la,
int * __restrict__ fs,
int * __restrict__ nm
)
{
	const int thread_id = currThreadID();

	for (int tid = thread_id; tid < array_size; tid += grid_size())
	{
		if (nm[tid] == 1)
		{
			pa[tid] = -1;
			ta[tid] = -1;
			la[tid] = -1;
			fs[tid] = -1;
			nm[tid] = 0;
		}
	}
}

/*
 * We now redistribute the points...
 */

__global__
void redistribute_points
(
const int array_size,
const int num_tetra,
const int * __restrict__ nm,
const int * __restrict__ nominated_tetra,
const tetra * __restrict__ mesh,
const int * __restrict__ fract_offsets,
const int * __restrict__ assoc_offsets,
const point * __restrict__ point_set,
const float * __restrict__ predConsts,
int * __restrict__ pa,
int * __restrict__ ta,
int * __restrict__ la,
int * __restrict__ fs
)
{
	const int thread_id = currThreadID();

	for (int tid = thread_id; tid < array_size; tid += grid_size())
	{
		if (nm[tid] == 0)
		{
			const int ta_id = ta[tid];

			const int tuple_id = nominated_tetra[ta_id];

			if (tuple_id != -1)
			{
				assert(tuple_id >= 0);

				// read in tetrahedra...
				tetra t[4];
				int back = 0;

				t[0] = load_tetra(mesh, ta[tuple_id]);
				++back;

				// get size of fracture...
				const int fract_size = fs[tuple_id];
				const int fract_offset = num_tetra + (tuple_id > 0 ? fract_offsets[tuple_id - 1] : 0);

				// read in remaining...
				for (int i = 1; i < fract_size; ++i, ++back)
				{
					const int idx = fract_offset + i - 1;
					t[i] = load_tetra(mesh, idx);
				}

				assert(back == fract_size);

				// load in points...
				const int pa_id = pa[tid];
				const point p = load_point(point_set, pa_id);

				point a, b, c, d;

				// now begin testing...
				const int assoc_offset = (tid > 0 ? assoc_offsets[tid -1] : 0) + array_size;

				for (int i = 0; i < back; ++i)
				{
					const tetra tmp = t[i];

					a = load_point(point_set, tmp.v[0]);
					b = load_point(point_set, tmp.v[1]);
					c = load_point(point_set, tmp.v[2]);
					d = load_point(point_set, tmp.v[3]);

					// calculate orientation values
					const real ort0 = orient3dFast(predConsts, (float* ) &d, (float* ) &c, (float* ) &b, (float* ) &p);
					const real ort1 = orient3dFast(predConsts, (float* ) &a, (float* ) &c, (float* ) &d, (float* ) &p);
					const real ort2 = orient3dFast(predConsts, (float* ) &a, (float* ) &d, (float* ) &b, (float* ) &p);
					const real ort3 = orient3dFast(predConsts, (float* ) &a, (float* ) &b, (float* ) &c, (float* ) &p);

//#define print
#ifdef print
					printf("Testing tetrahedron %d:\n"
						   "a: %.00f %.00f %.00f\n"
						   "b: %.00f %.00f %.00f\n"
						   "c: %.00f %.00f %.00f\n"
						   "d: %.00f %.00f %.00f\n\n"
						   "Against point %d:\n"
						   "%.00f %.00f %.00f\n\n"
						   "Orientations are :\n"
						   "%.00f %.00f %.00f %.00f\n\n",
						   i,
						   a.c[0], a.c[1], a.c[2],
						   b.c[0], b.c[1], b.c[2],
						   c.c[0], c.c[1], c.c[2],
						   d.c[0], d.c[1], d.c[2],
						   pa_id,
						   p.c[0], p.c[1], p.c[2],
						   ort0, ort1, ort2, ort3);
#undef print
#endif

					// if any of the orientations are negative,
					// we can continue
					if (ort0 < 0 || ort1 < 0 || ort2 < 0 || ort3 < 0)
						continue;

					assert(ort0 + ort1 + ort2 + ort3 != 0);

					// assign location value
					int loc = 0;
					loc |= ((ort0 > 0) << 0);
					loc |= ((ort1 > 0) << 1);
					loc |= ((ort2 > 0) << 2);
					loc |= ((ort3 > 0) << 3);

					// get fracture size (number of set bits)
					const int fract_size = __popc(loc);

					assert(fract_size > 0);

					// write back to memory
					const int idx = assoc_offset + i;

					assert(idx != tid);
					assert(pa[idx] == -1 && ta[idx] == -1 && la[idx] == -1 && fs[idx] == -1);

					pa[idx] = pa_id;
					ta[idx] = (i > 0 ? fract_offset + i - 1 : ta_id);
					la[idx] = loc;
					fs[idx] = fract_size;

					//assert(pa[idx] != -1 && ta[idx] != -1 && la[idx] != -1 && fs[idx] != -1);
				}

				// now we must clear the previous association data
				pa[tid] = -1;
				ta[tid] = -1;
				la[tid] = -1;
				fs[tid] = -1;
			}
		}
	}
}

/*
 * We need a way to get the potential number of
 * new associations as we redistribute the point sets,
 * i.e. this is one step to redistribution of points
 */

__global__
void get_association_offsets
(
const int array_size,
const int * __restrict__ nm,
const int * __restrict__ ta,
const int * __restrict__ nominated_tetra,
const int * __restrict__ fs,
int * __restrict__ assoc_offsets
)
{
	const int thread_id = currThreadID();

	for (int tid = thread_id; tid < array_size; tid += grid_size())
	{
		// initialization value
		assoc_offsets[tid] = 0;

		// want all tuples NOT nominated...
		if (nm[tid] == 0)
		{
			// get id of current tetrahedron
			// in the tuple
			const int ta_id = ta[tid];

			// if this tetrahedron was fractured...
			const int tuple_id = nominated_tetra[ta_id];
			if (tuple_id != -1)
			{
				const int fract_size = fs[tuple_id];
				assoc_offsets[tid] = fract_size;

//#define print
#ifdef print
				printf("assoc_offsets[%d] = %d\n", tid, assoc_offsets[tid]);
#undef print
#endif
			}
		}
	}
}

/*
 * Routine to fracture some tetrahedra!
 */

__global__
void fracture
(
const int array_size,
const int * __restrict__ nm,
const int * __restrict__ ta,
const int * __restrict__ pa,
const int * __restrict__ la,
const int * __restrict__ fract_offsets,
const int num_tetra,
tetra * __restrict__ mesh
)
{
	const int thread_id = currThreadID();

	for (int tid = thread_id; tid < array_size; tid += grid_size())
	{
		// if the tuple is nominated...
		if (nm[tid] == 1)
		{
			// read in relevant indices
			const int ta_id = ta[tid];
			const int pa_id = pa[tid];
			const int loc = la[tid];

			//printf("Going to insert point : %d (la = %d)\n", pa_id, loc);

			// allocate relevant storage
			tetra new_tetra[4];
			int back = 0;

			// read in tetrahedron
			const tetra t = load_tetra(mesh, ta_id);

			// face indexing...
			const int face_id[4][3] = { { 3, 2, 1 },
										{ 0, 2, 3 },
										{ 0, 3, 1 },
									    { 0, 1, 2 } };

			// get stack-local copies of new tetrahedra
			for (int i = 0; i < 4; ++i)
			{
				// if the bit is set...
				if (loc & (1 << i))
				{
					const tetra tmp = { t.v[face_id[i][0]],
										t.v[face_id[i][1]],
										t.v[face_id[i][2]],
										pa_id };

					new_tetra[back] = tmp;
					++back;
				}
			}

			// write back to mesh
			store_tetra(mesh, ta_id, new_tetra[0]);

			const int fract_begin = num_tetra + (tid > 0 ? fract_offsets[tid - 1] : 0);
			const int fract_end = num_tetra + fract_offsets[tid];

			for (int i = fract_begin, j = 1; i < fract_end; ++i, ++j)
			{
				assert(j < back); // no out-of-bounds
				store_tetra(mesh, i, new_tetra[j]);
			}
		}
	}
}

/*
 * Routine to gather offsets for writing into
 * the mesh array
 */

template<typename T>
struct dot : public thrust::unary_function<thrust::tuple<T, T>, T>
{
    __host__ __device__
    T operator()(thrust::tuple<T, T> v1)
    {
        return (thrust::get<0>(v1) - 1) * thrust::get<1>(v1);
    }
};

int* gather_offsets
(
const int array_size,
int * __restrict__ fs,
int * __restrict__ nm
)
{
	// want to compute prefix sum
	int *fract_offsets = device_malloc<int>(array_size);
	const int init = 0;
	cudaMemset(fract_offsets, init, array_size * sizeof(*fract_offsets));

	// fracture offsets are determined by the prefix sum
	// of all the nominated fracture sites in the
	// main tuple (ta, pa, la, fs)
	thrust::inclusive_scan(
			// begin iterator
			thrust::make_transform_iterator(
					thrust::make_zip_iterator(
							thrust::make_tuple(
									thrust::device_ptr<int>(fs),
									thrust::device_ptr<int>(nm))),
					dot<int>()),
			// end iterator
			thrust::make_transform_iterator(
					thrust::make_zip_iterator(
							thrust::make_tuple(
									thrust::device_ptr<int>(fs + array_size),
									thrust::device_ptr<int>(nm + array_size))),
					dot<int>()),
			// write iterator
			thrust::device_ptr<int>(fract_offsets));

	return fract_offsets;
}

/*
 * This routine uses a series of locks to
 * determine the intersection of bucket contents
 * Currently the time it takes a thread to get
 * to the atomicCAS() call is what determines which
 * buckets "win". I like the idea of using the warp
 * scheduler to determine the set of non-intersecting
 * buckets even though it introduces a stark lack of
 * control over the code.
 */

__global__
void resolve_conflicts
(
const int num_buckets,
const int * __restrict__ bucket_starts,
const int * __restrict__ ta,
int * __restrict__ locked_tetra,
int * __restrict__ nominated_tetra,
int * __restrict__ nm
)
{
	const int thread_id = currThreadID();

	// one thread per bucket
	for (int tid = thread_id; tid < num_buckets; tid += grid_size())
	{
		// iterate the bucket
		const int begin = (tid > 0 ? bucket_starts[tid - 1] : 0);
		const int end = bucket_starts[tid];

		for (int i = begin; i < end; ++i)
		{
			// get id of tetrahedron
			const int ta_id = ta[i];

			// check the lock...
			const int old = atomicCAS(locked_tetra + ta_id, 0, 1);

			// use this to filter out "bad" threads (i.e. the bucket
			// contents intersect)
			if (old == 1)
				return;
		}

		// now, these are the only threads alive
		for (int i = begin; i < end; ++i)
		{
			const int ta_id = ta[i];

			assert(nominated_tetra[ta_id] == -1);

			// we store the location of the nominated tetrahedron
			// in ta and index it by the tetrahedron id
			// nominated_tet[0] = 17; <-- tetrahedron 0 has assocation
			// array data stored at index 17 (pa[17], la[17], fs[17], ta[17] == 0)
			nominated_tetra[ta_id] = i;

			// then write to special nomination array
			nm[i] = 1;

//#define print
#ifdef print
			printf("nominated_tets[%d] = %d\n", ta_id, nominated_tetra[ta_id]);
#undef print
#endif
		}
	}
}

/*
 * Slightly more complex initialization routine for the
 * location associations and fracture sizes
 */

__global__
void init_la_and_fs
(
const int num_grid_points,
const int * __restrict__ ta,
const int * __restrict__ pa,
const point * __restrict__ point_set,
const tetra * __restrict__ mesh,
const float * __restrict__ predConsts,
int * __restrict__ la,
int * __restrict__ fs
)
{
	// get current thread id
	const int thread_id = currThreadID();

	// stride through data with grid step
	for (int tid = thread_id; tid < num_grid_points; tid += grid_size())
	{
		// read in tet and pt id's
		const int ta_id = ta[tid];
		const int pa_id = pa[tid];

		// load in the pt and tetrahedron
		const tetra t = load_tetra(mesh, ta_id);
		const point p = load_point(point_set, pa_id);

		// load in the 4 vertices as point structures
		const point a = load_point(point_set, t.v[0]);
		const point b = load_point(point_set, t.v[1]);
		const point c = load_point(point_set, t.v[2]);
		const point d = load_point(point_set, t.v[3]);

		// calculate orientation values
		const real ort0 = orient3dFast(predConsts, (float* ) &d, (float* ) &c, (float* ) &b, (float* ) &p);
		const real ort1 = orient3dFast(predConsts, (float* ) &a, (float* ) &c, (float* ) &d, (float* ) &p);
		const real ort2 = orient3dFast(predConsts, (float* ) &a, (float* ) &d, (float* ) &b, (float* ) &p);
		const real ort3 = orient3dFast(predConsts, (float* ) &a, (float* ) &b, (float* ) &c, (float* ) &p);

		assert(ort0 >= 0 && ort1 >= 0 && ort2 >= 0 && ort3 >= 0);

		// finally, calculate the value of la[i]
		int loc = 0;

		loc |= ((ort0 > 0) << 0);
		loc |= ((ort1 > 0) << 1);
		loc |= ((ort2 > 0) << 2);
		loc |= ((ort3 > 0) << 3);

		// fracture size is simply the number of the set bits in la[i]
		// __popc() = integer intrinsic that returns number of set bits
		const int fract_size = __popc(loc);

		// assignment
		la[tid] = loc;
		fs[tid] = fract_size;

//#define print
#ifdef print
			printf("Testing tetrahedron :\n"
				   "%.00f %.00f %.00f\n"
				   "%.00f %.00f %.00f\n"
				   "%.00f %.00f %.00f\n"
				   "%.00f %.00f %.00f\n\n"
				   "against point :\n"
				   "%.00f %.00f %.00f\n\n"
				   "Orientations are :\n"
				   "%.00f, %.00f, %.00f, %.00f\n\n"
				   "la = %d, fs = %d\n\n",
				   a.c[0], a.c[1], a.c[2],
				   b.c[0], b.c[1], b.c[2],
				   c.c[0], c.c[1], c.c[2],
				   d.c[0], d.c[1], d.c[2],
				   p.c[0], p.c[1], p.c[2],
				   ort0, ort1, ort2, ort3,
				   loc, fract_size);
#undef print
#endif
	}
}

/*
 * Routine used to initialize the tetrahedron/point relations (deprecated)
 */

__global__
void init_ta_and_pa
(
int * __restrict__ ta,
int * __restrict__ pa,
const int num_grid_points
)
{
	const int thread_id = currThreadID();

	for (int tid = thread_id; tid < num_grid_points; tid += grid_size())
	{
		// All points are initially in the root
		ta[tid] = 0;

		// pa is initially the sequence of grid points
		pa[tid] = tid;
	}
}

/*
 * Pair of routines to initialize the arrays
 */

__global__
void init
(
const int num_tetra,
int * __restrict__ locked_tetra,
int * __restrict__ nominated_tetra
)
{
	const int thread_id = currThreadID();

	for (int tid = thread_id; tid < num_tetra; tid += grid_size())
	{
		locked_tetra[tid] = 0; // initially unlocked
		nominated_tetra[tid] = -1; // impossible value
	}
}

__global__
void init
(
const int array_capacity,
const int num_grid_points,
int * __restrict__ pa,
int * __restrict__ ta,
int * __restrict__ la,
int * __restrict__ fs,
int * __restrict__ nm
)
{
	const int thread_id = currThreadID();

	for (int tid = thread_id; tid < array_capacity; tid += grid_size())
	{
		pa[tid] = (tid < num_grid_points ? tid : -1);
		ta[tid] = (tid < num_grid_points ? 0 : -1);
		la[tid] = -1;
		fs[tid] = -1;
		nm[tid] = 0;
	}
}

/*
 * Core triangulation routine.
 * This is the true point of reading. From here on, most
 * functions are defined in a stack-like manner so the first
 * is on the bottom and the top is the last "big" function used.
 */

void mesh::triangulate(void)
{
	// Build main association arrays
	// ta = tetrahedron association
	// pa = point association
	// la = location association
	// fs = fracture size
	// Point pa[i] is in or on tetrahedron ta[i] with
	// a location encoded by la[i] and a potential
	// fracture size of fs[i] (i.e. 1-to-4 flip,
	// 1-to-3 flip, or 1-to-2 flip)
	int array_size = num_grid_points;
	const int array_capacity = 8 * array_size;

	int *ta = device_malloc<int>(array_capacity);
	int *pa = device_malloc<int>(array_capacity);
	int *la = device_malloc<int>(array_capacity);
	int *fs = device_malloc<int>(array_capacity);

	// Simple array to store nomination info
	int *nm = device_malloc<int>(array_capacity);

	// Init predicate stuff (the gFlip way)
	float *predConsts = cuNew< RealType >( DPredicateBoundNum );
    kerInitPredicate<<< 1, 1 >>>(predConsts);

    // Initialize the arrays
    init<<<bpg, tpb>>>
    (
    array_capacity,
    num_grid_points,
    pa,
    ta,
    la,
    fs,
    nm
    );

    // We then initialize fs and la which is a tad more
    // computationally expensive
    init_la_and_fs<<<bpg, tpb>>>
	(
	num_grid_points,
	ta,
	pa,
	device_points,
	device_mesh,
	predConsts,
	la,
	fs
	);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);

	// Begin triangulation loop
	while (array_size != 0)
	{
		// Allocate memory for locks and storage (explained later
		// in kernel)
		int *locked_tetra = device_malloc<int>(num_tetra);
		int *nominated_tetra = device_malloc<int>(num_tetra);

	    init<<<bpg, tpb>>>
	    (
	    num_tetra,
	    locked_tetra,
	    nominated_tetra
	    );

	    // Sorted data is required for bucket hashing to work
	    // Could possibly overlap this with previous two kernels
	    // but payoff might not be large
	    thrust::sort_by_key(thrust::device_ptr<int>(pa),
	    		            thrust::device_ptr<int>(pa + array_size),
	    		            thrust::make_zip_iterator(thrust::make_tuple(thrust::device_ptr<int>(ta),
	    		            		                                     thrust::device_ptr<int>(la),
	    		            		                                     thrust::device_ptr<int>(fs))));

	    cudaDeviceSynchronize();

	    // We now need the number of buckets. In this case, we say
	    // the total number of buckets is the largest value in pa
	    // + 1 (i.e. maximum value determines number of buckets)
	    // This part of the routine probably a huge performance bottleneck
	    int num_buckets = 0;
	    cudaMemcpy(&num_buckets, pa + array_size - 1, sizeof(int), cudaMemcpyDeviceToHost);
	    ++num_buckets;

	    // Build array to store start of buckets
	    int *bucket_starts = device_malloc<int>(num_buckets);

	    // Get start of each bucket in pa
	    find_boundaries<<<bpg, tpb>>>
	    (
	    array_size, // number of keys
	    num_buckets,
	    pa,
	    bucket_starts
	    );

	    // Nominate a series of buckets
	    resolve_conflicts<<<bpg, tpb>>>
	    (
	    num_buckets,
	    bucket_starts,
	    ta,
	    locked_tetra,
	    nominated_tetra,
	    nm
	    );

	    cudaDeviceSynchronize();

	    // Fracture tetrahedra
	    int *fract_offsets = gather_offsets(array_size, fs, nm);
	    int *assoc_offsets = device_malloc<int>(array_size);

	    fracture<<<bpg, tpb>>>
	    (
	    array_size,
	    nm,
	    ta,
	    pa,
	    la,
	    fract_offsets,
	    num_tetra,
	    device_mesh
	    );

	    // Get potential new association sizes for
	    // redistribution of points
	    get_association_offsets<<<bpg, tpb>>>
	    (
	    array_size,
	    nm,
	    ta,
	    nominated_tetra,
	    fs,
	    assoc_offsets
	    );

	    cudaDeviceSynchronize();

	    // Use prefix sum to get association write indices
	    thrust::inclusive_scan(thrust::device_ptr<int>(assoc_offsets),
	    				       thrust::device_ptr<int>(assoc_offsets + array_size),
	    				       thrust::device_ptr<int>(assoc_offsets));

	    // Attempt to redistribute the points correctly...
	    redistribute_points<<<bpg, tpb>>>
	    (
	    array_size,
	    num_tetra,
	    nm,
	    nominated_tetra,
	    device_mesh,
	    fract_offsets,
	    assoc_offsets,
	    device_points,
	    predConsts,
	    pa,
	    ta,
	    la,
	    fs
	    );

	    // Remove aleady inserted points
	    clear_inserted<<<bpg, tpb>>>
	    (
	    array_size,
	    pa,
	    ta,
	    la,
	    fs,
	    nm
	    );

	    // Sort everything so that positives are up front
	    thrust::sort(thrust::make_zip_iterator(
	                    thrust::make_tuple(thrust::device_ptr<int>(pa),
	                                       thrust::device_ptr<int>(ta),
	                                       thrust::device_ptr<int>(la),
	                                       thrust::device_ptr<int>(fs),
	                                       thrust::device_ptr<int>(nm))),
	                 thrust::make_zip_iterator(
	                    thrust::make_tuple(thrust::device_ptr<int>(pa + array_capacity),
	                    				   thrust::device_ptr<int>(ta + array_capacity),
	                    				   thrust::device_ptr<int>(la + array_capacity),
	                    				   thrust::device_ptr<int>(fs + array_capacity),
	                    				   thrust::device_ptr<int>(nm + array_capacity))),
	                 tuple_comp<int>());

	    cudaDeviceSynchronize();

	    // We now need the new size of the association arrays...
	    const int old_array_size = array_size;

	    const auto iter = thrust::find(thrust::device_ptr<int>(pa),
	    		                       thrust::device_ptr<int>(pa + array_capacity),
	    		                       -1);

	    array_size = thrust::distance(thrust::device_ptr<int>(pa), iter);

	    // And now the number of new tetrahedra added to the mesh
	    const int num_new_tetra = thrust::device_ptr<int>(fract_offsets)[old_array_size - 1];
	    num_tetra += num_new_tetra;

//#define print
#ifdef print
	    std::cout << "Total number of tetrahedra in the mesh after insertion round : " << num_tetra << std::endl;
	    std::cout << "New size of the association arrays : " << array_size << std::endl;
#undef print
#endif

	    // Massive visual data dump
	    /*array_dump<<<bpg, tpb>>>
	    (
	    array_capacity,
	    pa,
	    ta,
	    la,
	    fs
	    );*/

	    //cudaDeviceSynchronize();

	    device_free(assoc_offsets);
	    device_free(fract_offsets);
	    device_free(bucket_starts);
	    device_free(nominated_tetra);
	    device_free(locked_tetra);
	}

	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);

	std::cout << "Tetrahedralization took " << milliseconds / 1000 << " seconds" << std::endl;
	std::cout << "Rate : " << num_grid_points / (milliseconds / 1000) << " points/second" << std::endl;

#define debug
#ifdef debug
	// Check to make sure all the orientations are correct
	assert_no_flat_tetra<<<bpg, tpb>>>
	(
	device_mesh,
	device_points,
	predConsts,
	num_tetra
	);

	cudaDeviceSynchronize();
#undef debug
#endif

    cuDelete(&predConsts);
    device_free(nm);
	device_free(ta);
	device_free(pa);
	device_free(la);
	device_free(fs);
}
