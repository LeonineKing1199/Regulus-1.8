/*
 * bucket_hash.cu
 *
 *  Created on: Apr 14, 2015
 *      Author: christian
 */




#include "structures.h"

/*
 * Awesome routine from a PhD thesis from UC Davis.
 * Title is "Efficient Hash Tables on the GPU" by
 * Dan Anthony Feliciano Alcantara.
 */

__global__
void find_boundaries
(
const int  num_keys,
const int  num_buckets,
const int * __restrict__ which_bucket,
int * __restrict__ bucket_starts
)
{
    const int thread_id = currThreadID();

    for (int tid = thread_id; tid < num_keys; tid += grid_size())
    {
        // get start and end of each bucket
        const int begin = (tid > 0 ? which_bucket[tid - 1] : 0);
        const int end   = which_bucket[tid];

        // if bucket has length > 0...
        if (begin != end)
        {
            for (int i = begin; i < end; ++i)
            {
                // sets bucket starts value to index of bucket
                bucket_starts[i] = tid;
            }
        }

        // last thread writes number of elements to the rest
        // of the array
        if (tid == num_keys - 1)
        {
            for (int i = end; i < num_buckets; ++i)
            {
                bucket_starts[i] = num_keys;
            }
        }
    }
}
