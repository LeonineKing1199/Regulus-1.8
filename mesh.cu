/*
 * mesh.cu
 *
 *  Created on: Apr 9, 2015
 *      Author: christian
 */




#include "structures.h"

mesh::mesh(const int gl)
{
	grid_length = gl;

	// May be extraneous but we check to make sure that
	// the grid length does not overflow signed 4 bytes ints
	{
		long int proposed_num_points = (long int ) grid_length * grid_length * grid_length;
		assert(proposed_num_points <= INT_MAX);

		num_grid_points = proposed_num_points;
	}

	// Number of points required for the root tetrahedron
	num_root_points = 4;

	// Get the number of bytes for allocations
	//const size_t point_bytes = (num_grid_points + num_root_points) * sizeof(point);
	//const size_t tetra_bytes = 8 * num_grid_points * sizeof(tetra);

	// Allocate memory
	//cudaMallocHost(&user_points, point_bytes);
	//cudaMallocHost(&user_mesh, tetra_bytes);

	//cudaMalloc(&device_points, point_bytes);
	//cudaMalloc(&device_mesh, tetra_bytes);

	// Modified code to use convenience wrappers
	const  int est_num_tetra = 8 * num_grid_points;

	user_points = host_malloc<point>(num_grid_points + num_root_points);
	user_mesh = host_malloc<tetra>(est_num_tetra);

	device_points = device_malloc<point>(num_grid_points + num_root_points);
	device_mesh = device_malloc<tetra>(est_num_tetra);

	// Initialize data
	build_domain();
}

mesh::~mesh(void)
{
	//cudaFreeHost(user_points);
	//cudaFreeHost(user_mesh);

	//cudaFree(device_points);
	//cudaFree(device_mesh);

	host_free(user_points);
	host_free(user_mesh);

	device_free(device_points);
	device_free(device_mesh);
}

void mesh::build_domain(void)
{
	// Assign Cartesian distribution
	for (int i = 0; i < num_grid_points; ++i)
	{
        const float x = (float ) (i / (grid_length * grid_length)),
                    y = (float ) ((i / grid_length) % grid_length),
                    z = (float ) (i % grid_length);

        point tmp = { x, y, z };
        user_points[i] = tmp;
	}

	// Minimum length required to encompass the entire grid
	root_edge_length = (grid_length - 1) * 3;

	const point rt_pt0 = { 0, 0, 0 };
	const point rt_pt1 = { (float ) root_edge_length, 0, 0 };
	const point rt_pt2 = { 0, (float ) root_edge_length, 0 };
	const point rt_pt3 = { 0, 0, (float ) root_edge_length };

	// Write root points to the end of the buffer
	user_points[num_grid_points + 0] = rt_pt0;
	user_points[num_grid_points + 1] = rt_pt1;
	user_points[num_grid_points + 2] = rt_pt2;
	user_points[num_grid_points + 3] = rt_pt3;

	// Finally, build root tetrahedron
	const tetra tmp = { num_grid_points + 0,
						num_grid_points + 1,
						num_grid_points + 2,
						num_grid_points + 3 };

	user_mesh[0] = tmp;
	num_tetra = 1;

	// Copy everything back to the device
	cudaMemcpy(device_points,
			   user_points,
			   (num_grid_points + num_root_points) * sizeof(point),
			   cudaMemcpyHostToDevice);

	cudaMemcpy(device_mesh,
			   user_mesh,
			   num_tetra * sizeof(tetra),
			   cudaMemcpyHostToDevice);
}

void mesh::print_user_points(void) const
{
	std::cout << "Printing point set...\n";

	for (int i = 0; i < num_grid_points + num_root_points; ++i)
	{
		const point tmp = user_points[i];
		std::cout << tmp.c[0] << ", " << tmp.c[1] << ", " << tmp.c[2] << "\n";
	}

	std::cout << std::endl;
}
