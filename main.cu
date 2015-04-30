/*
 * main.cu
 *
 *  Created on: Apr 9, 2015
 *      Author: christian
 */




#include "structures.h"

int main(void)
{
	// Grid length
	const int gl = 126;

	mesh cart_mesh(gl);
	cart_mesh.triangulate();

	//cart_mesh.print_user_points();

	return 0;
}
