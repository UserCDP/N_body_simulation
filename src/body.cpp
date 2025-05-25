#include "body.h"
#include "simulation.h"

int main() {
	// Inputs:
	const int threads = 16;  // The number of threads
	const double step_time = 3600.0;  // Time, in secs, in between each time step
	const int step_sync = 24;  // The number of steps before calculation results must be updated
	const int n = 8;  // Number of bodies
	const double[n * 2] initial_body_positions;  // The position of bodies at time step t=0: {x1, y1, x2, y2, ...}
	const double[n * 2] initial_body_velocities;  // The velocities of bodies at time step t=0: {vx1, vy1, vx2, vy2, ...}

	int t = 0;  // The current time step, multiply with step_time for current time in secs.

	// The matrix of forces between body i and j.
	// forces[(i * n + j) * 2] for x, forces[(i * n + j) * 2 + 1] for y
	// Additionally, forces[(i * n + j) * 2] == forces[(j * n + i) * 2] since matrix is symmetrical
	double* forces = malloc(n * n * 2 * sizeof(double));

	// Same as with matrix of forces without multiplying by 2 
	// (no x or y distinction, just the norm value between body i and j). (Still symmetrical)
	double* distances = malloc(n * n * sizeof(double));

	Body* bodies = malloc(n * sizeof(Body));


	// TODO Initialize Body instances as well as the forces and distances matrices by calling functions in another file (initialization.h and .cpp?)

	// TODO Sequential Simulation (call function in body.cpp, do not clog main()!)

	// TODO Concurrent Simulation (call function in body.cpp, DO NOT clog main()! Really don't for concurrent pls)


	free(forces);
	free(distances);
	return 0;
}