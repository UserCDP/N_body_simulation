#include "simulation.h"
#include "visualizer.h"
#include <gtkmm/application.h>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <thread>
#include <stdio.h>

const double G_CONST = 0.000000000066743;

void store_positions(Body *bodies, int n, double *all_positions, int current_time_step)
{
    for (int i = 0; i < n; i++)
    {
        all_positions[(current_time_step * n + i) * 2] = bodies[i].position[0];
        all_positions[(current_time_step * n + i) * 2 + 1] = bodies[i].position[1];
    }
}

/**
 * @brief Runs a sequential N-body simulation for a given number of bodies and time steps.
 *
 * This function simulates the motion of 'n' bodies under mutual gravitational attraction
 * over a specified number of time steps. At each step, it computes the pairwise gravitational
 * forces between all bodies, updates their velocities and positions accordingly, and stores
 * the positions for each time step.
 *
 * @param bodies            Pointer to an array of Body objects representing the bodies in the simulation.
 * @param n                 The number of bodies in the simulation.
 * @param forces            Pointer to a pre-allocated array to store computed forces between bodies.
 *                          The array should have size [n * n * 2] (for x and y components).
 * @param all_positions     Pointer to a pre-allocated array to store the positions of all bodies at each time step.
 *                          The array should have size [n * total_time_steps * 2].
 * @param step_time         The time interval (in seconds) between each simulation step.
 * @param total_time_steps  The total number of simulation steps to perform.
 */
void sequential_simulation(Body *bodies, int n, double *forces, double *all_positions, double step_time, int total_time_steps)
{
    for (int current_time_step = 0; current_time_step < total_time_steps; current_time_step++)
    {
        // Update forces
        for (int i = 0; i < n; i++)
        {
            for (int j = i; j < n; j++)
            {
                double distance_x = (bodies[i].position[0] - bodies[j].position[0]) * (bodies[i].position[0] - bodies[j].position[0]);
                double distance_y = (bodies[i].position[1] - bodies[j].position[1]) * (bodies[i].position[1] - bodies[j].position[1]);
                double force = G_CONST * bodies[i].mass * bodies[j].mass / (distance_x + distance_y);
                forces[(i * n + j) * 2] = force * sqrt(distance_x / (distance_x + distance_y));
                forces[(i * n + j) * 2 + 1] = force * sqrt(distance_y / (distance_x + distance_y));
                forces[(j * n + i) * 2] = forces[(i * n + j) * 2];
                forces[(j * n + i) * 2 + 1] = forces[(i * n + j) * 2 + 1];
            }
        }
        // Update velocities and positions
        for (int i = 0; i < n; i++)
        {
            double force_x = 0.0;
            double force_y = 0.0;
            for (int j = 0; j < n; j++)
            {
                if (j == i)
                {
                    continue;
                }
                force_x += forces[(i * n + j) * 2];
                force_y += forces[(i * n + j) * 2 + 1];
            }
            bodies[i].velocity[0] += force_x / bodies[i].mass * step_time;
            bodies[i].velocity[1] += force_y / bodies[i].mass * step_time;
            bodies[i].position[0] += bodies[i].velocity[0] * step_time;
            bodies[i].position[1] += bodies[i].velocity[1] * step_time;
        }
        // Store new positions
        store_positions(bodies, n, all_positions, current_time_step);
    }
}

int main()
{
    // Inputs:
    const int threads = 16;                                   // The number of threads
    const double step_time = 3600.0;                          // Time, in secs, in between each time step
    const int total_time_steps = 1024;                        // The number of steps to simulate
    const int n = 8;                                          // Number of bodies
    const std::vector<double> initial_body_positions[n * 2];  // The position of bodies at time step t=0: {x1, y1, x2, y2, ...}
    const std::vector<double> initial_body_velocities[n * 2]; // The velocities of bodies at time step t=0: {vx1, vy1, vx2, vy2, ...}

    int t = 0; // The current time step, multiply with step_time for current time in secs.

    // The matrix of forces between body i and j.
    // forces[(i * n + j) * 2] for x, forces[(i * n + j) * 2 + 1] for y
    // Additionally, forces[(i * n + j) * 2] == forces[(j * n + i) * 2] since matrix is symmetrical
    double *forces = (double *)malloc(n * n * 2 * sizeof(double));

    // Stores position of all bodies for every time step
    // all_positions[(i * n + j) * 2 (+1)] to access position_x of body j for timestep i (+1 for position_y)
    double *all_positions = (double *)malloc(n * total_time_steps * 2 * sizeof(double));

    Body *bodies = (Body *)malloc(n * sizeof(Body));

    // Initialize Body instances
    for (int i = 0; i < n; i++)
    {
        new (&bodies[i]) Body();
    }

    // std::thread simulation_thread(run_simulation, bodies, n, forces, all_positions, step_time, total_time_steps, threads);
    // simulation_thread.join();

    auto start = std::chrono::steady_clock::now();
    sequential_simulation(bodies, n, forces, all_positions, step_time, total_time_steps);
    auto finish = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

    printf("Time elapsed: %d\n", (int)elapsed);

    std::thread simulation_thread([&]() {
        run_simulation(bodies, forces, distances, n, step_time, step_sync, threads, &bodies_mutex);
    });

    auto app = Gtk::Application::create(argc, argv, "org.simulation.nbody");
    Visualizer vis(&bodies, &bodies_mutex);
    app->run(vis);

    // TODO Concurrent Simulation (call function in main.cpp, DO NOT clog main()!)

    free(forces);
    free(all_positions);
    free(bodies);
    return 0;
}