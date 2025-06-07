#include "simulation.h"
#include "visualizer.h"
#include <gtkmm/application.h>
#include "barnes_hutt.h"
#include <vector>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <thread>
#include <stdio.h>

const double G_CONST = 0.000000000066743;

void initialize_solar_system(Body* bodies) {
    // Values from https://nssdc.gsfc.nasa.gov/planetary/factsheet/
    // n must be equal to 10
    new (&bodies[0]) Body(1.989 * 1e30, 0.0 * 1e9, 0.0, 0.0, 0.0);
    new (&bodies[1]) Body(0.330 * 1e24, 57.90 * 1e9, 0.0, 0.0, 47400.0);
    new (&bodies[2]) Body(4.870 * 1e24, 108.2 * 1e9, 0.0, 0.0, 35000.0);
    new (&bodies[3]) Body(5.970 * 1e24, 149.6 * 1e9, 0.0, 0.0, 29800.0);
    new (&bodies[4]) Body(0.642 * 1e24, 228.0 * 1e9, 0.0, 0.0, 24100.0);
    new (&bodies[5]) Body(189.8 * 1e25, 778.5 * 1e9, 0.0, 0.0, 13100.0);
    new (&bodies[6]) Body(568.0 * 1e24, 1432 * 1e9, 0.0, 0.0, 9700.0);
    new (&bodies[7]) Body(86.80 * 1e24, 2867 * 1e9, 0.0, 0.0, 6800.0);
    new (&bodies[8]) Body(102.0 * 1e24, 4515 * 1e9, 0.0, 0.0, 5400.0);
    new (&bodies[9]) Body(0.073 * 1e24, 150.0 * 1e9, 0.0, 0.0, 30800.0); // Moon!
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
                double dx = bodies[j].position[0] - bodies[i].position[0];
                double dy = bodies[j].position[1] - bodies[i].position[1];
                double dist_sq = dx * dx + dy * dy + 1e-6;
                double dist = std::sqrt(dist_sq);
                double force = G_CONST * bodies[i].mass * bodies[j].mass / dist_sq;
                forces[(i * n + j) * 2] = force * dx / dist;
                forces[(i * n + j) * 2 + 1] = force * dy / dist;
                forces[(j * n + i) * 2] = -forces[(i * n + j) * 2];
                forces[(j * n + i) * 2 + 1] = -forces[(i * n + j) * 2 + 1];

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

void barnes_hut_simulation(Body* bodies, int n, double* all_positions, double step_time, int total_time_steps, double theta, BHTree*& last_tree) {
    for (int t = 0; t < total_time_steps; t++) {
        Quad root_quad(0.5e11, 0.5e11, 2e11); // adjust bounding box as needed
        BHTree tree(root_quad);
        for (int i = 0; i < n; ++i)
            tree.insert(&bodies[i]);

        for (int i = 0; i < n; ++i)
            tree.update_force(&bodies[i], theta, G_CONST);

        for (int i = 0; i < n; ++i) {
            bodies[i].position[0] += bodies[i].velocity[0] * step_time;
            bodies[i].position[1] += bodies[i].velocity[1] * step_time;
        }

        store_positions(bodies, n, all_positions, t);
    }
}

int main(int argc, char** argv)
{
    // Inputs:
    const int threads = 16; // The number of threads
    const double step_time = 864000.0; // Time, in secs, in between each time step
    const int total_time_steps = 1000; // The number of steps to simulate
    const int n = 1000; // Number of bodies
    const std::vector<double> initial_body_positions[n * 2]; // The position of bodies at time step t=0: {x1, y1, x2, y2, ...}
    const std::vector<double> initial_body_velocities[n * 2]; // The velocities of bodies at time step t=0: {vx1, vy1, vx2, vy2, ...}

    // Allocate memory
    // The matrix of forces between body i and j.
    // forces[(i * n + j) * 2] for x, forces[(i * n + j) * 2 + 1] for y
    // Additionally, forces[(i * n + j) * 2] == forces[(j * n + i) * 2] since matrix is symmetrical
    // Stores position of all bodies for every time step
    // all_positions[(i * n + j) * 2 (+1)] to access position_x of body j for timestep i (+1 for position_y)
    double* bis_forces = new double[n * n * 2];                         // Force matrix
    double* bis_all_positions = new double[n * total_time_steps * 2];   // Position buffer
    Body* bis_bodies = new Body[n];                                     // Array of Body instances
    
    // Initialize bodies
    for (int i = 0; i < n; ++i) {
        new (&bis_bodies[i]) Body();  // Placement new to construct in-place
    }
    // Run the sequential simulation
    auto start = std::chrono::steady_clock::now();
    sequential_simulation(bis_bodies, n, bis_forces, bis_all_positions, step_time, total_time_steps);
    auto end = std::chrono::steady_clock::now();
    auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Sequential simulation time: " << elapsed_us << " μs\n";


    // Reallocate memory for concurrent
    double* forces = new double[n * n * 2];                         // Force matrix
    double* all_positions = new double[n * total_time_steps * 2];   // Position buffer
    Body* bodies = new Body[n];
    for (int i = 0; i < n; ++i) {
        new (&bodies[i]) Body();  // Placement new to construct in-place
    }
    // Run concurrent simulation
    start = std::chrono::steady_clock::now();
    thread_creating_simulation(bodies, n, forces, all_positions, step_time, total_time_steps, threads);
    end = std::chrono::steady_clock::now();
    elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Thread Creating concurrent simulation time: " << elapsed_us << " μs\n";

    // Reallocate memory for Barnes-Hutt simulation
    double* bh_positions = new double[n * total_time_steps * 2];
    Body* bh_bodies = new Body[n];
    for (int i = 0; i < n; ++i) {
        new (&bh_bodies[i]) Body(bis_bodies[i]);
    }
    std::cout << "First body position: " << bh_bodies[0].position[0] << ", " << bh_bodies[0].position[1] << std::endl; //debug
    std::cout << "First body velocity: " << bh_bodies[0].velocity[0] << ", " << bh_bodies[0].velocity[1] << std::endl; //debug

    BHTree* final_tree = nullptr;
    // Run Barnes-Hut simulation
    start = std::chrono::steady_clock::now();
    barnes_hut_simulation(bh_bodies, n, bh_positions, step_time, total_time_steps, 0.5, final_tree); // theta = 0.5 is standard
    end = std::chrono::steady_clock::now();
    elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Barnes-Hut simulation time: " << elapsed_us << " μs\n";
    // Compare Barnes-Hut and threaded simulation output 
    double diff_sum = 0.0;
    for (int i = 0; i < n; i++) {
        double dx = bh_positions[(total_time_steps - 1) * n * 2 + i * 2] -
                all_positions[(total_time_steps - 1) * n * 2 + i * 2];
        double dy = bh_positions[(total_time_steps - 1) * n * 2 + i * 2 + 1] -
                all_positions[(total_time_steps - 1) * n * 2 + i * 2 + 1];
        diff_sum += std::sqrt(dx * dx + dy * dy);
    }
    std::cout << "Average final position difference per body: " << diff_sum / n << " meters\n";
    std::cout << "Final x of body 0: " << bh_positions[(total_time_steps - 1) * n * 2 + 0] << std::endl; //test

    // Start GTK application with the visualizer
    auto app = Gtk::Application::create(argc, argv, "org.simulation.nbody");
    Visualizer vis(all_positions, n, total_time_steps);
    //Visualizer vis(bh_positions, n, total_time_steps, final_tree); //visualize the barnes simulation
    app->run(vis);

    // Free resources
    delete[] bis_forces;
    delete[] bis_all_positions;
    delete[] bis_bodies;
    // Free memory
    delete[] bh_positions;
    delete[] bh_bodies;
    delete final_tree;
    delete[] forces;
    delete[] all_positions;
    delete[] bodies;
    return 0;
}