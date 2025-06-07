#include "simulation.h"
#include <cmath>
#include <vector>
#include <mutex>
#include <thread>
#include <chrono>
const double G_CONST = 6.67430e-11;

/**
 * @brief Stores the positions of bodies at the current time step in the all_positions array.
 *
 * This function takes an array of Body structures and stores their positions in a flat array
 * for visualization or further processing.
 *
 * @param bodies Pointer to the array of Body structures representing the bodies in the simulation.
 * @param n The number of bodies in the simulation.
 * @param all_positions Pointer to the array where positions will be stored.
 * @param current_time_step The current time step index, used to determine where to store the positions.
 * The positions are stored in a flat array format, where each body's position is stored as two consecutive elements
 * (x and y coordinates) for each time step.
 */
void store_positions(Body *bodies, int n, double *all_positions, int current_time_step)
{
    for (int i = 0; i < n; i++)
    {
        all_positions[(current_time_step * n + i) * 2] = bodies[i].position[0];
        all_positions[(current_time_step * n + i) * 2 + 1] = bodies[i].position[1];
    }
}

/**
 * @brief Runs an N-body simulation using mutexes to synchronize access to shared resources.
 *
 * This function simulates the motion of bodies under mutual gravitational attraction,
 * updating their positions and velocities in a loop.
 *
 * @param bodies Pointer to the array of Body structures representing the bodies in the simulation.
 * @param forces Pointer to the array where computed force vectors between bodies will be stored.
 * @param distances Pointer to the array where distances between bodies will be stored (not used in this implementation).
 * @param n The number of bodies in the simulation.
 * @param step_time The time increment for each simulation step.
 * @param step_sync The synchronization step (not used in this implementation).
 * @param threads The number of threads to use for parallel computation (not used in this implementation).
 * @param mutex Pointer to a mutex used to synchronize access to shared resources.
 * @note This function runs indefinitely in a loop, simulating the motion of bodies until externally stopped.
 * @note The mutex is used to ensure that only one thread can access and modify the bodies' positions and velocities at a time.
 */
void mutex_using_simulations(Body *bodies, double *forces, double *distances,
                             int n, double step_time, int step_sync, int threads,
                             std::mutex *mutex)
{

    while (true)
    {
        {
            std::lock_guard<std::mutex> lock(*mutex);

            for (int i = 0; i < n; ++i)
            {
                double fx = 0.0, fy = 0.0;

                for (int j = 0; j < n; ++j)
                {
                    if (i == j)
                        continue;

                    double dx = bodies[j].position[0] - bodies[i].position[0];
                    double dy = bodies[j].position[1] - bodies[i].position[1];
                    double dist_sq = dx * dx + dy * dy + 1e-6;
                    double dist = std::sqrt(dist_sq);
                    double F = G_CONST * bodies[i].mass * bodies[j].mass / dist_sq;

                    fx += F * dx / dist;
                    fy += F * dy / dist;
                }

                bodies[i].velocity[0] += fx / bodies[i].mass * step_time;
                bodies[i].velocity[1] += fy / bodies[i].mass * step_time;
            }

            for (int i = 0; i < n; ++i)
            {
                bodies[i].position[0] += bodies[i].velocity[0] * step_time;
                bodies[i].position[1] += bodies[i].velocity[1] * step_time;
            }
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(33));
    }
}

/**
 * @brief Calculates the gravitational forces between bodies in a multithreaded manner.
 *
 * This function computes the gravitational forces acting on each body due to every other body
 * in the simulation. It uses a pairwise approach, where the force between each pair of bodies is calculated
 * and stored in a shared forces array. The forces are stored in a flat array format, where each pair of bodies
 * (i, j) has its force components stored as:
 * forces[(i * n + j) * 2] for the x-component and forces[(i * n + j) * 2 + 1] for the y-component.
 * The forces are symmetric, meaning that forces[(i * n + j) * 2] == -forces[(j * n + i) * 2].
 *
 * @param bodies Pointer to the array of Body structures representing the bodies in the simulation.
 * @param n The number of bodies in the simulation.
 * @param begin The starting index for the thread's computation.
 * @param end The ending index for the thread's computation.
 * @param forces Pointer to the array where computed forces between bodies will be stored.
 * The array should have size [n * n * 2] (for x and y components).
 */
void force_calculating_thread(Body *bodies, int n, int begin, int end, double *forces)
{
    for (int i = begin; i < end; i++)
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
}

/**
 * @brief Updates the positions and velocities of bodies based on the calculated forces.
 *
 * This function iterates over a range of bodies and updates their velocities and positions
 * based on the forces acting on them. The forces are assumed to be pre-calculated and stored
 * in the forces array, where each body's force is calculated as the sum of forces from all other bodies.
 * The positions are updated using the formula:
 * position += velocity * step_time, where velocity is updated based on the force divided by the body's mass.
 *
 * @param bodies Pointer to the array of Body structures representing the bodies in the simulation.
 * @param n The number of bodies in the simulation.
 * @param begin The starting index for the thread's computation.
 * @param end The ending index for the thread's computation.
 * @param forces Pointer to the array where computed forces between bodies are stored.
 * @param all_positions Pointer to the array where positions of all bodies at each time step are stored.
 * @param step_time The time increment for each simulation step.
 */
void update_calculating_thread(Body *bodies, int n, int begin, int end, double *forces, double *all_positions, double step_time)
{
    for (int i = begin; i < end; i++)
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
}

/**
 * @brief Runs a multithreaded N-body simulation by creating threads for force calculation and position updates.
 *
 * This function simulates the motion of 'n' bodies under mutual gravitational attraction
 * over a specified number of time steps. It divides the workload among multiple threads,
 * calculating the pairwise gravitational forces between bodies and updating their velocities and positions accordingly.
 * It uses a vector of threads to manage the concurrent execution of force calculations and position updates.
 *
 * @param bodies Pointer to an array of Body objects representing the bodies in the simulation.
 * @param n The number of bodies in the simulation.
 * @param forces Pointer to a pre-allocated array to store computed forces between bodies.
 *                          The array should have size [n * n * 2] (for x and y components).
 * @param all_positions Pointer to a pre-allocated array to store the positions of all bodies at each time step.
 *                          The array should have size [n * total_time_steps * 2].
 * @param step_time The time interval (in seconds) between each simulation step.
 * @param total_time_steps The total number of simulation steps to perform.
 * @param threads The number of threads to use for parallel computation.
 * @note This function creates a number of threads equal to `threads - 1`, as one thread is used for the main simulation loop.
 * Each thread is responsible for calculating forces or updating positions for a subset of bodies.
 */
void thread_creating_simulation(Body *bodies, int n, double *forces, double *all_positions, double step_time, int total_time_steps, int threads)
{
    for (int current_time_step = 0; current_time_step < total_time_steps; current_time_step++)
    {
        std::vector<std::thread> workers(threads - 1);
        // Update forces
        int total_operations = (n * n + n) / 2;
        int operations_to_perform = 0;
        int current_begin = 0;
        int current_end = 0;
        int thread_number = 0;
        while (thread_number < threads - 1)
        {
            current_end++;
            operations_to_perform += n - current_end + 1;
            if (operations_to_perform >= total_operations / (threads - 1))
            {
                operations_to_perform -= total_operations / (threads - 1);
                workers[thread_number] = std::thread(&force_calculating_thread, bodies, n, current_begin, current_end, forces);
                current_begin = current_end;
                thread_number++;
            }
        }
        for (int i = 0; i < threads - 1; i++)
        {
            workers[i].join();
        }
        // Update bodies
        current_begin = 0;
        current_end = 0;
        int block_size = n / (threads - 1) + 1;
        for (int i = 0; i < threads - 1; i++)
        {
            current_end += block_size;
            if (current_end > n)
            {
                current_end = n;
            }
            workers[i] = std::thread(&update_calculating_thread, bodies, n, current_begin, current_end, forces, all_positions, step_time);
            current_begin = current_end;
        }
        for (int i = 0; i < threads - 1; i++)
        {
            workers[i].join();
        }
        store_positions(bodies, n, all_positions, current_time_step);
    }
}