#include "simulation.h"
#include <cmath>
#include <vector>
#include <mutex>
#include <thread>
#include <chrono>
const double G_CONST = 6.67430e-11;

void store_positions(Body* bodies, int n, double* all_positions, int current_time_step)
{
    for (int i = 0; i < n; i++)
    {
        all_positions[(current_time_step * n + i) * 2] = bodies[i].position[0];
        all_positions[(current_time_step * n + i) * 2 + 1] = bodies[i].position[1];
    }
}

void mutex_using_simulations(Body* bodies, double* forces, double* distances,
                             int n, double step_time, int step_sync, int threads,
                             std::mutex* mutex)
{

    while (true) {
        {
            std::lock_guard<std::mutex> lock(*mutex);

            for (int i = 0; i < n; ++i) {
                double fx = 0.0, fy = 0.0;

                for (int j = 0; j < n; ++j) {
                    if (i == j) continue;

                    double dx = bodies[j].position[0] - bodies[i].position[0];
                    double dy = bodies[j].position[1] - bodies[i].position[1];
                    double dist_sq = dx*dx + dy*dy + 1e-6;
                    double dist = std::sqrt(dist_sq);
                    double F = G_CONST * bodies[i].mass * bodies[j].mass / dist_sq;

                    fx += F * dx / dist;
                    fy += F * dy / dist;
                }

                bodies[i].velocity[0] += fx / bodies[i].mass * step_time;
                bodies[i].velocity[1] += fy / bodies[i].mass * step_time;
            }

            for (int i = 0; i < n; ++i) {
                bodies[i].position[0] += bodies[i].velocity[0] * step_time;
                bodies[i].position[1] += bodies[i].velocity[1] * step_time;
            }
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(33));
    }
}


// The following functions create and join threads EVERY time step
void force_calculating_thread(Body* bodies, int n, int begin, int end, double* forces) {
    for (int i = begin; i < end; i++) {
        for (int j = i; j < n; j++) {
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

void update_calculating_thread(Body* bodies, int n, int begin, int end, double* forces, double* all_positions, double step_time) {
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

void thread_creating_simulation(Body* bodies, int n, double* forces, double* all_positions, double step_time, int total_time_steps, int threads) {
    for (int current_time_step = 0; current_time_step < total_time_steps; current_time_step++) {
        std::vector<std::thread> workers(threads - 1);
        // Update forces
        int total_operations = (n * n + n) / 2;
        int operations_to_perform = 0;
        int current_begin = 0;
        int current_end = 0;
        int thread_number = 0;
        while (thread_number < threads - 1) {
            current_end++;
            operations_to_perform += n - current_end + 1;
            if (operations_to_perform >= total_operations / (threads - 1)) {
                operations_to_perform -= total_operations / (threads - 1);
                workers[thread_number] = std::thread(&force_calculating_thread, bodies, n, current_begin, current_end, forces);
                current_begin = current_end;
                thread_number++;
            }
        }
        for (int i = 0; i < threads - 1; i++) {
            workers[i].join();
        }
        // Update bodies
        current_begin = 0;
        current_end = 0;
        int block_size = n / (threads - 1) + 1;
        for (int i = 0; i < threads - 1; i++) {
            current_end += block_size;
            if (current_end > n) {
                current_end = n;
            }
            workers[i] = std::thread(&update_calculating_thread, bodies, n, current_begin, current_end, forces, all_positions, step_time);
            current_begin = current_end;
        }
        for (int i = 0; i < threads - 1; i++) {
            workers[i].join();
        }
        store_positions(bodies, n, all_positions, current_time_step);
    }
}