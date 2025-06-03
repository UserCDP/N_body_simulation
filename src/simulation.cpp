#include "simulation.h"
#include <cmath>
#include <mutex>
#include <thread>
#include <chrono>

void run_simulation(Body* bodies, double* forces, double* distances,
                    int n, double step_time, int step_sync, int threads,
                    std::mutex* mutex)
{
    const double G = 6.67430e-11;

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
                    double F = G * bodies[i].mass * bodies[j].mass / dist_sq;

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
