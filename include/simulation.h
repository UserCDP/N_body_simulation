#pragma once
#include "body.h"
#include <mutex>

void store_positions(Body* bodies, int n, double* all_positions, int current_time_step);

void mutex_using_simulation(Body* bodies, double* forces, double* distances,
                    int n, double step_time, int step_sync, int threads,
                    std::mutex* mutex); // TODO Wrong inputs

void thread_creating_simulation(Body* bodies, int n, double* forces, double* all_positions, double step_time, int total_time_steps, int threads);