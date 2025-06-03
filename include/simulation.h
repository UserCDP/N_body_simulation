#pragma once
#include "body.h"
#include <mutex>
void run_simulation(Body* bodies, double* forces, double* distances,
                    int n, double step_time, int step_sync, int threads,
                    std::mutex* mutex);