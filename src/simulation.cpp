#include "simulation.h"
#include <iostream>
#include <cmath> // For mathematical calculations

// Constructor: Initialize the simulation with bodies and a time step
Simulation::Simulation(const std::vector<Body> &initialBodies, double timeStep)
    : bodies(initialBodies), timeStep(timeStep) {}

// Method to update the simulation state
void Simulation::update()
{
    // Calculate forces between bodies
    for (size_t i = 0; i < bodies.size(); ++i)
    {
        for (size_t j = i + 1; j < bodies.size(); ++j)
        {
            // Calculate the gravitational force between bodies[i] and bodies[j]
            double dx = bodies[j].position[0] - bodies[i].position[0];
            double dy = bodies[j].position[1] - bodies[i].position[1];
            double distance = std::sqrt(dx * dx + dy * dy);

            // Avoid division by zero
            if (distance == 0)
                continue;

            // Gravitational constant (example value)
            const double G = 6.67430e-11;

            // Force magnitude
            double force = G * bodies[i].mass * bodies[j].mass / (distance * distance);

            // Force components
            double fx = force * dx / distance;
            double fy = force * dy / distance;

            // Apply forces to the bodies
            bodies[i].velocity[0] += fx / bodies[i].mass * timeStep;
            bodies[i].velocity[1] += fy / bodies[i].mass * timeStep;

            bodies[j].velocity[0] -= fx / bodies[j].mass * timeStep;
            bodies[j].velocity[1] -= fy / bodies[j].mass * timeStep;
        }
    }

    // Update positions of all bodies
    for (Body &body : bodies)
    {
        body.position[0] += body.velocity[0] * timeStep;
        body.position[1] += body.velocity[1] * timeStep;
    }
}

// Method to output the current state of the simulation
void Simulation::output(int step) const
{
    std::cout << "Step " << step << ":\n";
    for (const Body &body : bodies)
    {
        std::cout << "Mass: " << body.mass
                  << ", Position: (" << body.position[0] << ", " << body.position[1] << ")"
                  << ", Velocity: (" << body.velocity[0] << ", " << body.velocity[1] << ")\n";
    }
}