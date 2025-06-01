#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include "body.h" // Include the Body class

class Simulation
{
private:
    std::vector<Body> bodies; // List of bodies in the simulation
    double timeStep;          // Time step for the simulation

public:
    // Constructor
    Simulation(const std::vector<Body> &initialBodies, double timeStep);

    // Method to update the simulation state
    void update();

    // Method to output the current state (optional)
    void output(int step) const;

    // Additional methods (e.g., for initialization or visualization)
};

#endif // SIMULATION_H