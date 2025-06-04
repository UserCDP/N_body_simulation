#include "body.h"
#include <cstdlib> // For random number generation
#include <ctime>   // For seeding random number generator

// Default constructor: Random initialization
Body::Body()
	: mass(static_cast<double>(rand()) / RAND_MAX * 1e28 + 1.0),
	  position{static_cast<double>(rand()) / RAND_MAX * 1e10, static_cast<double>(rand()) / RAND_MAX * 1e10},
	  velocity{((static_cast<double>(rand()) / RAND_MAX) - 0.5) * 1e5, ((static_cast<double>(rand()) / RAND_MAX) - 0.5) * 1e5 }
{
}

// Constructor: Random position and velocity, given mass
Body::Body(double mass)
	: mass(mass),
	  position{ static_cast<double>(rand()) / RAND_MAX * 1e10, static_cast<double>(rand()) / RAND_MAX * 1e10 },
	  velocity{ ((static_cast<double>(rand()) / RAND_MAX) - 0.5) * 1e5, ((static_cast<double>(rand()) / RAND_MAX) - 0.5) * 1e5 }
{
}

// Constructor: Specified mass and position, velocity initialized to 0
Body::Body(double mass, double x_position, double y_position)
	: mass(mass),
	  position{x_position, y_position},
	  velocity{0.0, 0.0} // Velocity set to 0
{
}

// Constructor: Fully specified mass, position, and velocity
Body::Body(double mass, double x_position, double y_position, double x_velocity, double y_velocity)
	: mass(mass),
	  position{x_position, y_position},
	  velocity{x_velocity, y_velocity}
{
}