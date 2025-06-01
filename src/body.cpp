#include "body.h"
#include <cstdlib> // For random number generation
#include <ctime>   // For seeding random number generator

// Default constructor: Random initialization
Body::Body()
	: mass(static_cast<double>(rand()) / RAND_MAX * 10.0 + 1.0),														 // Random mass between 1 and 10
	  position{static_cast<double>(rand()) / RAND_MAX * 100.0, static_cast<double>(rand()) / RAND_MAX * 100.0},			 // Random position in a 100x100 area
	  velocity{static_cast<double>(rand()) / RAND_MAX * 10.0 - 5.0, static_cast<double>(rand()) / RAND_MAX * 10.0 - 5.0} // Random velocity between -5 and 5
{
}

// Constructor: Random position and velocity, given mass
Body::Body(double mass)
	: mass(mass),
	  position{static_cast<double>(rand()) / RAND_MAX * 100.0, static_cast<double>(rand()) / RAND_MAX * 100.0},			 // Random position in a 100x100 area
	  velocity{static_cast<double>(rand()) / RAND_MAX * 10.0 - 5.0, static_cast<double>(rand()) / RAND_MAX * 10.0 - 5.0} // Random velocity between -5 and 5
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