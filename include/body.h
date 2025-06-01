class Body
{
public:
	Body();													 // Random initialization
	Body(double mass);										 // Random position and velocity
	Body(double mass, double x_position, double y_position); // Starting velocity of 0
	Body(double mass, double x_position, double y_position, double x_velocity, double y_velocity);
	const double mass;
	double position[2]; // X coordinate then Y
	double velocity[2]; // X coordinate then Y
};