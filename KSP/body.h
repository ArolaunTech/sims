#ifndef BODY_H
#define BODY_H

struct Body {
	double grav_parameter;
	double radius;
	double rotational_period;

	Body(double newGravParameter = 0, double newRadius = 0, double newRotationalPeriod = 0) : 
		grav_parameter(newGravParameter), 
		radius(newRadius), 
		rotational_period(newRotationalPeriod) {}
};

#endif