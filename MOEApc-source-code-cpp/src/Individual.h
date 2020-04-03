#ifndef Individual_H
#define Individual_H

class CIndividual
{
public:
	double * x; // decision vector

	double * f; // objective vector

	double rho; // distance to reference point

	double * theta; // polar coordinates (angles)	

public:
	// constructor
	CIndividual(void) 
	{
	}
	// destructor
	~CIndividual(void) 
	{
	}
};


#endif