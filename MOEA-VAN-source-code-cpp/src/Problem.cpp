#include<algorithm>
#include <iostream>
#include "Problem.h"


CProblem::CProblem(void)
{
	initialization();
}


CProblem::~CProblem(void)
{
	delete [] lb;
	delete [] ub;
}




void CProblem::initialization(void)
{
	// READ PROBLEM DEFINITION FROM FILE

	std::ifstream fin; // create ifstream object called fin	
	fin.open("param_problem.dat"); // opens the file

	std::string temp; // temp string
	fin >> temp >> name; // read problem name
	fin >> temp >> m; // read number of objectives
	fin >> temp >> n; // read number of variables
	
	fin.close(); // close the file

	// INITIALIZE LB, UB AND FEVAL
	lb = new double[n];
	ub = new double[n];

	// DTLZ test suite
	if(name.compare("DTLZ1") == 0)
	{
		std::fill(lb, lb+n, 0.0);
		std::fill(ub, ub+n, 1.0);
		feval = DTLZ::DTLZ1;
	}
	else if(name.compare("DTLZ2") == 0)
	{
		std::fill(lb, lb+n, 0.0);
		std::fill(ub, ub+n, 1.0);
		feval = DTLZ::DTLZ2;
	}
	else if(name.compare("DTLZ3") == 0)
	{
		std::fill(lb, lb+n, 0.0);
		std::fill(ub, ub+n, 1.0);
		feval = DTLZ::DTLZ3;
	}
	else if(name.compare("DTLZ4") == 0)
	{
		std::fill(lb, lb+n, 0.0);
		std::fill(ub, ub+n, 1.0);
		feval = DTLZ::DTLZ4;
	}
	else if(name.compare("DTLZ5") == 0)
	{
		std::fill(lb, lb+n, 0.0);
		std::fill(ub, ub+n, 1.0);
		feval = DTLZ::DTLZ5;
	}
	else if(name.compare("DTLZ6") == 0)
	{
		std::fill(lb, lb+n, 0.0);
		std::fill(ub, ub+n, 1.0);
		feval = DTLZ::DTLZ6;
	}
	else if(name.compare("DTLZ7") == 0)
	{
		std::fill(lb, lb+n, 0.0);
		std::fill(ub, ub+n, 1.0);
		feval = DTLZ::DTLZ7;
	}
	// LZ09 test suite
	else if(name.compare("LZ09_F1") == 0)
	{
		std::fill(lb, lb+n, 0.0);
		std::fill(ub, ub+n, 1.0);
		feval = LZ09::F1;
	}
	else if(name.compare("LZ09_F2") == 0)
	{
		std::fill(lb, lb+n, -1.0); lb[0]=0.0;
		std::fill(ub, ub+n, 1.0);
		feval = LZ09::F2;
	}
	else if(name.compare("LZ09_F3") == 0)
	{
		std::fill(lb, lb+n, -1.0); lb[0]=0.0;
		std::fill(ub, ub+n, 1.0);
		feval = LZ09::F3;
	}
	else if(name.compare("LZ09_F4") == 0)
	{
		std::fill(lb, lb+n, -1.0); lb[0]=0.0;
		std::fill(ub, ub+n, 1.0);
		feval = LZ09::F4;
	}
	else if(name.compare("LZ09_F5") == 0)
	{
		std::fill(lb, lb+n, -1.0); lb[0]=0.0;
		std::fill(ub, ub+n, 1.0);
		feval = LZ09::F5;
	}
	else if(name.compare("LZ09_F6") == 0)
	{
		std::fill(lb, lb+n, -2.0); lb[0]=0.0; lb[1]=0.0;
		std::fill(ub, ub+n, 2.0); ub[0]=1.0; ub[1]=1.0;
		feval = LZ09::F6;
	}
	else if(name.compare("LZ09_F7") == 0)
	{
		std::fill(lb, lb+n, 0.0);
		std::fill(ub, ub+n, 1.0);
		feval = LZ09::F7;
	}
	else if(name.compare("LZ09_F8") == 0)
	{
		std::fill(lb, lb+n, 0.0);
		std::fill(ub, ub+n, 1.0);
		feval = LZ09::F8;
	}
	else if(name.compare("LZ09_F9") == 0)
	{
		std::fill(lb, lb+n, -1.0); lb[0]=0.0;
		std::fill(ub, ub+n, 1.0);
		feval = LZ09::F9;
	}
	else if(name.compare("LGZ1") == 0)
	{
		std::fill(lb, lb+n, 0.0);
		std::fill(ub, ub+n, 1.0);
		feval = LGZ::F1;
	}
	else if(name.compare("LGZ2") == 0)
	{
		std::fill(lb, lb+n, 0.0);
		std::fill(ub, ub+n, 1.0);
		feval = LGZ::F2;
	}
	else if(name.compare("LGZ3") == 0)
	{
		std::fill(lb, lb+n, 0.0);
		std::fill(ub, ub+n, 1.0);
		feval = LGZ::F3;
	}
	else if(name.compare("LGZ4") == 0)
	{
		std::fill(lb, lb+n, 0.0);
		std::fill(ub, ub+n, 1.0);
		feval = LGZ::F4;
	}
	else if(name.compare("LGZ5") == 0)
	{
		std::fill(lb, lb+n, 0.0);
		std::fill(ub, ub+n, 1.0);
		feval = LGZ::F5;
	}
	else if(name.compare("LGZ6") == 0)
	{
		std::fill(lb, lb+n, 0.0);
		std::fill(ub, ub+n, 1.0);
		feval = LGZ::F6;
	}
	else if(name.compare("LGZ7") == 0)
	{
		std::fill(lb, lb+n, 0.0);
		std::fill(ub, ub+n, 1.0);
		feval = LGZ::F7;
	}

}


void CProblem::evaluate(CIndividual individual)
{

	feval(individual.x, n, individual.f, m);
}
