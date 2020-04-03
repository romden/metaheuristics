/*
Hai-Lin Liu, Fangqing Gu and Qingfu Zhang.
Decomposition of a Multiobjective Optimization Problem into a Number of Simple Multiobjective Subproblems
2013
*/

#include "LGZ.h"

namespace LGZ
{

	void F1(const double * x,  int n, double * f, int m)
	{
		double pi = 3.1415926535897932384626433832795;
		double tj;
		double g = 0.0;

		for(int j = 1; j < n; j++)
		{
			tj = x[j] - sin(0.5*pi*x[0]);
			g+= 2.0*sin(pi*x[0])*(-0.9*tj*tj + pow(abs(tj), 0.6));
		}

		f[0] = (1.0 + g)*x[0];
		f[1] = (1.0 + g)*(1- sqrt(x[0]));
	}


	void F2(const double * x,  int n, double * f, int m)
	{
		double pi = 3.1415926535897932384626433832795;
		double tj;
		double g = 0.0;

		for(int j = 1; j < n; j++)
		{
			tj = x[j] - sin(0.5*pi*x[0]);
			g+= 10.0*sin(pi*x[0])*(abs(tj) / (1+exp(5.0*abs(tj) ) ) );
		}

		f[0] = (1.0 + g)*x[0];
		f[1] = (1.0 + g)*(1- x[0]*x[0]);
	}


	void F3(const double * x,  int n, double * f, int m)
	{
		double pi = 3.1415926535897932384626433832795;
		double tj;
		double g = 0.0;

		for(int j = 1; j < n; j++)
		{
			tj = x[j] - sin(0.5*pi*x[0]);
			g+= 10.0*sin(pi*x[0]/2.0)*(abs(tj) / (1+exp(5.0*abs(tj) ) ) );
		}

		f[0] = (1.0 + g)*cos(pi*x[0]/2.0);
		f[1] = (1.0 + g)*sin(pi*x[0]/2.0);
	}


	void F4(const double * x,  int n, double * f, int m)
	{
		double pi = 3.1415926535897932384626433832795;
		double tj;
		double g = 0.0;

		for(int j = 1; j < n; j++)
		{
			tj = x[j] - sin(0.5*pi*x[0]);
			g+= 10.0*sin(pi*x[0])*(abs(tj) / (1+exp(5.0*abs(tj) ) ) );
		}
		//g+= 1.0;

		f[0] = (1.0 + g)*x[0];
		f[1] = (1.0 + g)*(1.0 - pow(x[0], 0.5)*cos(2*pi*x[0])*cos(2*pi*x[0]));
	}


	void F5(const double * x,  int n, double * f, int m)
	{
		double pi = 3.1415926535897932384626433832795;
		double tj;
		double g = 0.0;

		for(int j = 1; j < n; j++)
		{
			tj = x[j] - sin(0.5*pi*x[0]);
			g+= 2.0*abs(cos(pi*x[0]))*(-0.9*tj*tj + pow(abs(tj), 0.6));
		}

		f[0] = (1.0 + g)*x[0];
		f[1] = (1.0 + g)*(1- sqrt(x[0]));
	}


	void F6(const double * x,  int n, double * f, int m)
	{
		double pi = 3.1415926535897932384626433832795;
		double tj;
		double g = 0.0;

		for(int j = 2; j < n; j++)
		{
			tj = x[j] - x[0]*x[1];
			g+= 2.0*sin(pi*x[0])*(-0.9*tj*tj + pow(abs(tj), 0.6));
		}

		f[0] = (1.0 + g)*x[0]*x[1];
		f[1] = (1.0 + g)*x[0]*(1- x[1]);
		f[2] = (1.0 + g)*(1- x[0]);
	}


	void F7(const double * x,  int n, double * f, int m)
	{
		double pi = 3.1415926535897932384626433832795;
		double tj;
		double g = 0.0;

		for(int j = 2; j < n; j++)
		{
			tj = x[j] - x[0]*x[1];
			g+= 2.0*sin(pi*x[0])*(-0.9*tj*tj + pow(abs(tj), 0.6));
		}

		f[0] = (1.0 + g)*cos(pi*x[0]/2.0)*cos(pi*x[1]/2.0);
		f[1] = (1.0 + g)*cos(pi*x[0]/2.0)*sin(pi*x[1]/2.0);
		f[2] = (1.0 + g)*sin(pi*x[0]/2.0);
	}

}