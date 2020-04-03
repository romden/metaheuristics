
#include "DTLZ.h"
#include <math.h>

namespace DTLZ
{


	void DTLZ1(const double * x, int n, double * f, int m)
	{

		int k = n - m + 1; // number of position parameters
		double g = 0.0; // func

		double PI = 3.1415926535897932384626433832795;


		for (int i = n - k; i < n; i++){
			g += (x[i] - 0.5)*(x[i] - 0.5) - cos(20.0 * PI * ( x[i] - 0.5));
		}

		g = 100 * (k + g);
		for (int i = 0; i < m; i++){
			f[i] = (1.0 + g) * 0.5;
		}

		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < m - (i + 1); j++){
				f[i] *= x[j];
			}
			if (i != 0)
			{
				int aux = m - (i + 1);
				f[i] *= 1 - x[aux];
			} 
		}

	} // DTLZ1


	void DTLZ2(const double * x, int n, double * f, int m)
	{

		int k = n - m + 1; // number of position parameters
		double g = 0.0; // func

		double PI = 3.1415926535897932384626433832795;


		for (int i = n - k; i < n; i++){
			g += (x[i] - 0.5)*(x[i] - 0.5);
		}

		for (int i = 0; i < m; i++){
			f[i] = 1.0 + g;
		}

		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < m - (i + 1); j++){
				f[i] *= cos(x[j]*0.5*PI);
			}
			if (i != 0)
			{
				int aux = m - (i + 1);
				f[i] *= sin(x[aux]*0.5*PI);
			} //if
		}

	} // DTLZ2


	void DTLZ3(const double * x, int n, double * f, int m)
	{

		int k = n - m + 1; // number of position parameters
		double g = 0.0; // func

		double PI = 3.1415926535897932384626433832795;


		for (int i = n - k; i < n; i++){
			g += (x[i] - 0.5)*(x[i] - 0.5) - cos(20.0 * PI * (x[i] - 0.5));
		}

		g = 100.0 * (k + g);
		for (int i = 0; i < m; i++){
			f[i] = 1.0 + g;
		}

		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < m - (i + 1); j++){
				f[i] *= cos(x[j]*0.5*PI);
			}
			if (i != 0)
			{
				int aux = m - (i + 1);
				f[i] *= sin(x[aux]*0.5*PI);
			} 
		}

	} // DTLZ3


	void DTLZ4(const double * x, int n, double * f, int m)
	{

		int k = n - m + 1; // number of position parameters
		double g = 0.0; // func

		double PI = 3.1415926535897932384626433832795;


		double alpha = 100.0;
		for (int i = n - k; i < n; i++){
			g += (x[i] - 0.5)*(x[i] - 0.5);
		}

		for (int i = 0; i < m; i++){
			f[i] = 1.0 + g;
		}

		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < m - (i + 1); j++){
				f[i] *= cos(pow(x[j],alpha)*(PI/2.0));
			}
			if (i != 0)
			{
				int aux = m - (i + 1);
				f[i] *= sin(pow(x[aux],alpha)*(PI/2.0));
			} //if
		}

	} // DTLZ4


	void DTLZ5(const double * x, int n, double * f, int m)
	{

		int k = n - m + 1; // number of position parameters
		double g = 0.0; // func
		double PI = 3.1415926535897932384626433832795;
		double * theta = new double[m-1];
		double t;


		for (int i = n - k; i < n; i++){
			g += (x[i] - 0.5)*(x[i] - 0.5);
		}

		t = PI  / (4.0 * (1.0 + g));

		theta[0] = x[0] * PI / 2.0;
		for (int i = 1; i < (m-1); i++){
			theta[i] = t * (1.0 + 2.0 * g * x[i]);
		}

		for (int i = 0; i < m; i++){
			f[i] = 1.0 + g;
		}

		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < m - (i + 1); j++){
				f[i] *= cos(theta[j]);
			}
			if (i != 0)
			{
				int aux = m - (i + 1);
				f[i] *= sin(theta[aux]);
			} // if
		}

		delete [] theta;

	} // DTLZ5



	void DTLZ6(const double * x, int n, double * f, int m)
	{

		int k = n - m + 1; // number of position parameters
		double g = 0.0; // func
		double PI = 3.1415926535897932384626433832795;
		double * theta = new double[m-1];
		double t;


		for (int i = n - k; i < n; i++){
			g += pow(x[i],0.1);
		}

		t = PI  / (4.0 * (1.0 + g));

		theta[0] = x[0] * PI / 2.0;
		for (int i = 1; i < (m-1); i++){
			theta[i] = t * (1.0 + 2.0 * g * x[i]);
		}

		for (int i = 0; i < m; i++){
			f[i] = 1.0 + g;
		}

		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < m - (i + 1); j++){
				f[i] *= cos(theta[j]);
			}
			if (i != 0)
			{
				int aux = m - (i + 1);
				f[i] *= sin(theta[aux]);
			} 
		}

		delete [] theta;

	} // DTLZ6


	void DTLZ7(const double * x, int n, double * f, int m)
	{

		int k = n - m + 1; // number of position parameters
		double g = 0.0; // func

		double PI = 3.1415926535897932384626433832795;


		for (int i = n - k; i < n; i++){
			g += x[i] ;
		}

		g = 1.0 + (9.0 * g) / k;
		//Calculate the value of f1,f2,f3,...,fM-1 (take acount of vectors start at 0)
		for (int i = 0; i < m-1; i++){
			f[i] = x[i];
		}

		//Calculate fM
		double h = 0.0;
		for (int i = 0; i < m -1; i++){
			h += (f[i]/(1.0 + g))*(1.0 + sin(3.0 * PI * f[i]));
		}

		h = m - h;

		f[m-1] = (1.0 + g) * h;

	} // DTLZ7


}