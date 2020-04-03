#include "Util.h"


CUtil::CUtil(void)
{
	// initialize random seed
	srand (time(NULL)); 
}


CUtil::~CUtil(void)
{
}


// function returns integer in [0, upper) 
int CUtil::randi(int upper)
{
	return (std::rand() % upper);
}


// function returns a uniform number in [0, 1]
double CUtil::rand(void)
{
	return (std::rand() / double(RAND_MAX));
}


// function returns a uniform number in [lower, upper]
double CUtil::rand(double lower, double upper)
{
    return (lower + rand()*(upper - lower));
}



void CUtil::minFastSort(double * x, int * idx, int n, int m)
{
	int i, j, id;
	double temp;

	for (i = 0; i < m; i++) 
	{
		for (j = i + 1; j < n; j++) 
		{
			if (x[i] > x[j]) 
			{
				temp = x[i];
				x[i] = x[j];
				x[j] = temp;
				id = idx[i];
				idx[i] = idx[j];
				idx[j] = id;
			} // if
		}
	} // for

} // minFastSort


void CUtil::randomPermutation(int * perm, int size) 
{
	int * index = new int[size];
	bool * flag = new bool[size];

	for (int n = 0; n < size; n++) 
	{
		index[n] = n;
		flag[n] = true;
	}

	int start;
	int num = 0;
	while (num < size)
	{
		start = randi(size); // random integer

		while (true) 
		{
			if (flag[start]) 
			{
				perm[num] = index[start];
				flag[start] = false;
				num++;
				break;
			}
			if (start == (size - 1)) {
				start = 0;
			} else {
				start++;
			}
		}
	} // while

	delete[] index;
	delete[] flag;

} // randomPermutation



double CUtil::distance(double * vector1, double * vector2, int dim)
{
	double sum = 0;
	for (int n = 0; n < dim; n++) {
		sum += (vector1[n] - vector2[n]) * (vector1[n] - vector2[n]);
	}
	return sqrt(sum);
} // distVector


// scalar product
double CUtil::innerProduct(double * vector1, double * vector2, int dim)
{
	double res = 0;
	for (int n = 0; n < dim; n++) {
		res += vector1[n] * vector2[n];
	}
	return res;
} // distVector


// cos(angle) between two vectors
double CUtil::cosAngle(double * vector1, double * vector2, int dim)
{
	double norm1 = 0.0;
    double norm2 = 0.0;
	double scalar = 0.0;

	for (int i = 0; i < dim; i++)
	{
		scalar += vector1[i] * vector2[i];
		norm1 += vector1[i]*vector1[i];
		norm2 += vector2[i]*vector2[i];
	}

	return scalar/sqrt(norm1*norm2);
} // cosAngle


// norm of the vector
double CUtil::Lpnorm(double * vector, int dim, int p)
{
	// Lp norm when p = inf, Chebyshev norm
	if(p==0)
	{
		double norm = vector[0];
		for (int i = 0; i < dim; i++)
		{
			if(vector[i] > norm){
				norm = vector[i];
			}
		}
		return norm;
	}

	// Lp norm
	double norm = 0.0;
	for (int i = 0; i < dim; i++){
		norm += pow(vector[i], p);
	}

	return pow(norm, 1.0/p);
}