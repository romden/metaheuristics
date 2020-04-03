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


// function to calculate Lp norm
double CUtil::Lp_norm(double * v, int m, double p)
{
	double	norm = 0;
	for(int j = 0; j < m; j++){		
		norm += pow(v[j], p); 
	}
	return pow(norm, 1.0/p);
}


// function to calculate Euclidean distance between two points
double CUtil::distance(double * point1, double * point2, int dim)
{
	double sum = 0;
	for (int n = 0; n < dim; n++) {
		sum += (point1[n] - point2[n]) * (point1[n] - point2[n]);
	}
	return sqrt(sum);
} // distance


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