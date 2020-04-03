
#ifndef LZ09_H
#define LZ09_H

#include <vector>
#include <cmath>



namespace LZ09
{
	void F1(const double * x, int n, double * f, int m);
	void F2(const double * x, int n, double * f, int m);
	void F3(const double * x, int n, double * f, int m);
	void F4(const double * x, int n, double * f, int m);
	void F5(const double * x, int n, double * f, int m);
	void F6(const double * x, int n, double * f, int m);
	void F7(const double * x, int n, double * f, int m);
	void F8(const double * x, int n, double * f, int m);
	void F9(const double * x, int n, double * f, int m);
}
#endif











// second version - defferent from matlab implementations

//#ifndef LZ09_H
//#define LZ09_H
//
//#include <vector>
//#include <cmath>
//
//
//namespace LZ09
//{
//
//	void alphafunction(double alpha[], std::vector<double>&x, int dim, int type);
//	double betafunction(std::vector <double>x, int type);
//	double psfunc2(double &x, double &t1, int dim, int type, int css);
//	double psfunc3(double &x, double &t1, double &t2, int dim, int type);
//	void LZ09functions(std::vector<double> &x_var, std::vector <double> &y_obj);
//
//	void F1(const double * x, int n, double * f, int m);
//	void F2(const double * x, int n, double * f, int m);
//	void F3(const double * x, int n, double * f, int m);
//	void F4(const double * x, int n, double * f, int m);
//	void F5(const double * x, int n, double * f, int m);
//	void F6(const double * x, int n, double * f, int m);
//	void F7(const double * x, int n, double * f, int m);
//	void F8(const double * x, int n, double * f, int m);
//	void F9(const double * x, int n, double * f, int m);
//}
//#endif