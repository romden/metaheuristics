#include "LZ09.h"


namespace LZ09
{

	const double PI = 3.1415926535897932384626433832795;


	void F1(const double * x, int n, double * f, int m)
	{

		unsigned int j, count1, count2;
		double sum1, sum2, yj;

		sum1 = sum2 = 0.0;
		count1 = count2 = 0;

		for ( j = 2; j <= n; j++ )
		{
			yj = x[j - 1] - pow(x[0],  0.5*(1.0 + 3.0*(j - 2.0)/(n - 2.0) ));
			yj = yj * yj;
			if ( j % 2 == 0 )
			{
				sum2 += yj;
				count2++;
			}
			else
			{
				sum1 += yj;
				count1++;
			}
		}

		f[0] = x[0] + 2.0 * sum1 / ( double )count1;
		f[1] = 1.0 - sqrt( x[0] ) + 2.0 * sum2 / ( double )count2;

	} // LZ09_F1



	void F2(const double * x, int n, double * f, int m)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;

		sum1 = sum2 = 0.0;
		count1 = count2 = 0;
		for ( j = 2; j <= n; j++ )
		{
			yj = x[j - 1] - sin( 6.0 * PI * x[0] + j * PI / n );
			yj = yj * yj;
			if ( j % 2 == 0 )
			{
				sum2 += yj;
				count2++;
			}
			else
			{
				sum1 += yj;
				count1++;
			}
		}
		f[0] = x[0] + 2.0 * sum1 / ( double )count1;
		f[1] = 1.0 - sqrt( x[0] ) + 2.0 * sum2 / ( double )count2;

	} // LZ09_F2


	void F3(const double * x, int n, double * f, int m)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;

		sum1 = sum2 = 0.0;
		count1 = count2 = 0;
		for ( j = 2; j <= n; j++ )
		{ 
			if ( j % 2 == 0 )
			{
				yj = x[j - 1] - 0.8*x[0]*sin( 6.0 * PI * x[0] + j * PI / n );
				sum2 += yj * yj;
				count2++;
			}
			else
			{
				yj = x[j - 1] - 0.8*x[0]*cos( 6.0 * PI * x[0] + j * PI / n );
				sum1 += yj * yj;
				count1++;
			}
		}
		f[0] = x[0] + 2.0 * sum1 / ( double )count1;
		f[1] = 1.0 - sqrt( x[0] ) + 2.0 * sum2 / ( double )count2;

	} // LZ09_F3


	void F4(const double * x, int n, double * f, int m)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;

		sum1 = sum2 = 0.0;
		count1 = count2 = 0;
		for ( j = 2; j <= n; j++ )
		{ 
			if ( j % 2 == 0 )
			{
				yj = x[j - 1] - 0.8*x[0]*sin( 6.0 * PI * x[0] + j * PI / n );
				sum2 += yj * yj;
				count2++;
			}
			else
			{
				yj = x[j - 1] - 0.8*x[0]*cos( (6.0 * PI * x[0] + j * PI / n)/3 );
				sum1 += yj * yj;
				count1++;
			}
		}
		f[0] = x[0] + 2.0 * sum1 / ( double )count1;
		f[1] = 1.0 - sqrt( x[0] ) + 2.0 * sum2 / ( double )count2;

	} // LZ09_F4




	void F5(const double * x, int n, double * f, int m)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;

		sum1 = sum2 = 0.0;
		count1 = count2 = 0;
		for ( j = 2; j <= n; j++ )
		{
			if ( j % 2 == 0 )
			{
				yj = x[j - 1] - 0.3 * x[0] * ( x[0] * cos( 24.0 * PI * x[0] + 4.0 * j * PI/n ) + 2.0 )
					* sin( 6.0 * PI * x[0] + j * PI/n );
				sum2 += yj * yj;
				count2++;
			}
			else
			{
				yj = x[j - 1] - 0.3 * x[0] * ( x[0] * cos( 24.0 * PI * x[0] + 4.0 * j * PI/n ) + 2.0 )
					* cos( 6.0 * PI * x[0] + j * PI/n );
				sum1 += yj * yj;
				count1++;
			}
		}
		f[0] = x[0] + 2.0 * sum1 / ( double )count1;
		f[1] = 1.0 - sqrt( x[0] ) + 2.0 * sum2 / ( double )count2;

	} // LZ09_F5


	void F6(const double * x, int n, double * f, int m)
	{
		unsigned int j, count1, count2, count3;
		double sum1, sum2, sum3, yj;

		sum1 = sum2 = sum3 = 0.0;
		count1 = count2 = count3 = 0;
		for ( j = 3; j <= n; j++ )
		{
			yj = x[j - 1] - 2.0 * x[1] * sin( 2.0 * PI * x[0] + j * PI / n );
			if ( j % 3 == 1 )
			{
				sum1 += yj * yj;
				count1++;
			}
			else if ( j % 3 == 2 )
			{
				sum2 += yj * yj;
				count2++;
			}
			else
			{
				sum3 += yj * yj;
				count3++;
			}
		}
		f[0] = cos( 0.5 * PI * x[0] ) * cos( 0.5 * PI * x[1] ) + 2.0 * sum1 / ( double )count1;
		f[1] = cos( 0.5 * PI * x[0] ) * sin( 0.5 * PI * x[1] ) + 2.0 * sum2 / ( double )count2;
		f[2] = sin( 0.5 * PI * x[0] ) + 2.0 * sum3 / ( double )count3;

	} // LZ09_F6



	void F7(const double * x, int n, double * f, int m)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, prod1, prod2, yj, pj;

		sum1 = sum2 = 0.0;
		count1 = count2 = 0;
		prod1 = prod2 = 1.0;
		for ( j = 2; j <= n; j++ )
		{
			yj = x[j - 1] - pow( x[0], 0.5 * ( 1.0 + 3.0 * ( j - 2.0 ) / ( n - 2.0 ) ) );
			pj = cos( 8.0 * yj * PI );
			if ( j % 2 == 0 )
			{
				sum2 += 4.0*yj * yj - pj +1.0;
				count2++;
			}
			else
			{
				sum1 += 4.0*yj * yj - pj +1.0;
				count1++;
			}
		}
		f[0] = x[0] + 2.0 * sum1 / ( double )count1;
		f[1] = 1.0 - sqrt( x[0] ) + 2.0 * sum2 / ( double )count2;

	} // LZ09_F7



	void F8(const double * x, int n, double * f, int m)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, prod1, prod2, yj, pj;

		sum1 = sum2 = 0.0;
		count1 = count2 = 0;
		prod1 = prod2 = 1.0;
		for ( j = 2; j <= n; j++ )
		{
			yj = x[j - 1] - pow( x[0], 0.5 * ( 1.0 + 3.0 * ( j - 2.0 ) / ( n - 2.0 ) ) );
			pj = cos( 20.0 * yj * PI / sqrt( j + 0.0 ) );
			if ( j % 2 == 0 )
			{
				sum2 += yj * yj;
				prod2 *= pj;
				count2++;
			}
			else
			{
				sum1 += yj * yj;
				prod1 *= pj;
				count1++;
			}
		}
		f[0] = x[0] + 2.0 * ( 4.0 * sum1 - 2.0 * prod1 + 2.0 ) / ( double )count1;
		f[1] = 1.0 - sqrt( x[0] ) + 2.0 * ( 4.0 * sum2 - 2.0 * prod2 + 2.0 ) / ( double )count2;

	} // LZ09_F8


	void F9(const double * x, int n, double * f, int m)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;

		sum1 = sum2 = 0.0;
		count1 = count2 = 0;
		for ( j = 2; j <= n; j++ )
		{
			yj = x[j - 1] - sin( 6.0 * PI * x[0] + j * PI / n );
			yj = yj * yj;
			if ( j % 2 == 0 )
			{
				sum2 += yj;
				count2++;
			}
			else
			{
				sum1 += yj;
				count1++;
			}
		}
		f[0] = x[0] + 2.0 * sum1 / ( double )count1;
		f[1] = 1.0 - x[0]*x[0] + 2.0 * sum2 / ( double )count2;

	} // LZ09_F


} // LZ09 namespace
















//namespace LZ09
//{
//	int ptype;
//	int ltype;
//	int dtype;
//
//	int nvar;
//	int nobj;
//
//	void F1(const double * x, int n, double * f, int m)
//	{
//
//		ptype = 21;
//		ltype = 21;
//		dtype = 1;
//
//		nvar = n;
//		nobj = m;
//
//		std::vector<double>x_var(x, x+n);
//		std::vector<double>y_obj(f, f+m);
//
//		LZ09functions(x_var, y_obj);
//
//		for(std::vector<double>::iterator it = y_obj.begin(); it != y_obj.end(); ++it){
//			f[it - y_obj.begin()] = *it;
//		}
//
//	} // LZ09_F1
//
//
//	void F2(const double * x, int n, double * f, int m)
//	{
//
//		ptype = 21;
//		ltype = 22;
//		dtype = 1;
//
//		nvar = n;
//		nobj = m;
//
//		std::vector<double>x_var(x, x+n);
//		std::vector<double>y_obj(f, f+m);
//
//		LZ09functions(x_var, y_obj);
//
//		for(std::vector<double>::iterator it = y_obj.begin(); it != y_obj.end(); ++it){
//			f[it - y_obj.begin()] = *it;
//		}
//
//	} // LZ09_F2
//
//
//	void F3(const double * x, int n, double * f, int m)
//	{
//
//		ptype = 21;
//		ltype = 23;
//		dtype = 1;
//
//		nvar = n;
//		nobj = m;
//
//		std::vector<double>x_var(x, x+n);
//		std::vector<double>y_obj(f, f+m);
//
//		LZ09functions(x_var, y_obj);
//
//		for(std::vector<double>::iterator it = y_obj.begin(); it != y_obj.end(); ++it){
//			f[it - y_obj.begin()] = *it;
//		}
//
//	} // LZ09_F3
//
//
//	void F4(const double * x, int n, double * f, int m)
//	{
//
//		ptype = 21;
//		ltype = 24;
//		dtype = 1;
//
//		nvar = n;
//		nobj = m;
//
//		std::vector<double>x_var(x, x+n);
//		std::vector<double>y_obj(f, f+m);
//
//		LZ09functions(x_var, y_obj);
//
//		for(std::vector<double>::iterator it = y_obj.begin(); it != y_obj.end(); ++it){
//			f[it - y_obj.begin()] = *it;
//		}
//
//	} // LZ09_F4
//
//
//	void F5(const double * x, int n, double * f, int m)
//	{
//
//		ptype = 21;
//		ltype = 26;
//		dtype = 1;
//
//		nvar = n;
//		nobj = m;
//
//		std::vector<double>x_var(x, x+n);
//		std::vector<double>y_obj(f, f+m);
//
//		LZ09functions(x_var, y_obj);
//
//		for(std::vector<double>::iterator it = y_obj.begin(); it != y_obj.end(); ++it){
//			f[it - y_obj.begin()] = *it;
//		}
//
//	} // LZ09_F5
//
//
//	void F6(const double * x, int n, double * f, int m)
//	{
//
//		ptype = 31;
//		ltype = 32;
//		dtype = 1;
//
//		nvar = n;
//		nobj = m;
//
//		std::vector<double>x_var(x, x+n);
//		std::vector<double>y_obj(f, f+m);
//
//		LZ09functions(x_var, y_obj);
//
//		for(std::vector<double>::iterator it = y_obj.begin(); it != y_obj.end(); ++it){
//			f[it - y_obj.begin()] = *it;
//		}
//
//	} // LZ09_F6
//
//
//	void F7(const double * x, int n, double * f, int m)
//	{
//
//		ptype = 21;
//		ltype = 21;
//		dtype = 3;
//
//		nvar = n;
//		nobj = m;
//
//		std::vector<double>x_var(x, x+n);
//		std::vector<double>y_obj(f, f+m);
//
//		LZ09functions(x_var, y_obj);
//
//		for(std::vector<double>::iterator it = y_obj.begin(); it != y_obj.end(); ++it){
//			f[it - y_obj.begin()] = *it;
//		}
//
//	} // LZ09_F7
//
//
//	void F8(const double * x, int n, double * f, int m)
//	{
//
//		ptype = 21;
//		ltype = 21;
//		dtype = 4;
//
//		nvar = n;
//		nobj = m;
//
//		std::vector<double>x_var(x, x+n);
//		std::vector<double>y_obj(f, f+m);
//
//		LZ09functions(x_var, y_obj);
//
//		for(std::vector<double>::iterator it = y_obj.begin(); it != y_obj.end(); ++it){
//			f[it - y_obj.begin()] = *it;
//		}
//
//	} // LZ09_F8
//
//
//	void F9(const double * x, int n, double * f, int m)
//	{
//
//		ptype = 22;
//		ltype = 22;
//		dtype = 1;
//
//		nvar = n;
//		nobj = m;
//
//		std::vector<double>x_var(x, x+n);
//		std::vector<double>y_obj(f, f+m);
//
//		LZ09functions(x_var, y_obj);
//
//		for(std::vector<double>::iterator it = y_obj.begin(); it != y_obj.end(); ++it){
//			f[it - y_obj.begin()] = *it;
//		}
//
//	} // LZ09_F9
//
//
//
//	// control the PF shape
//	void alphafunction(double alpha[], std::vector<double>&x, int dim, int type)
//	{
//		double pi = 3.1415926535897932384626433832795;
//
//		if(dim==2)
//		{
//			if(type==21){
//				alpha[0] = x[0];
//				alpha[1] = 1 - sqrt(x[0]);
//			}
//
//			if(type==22){
//				alpha[0] = x[0];
//				alpha[1] = 1 - x[0]*x[0];
//			}
//
//			if(type==23){
//				alpha[0] = x[0];
//				alpha[1] = 1 - sqrt(alpha[0]) - alpha[0]*sin(10*alpha[0]*alpha[0]*pi);
//			}
//
//			if(type==24){
//				alpha[0] = x[0];
//				alpha[1] = 1 - x[0] - 0.05*sin(4*pi*x[0]);
//			}
//
//		}
//		else
//		{
//
//			if(type==31){
//				alpha[0] = cos(x[0]*pi/2)*cos(x[1]*pi/2);
//				alpha[1] = cos(x[0]*pi/2)*sin(x[1]*pi/2);
//				alpha[2] = sin(x[0]*pi/2);		
//			}
//
//			if(type==32){
//				alpha[0] = 1 - cos(x[0]*pi/2)*cos(x[1]*pi/2);
//				alpha[1] = 1 - cos(x[0]*pi/2)*sin(x[1]*pi/2);
//				alpha[2] = 1 - sin(x[0]*pi/2);		
//			}
//
//			if(type==33){
//				alpha[0] = x[0];
//				alpha[1] = x[1];
//				alpha[2] = 3 - (sin(3*pi*x[0]) + sin(3*pi*x[1])) - 2*(x[0] + x[1]);
//			}
//
//			if(type==34){
//				alpha[0] = x[0]*x[1];
//				alpha[1] = x[0]*(1 - x[1]);
//				alpha[2] = (1 - x[0]);		
//			}	
//		}
//
//	} // alphafunction
//
//
//	// control the distance
//	double betafunction(std::vector <double>x, int type)
//	{
//		double pi = 3.1415926535897932384626433832795;
//		double beta;
//		int dim = x.size();
//
//		if(dim==0) beta = 0;
//
//		if(type==1){
//			beta = 0;
//			for(int i=0; i<dim; i++){
//				beta+= x[i]*x[i];
//			}	   
//			beta = 2.0*beta/dim;
//		}
//
//		if(type==2){
//			beta = 0;
//			for(int i=0; i<dim; i++){
//				beta+= sqrt((double)(i+1))*x[i]*x[i];
//			}	   
//			beta = 2.0*beta/dim;
//		}
//
//		if(type==3){
//			double sum = 0, xx;
//			for(int i=0; i<dim; i++){
//				xx = 2*x[i];
//				sum+= (xx*xx - cos(4*pi*xx) + 1);			
//			}	
//			beta = 2.0*sum/dim;
//		}
//
//		if(type==4){
//			double sum = 0, prod = 1, xx;
//			for(int i=0; i<dim; i++){
//				xx  = 2*x[i];
//				sum+= xx*xx;
//				prod*=cos(10*pi*xx/sqrt((double)(i+1)));
//			}	    		
//			beta = 2.0*(sum - 2*prod + 2)/dim;	
//		}
//
//		return beta;
//
//	} // betafunction
//
//
//	// control the PS shape of 2-d instances
//	double psfunc2(double &x, double &t1, int dim, int type, int css){
//		// type:  the type of curve 
//		// css:   the class of index
//
//		double pi = 3.1415926535897932384626433832795;
//		double beta;
//
//		dim++;
//
//		if(type==21){
//			double xy   = 2*(x - 0.5);
//			beta = xy - pow(t1, 0.5*(nvar + 3*dim - 8)/(nvar - 2));
//		}	
//
//		if(type==22){
//			double theta = 6*pi*t1 + dim*pi/nvar;  
//			double xy    = 2*(x - 0.5);
//			beta = xy - sin(theta);
//		}
//
//		if(type==23){
//			double theta = 6*pi*t1 + dim*pi/nvar;
//			double ra    = 0.8*t1;
//			double xy    = 2*(x - 0.5);
//			if(css==1)
//				beta = xy - ra*cos(theta);
//			else{
//				beta = xy - ra*sin(theta);
//			}
//		}
//
//		if(type==24){
//			double theta = 6*pi*t1 + dim*pi/nvar;
//			double xy    = 2*(x - 0.5);
//			double ra    = 0.8*t1;
//			if(css==1)
//				beta = xy - ra*cos(theta/3);
//			else{
//				beta = xy - ra*sin(theta);
//			}
//		}
//
//		if(type==25){
//			double rho   = 0.8;
//			double phi   = pi*t1;
//			double theta = 6*pi*t1 + dim*pi/nvar;
//			double xy    = 2*(x - 0.5);
//			if(css==1)
//				beta = xy - rho*sin(phi)*sin(theta);
//			else if(css==2)
//				beta = xy - rho*sin(phi)*cos(theta);
//			else
//				beta = xy - rho*cos(phi);			
//		}
//
//		if(type==26){
//			double theta = 6*pi*t1 + dim*pi/nvar;
//			double ra    = 0.3*t1*(t1*cos(4*theta) + 2);
//			double xy    = 2*(x - 0.5);
//			if(css==1)
//				beta = xy - ra*cos(theta);
//			else{
//				beta = xy - ra*sin(theta);
//			}
//		}
//
//		return beta;
//
//	} // psfunc2
//
//
//	// control the PS shapes of 3-D instances
//	double psfunc3(double &x, double &t1, double &t2, int dim, int type){
//		// type:  the type of curve 
//		// css:   the class of index
//
//		double pi = 3.1415926535897932384626433832795;
//		double beta;
//
//		dim++;
//
//		if(type==31){
//			double xy  = 4*(x - 0.5);
//			double rate = 1.0*dim/nvar;
//			beta = xy - 4*(t1*t1*rate + t2*(1.0-rate)) + 2;
//		}
//
//		if(type==32){
//			double theta = 2*pi*t1 + dim*pi/nvar;
//			double xy    = 4*(x - 0.5);
//			beta = xy - 2*t2*sin(theta);	
//		}
//
//		return beta;
//
//	} // psfunc3
//
//
//	void LZ09functions(std::vector<double> &x_var, std::vector <double> &y_obj)
//	{
//
//		// 2-objective case
//		if(nobj==2)
//		{
//			if(ltype==21||ltype==22||ltype==23||ltype==24||ltype==26)
//			{
//				double g = 0, h = 0, a, b;
//				std::vector <double> aa;
//				std::vector <double> bb;
//				for(int n=1;n<nvar;n++)
//				{
//
//					if(n%2==0){
//						a = psfunc2(x_var[n],x_var[0],n,ltype,1);  // linkage
//						aa.push_back(a);
//					}
//					else
//					{
//						b = psfunc2(x_var[n],x_var[0],n,ltype,2);
//						bb.push_back(b);
//					}	
//
//				}
//
//				g = betafunction(aa,dtype);
//				h = betafunction(bb,dtype);
//
//				double alpha[2];
//				alphafunction(alpha,x_var,2,ptype);  // shape function
//				y_obj[0] = alpha[0] + h;
//				y_obj[1] = alpha[1] + g; 
//				aa.clear(); 
//				bb.clear();
//			}
//
//			if(ltype==25)
//			{
//				double g = 0, h = 0, a, b;
//				double e = 0, c;
//				std::vector <double> aa;
//				std::vector <double> bb;
//				for(int n=1;n<nvar;n++){
//					if(n%3==0){
//						a = psfunc2(x_var[n],x_var[0],n,ltype,1); 
//						aa.push_back(a);
//					}
//					else if(n%3==1)
//					{
//						b = psfunc2(x_var[n],x_var[0],n,ltype,2);
//						bb.push_back(b);
//					}	
//					else{
//						c = psfunc2(x_var[n],x_var[0],n,ltype,3);
//						if(n%2==0)    aa.push_back(c);			
//						else          bb.push_back(c);
//					}
//				}		
//				g = betafunction(aa,dtype);          // distance function
//				h = betafunction(bb,dtype);
//				double alpha[2];
//				alphafunction(alpha,x_var,2,ptype);  // shape function
//				y_obj[0] = alpha[0] + h;
//				y_obj[1] = alpha[1] + g; 
//				aa.clear(); 
//				bb.clear();
//			}
//		}
//
//
//		// 3-objective case
//		if(nobj==3)
//		{
//			if(ltype==31||ltype==32)
//			{
//				double g = 0, h = 0, e = 0, a;
//				std::vector <double> aa;
//				std::vector <double> bb;
//				std::vector <double> cc;
//				for(int n=2;n<nvar;n++)
//				{
//					a = psfunc3(x_var[n],x_var[0],x_var[1],n,ltype);
//					if(n%3==0)	    aa.push_back(a);
//					else if(n%3==1)	bb.push_back(a);
//					else            cc.push_back(a);
//				}
//
//				g = betafunction(aa,dtype);
//				h = betafunction(bb,dtype);
//				e = betafunction(cc,dtype);
//
//				double alpha[3];
//				alphafunction(alpha,x_var,3,ptype);  // shape function
//				y_obj[0] = alpha[0] + h;
//				y_obj[1] = alpha[1] + g; 
//				y_obj[2] = alpha[2] + e; 
//				aa.clear(); 
//				bb.clear();
//				cc.clear();
//			}
//		}
//
//	} // LZ09functions
//
//
//}// LZ09 namespace
