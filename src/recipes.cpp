#include "recipes.h"

/**
 * @brief Calculates second derivatives of an interpolating function.
 * @param x Vector of function arguments.
 * @param f Vector of function values given a vector of function arguments.
 * @param yp1 First derivative of function at x[0].
 * @param ypn First derivative of function at x[n-1].
 * @return Vector of second derivatives of interpolating function.
 */
std::vector<double> NumericalRecipes::spline(const std::vector<double>& x, const std::vector<double>& f, double yp1, double ypn){
	int n = x.size();
	std::vector<double> y2(n, 0.0);
	std::vector<double> u(n-1, 0.0);
	int i,k;
	double p,qn,sig,un;
	//The lower boundary condition is set either to be "natural"
	//or else to have a specified first derivative.
	if(yp1 > 0.99e30){
		y2[0] = 0.0;
		u[0] = 0.0;
	}
	else{
		y2[0] = -0.5;
		u[0] = (3.0/(x[1]-x[0]))*((f[1]-f[0])/(x[1]-x[0])-yp1);
	}
	//This is the decomposition loop of the tridiagonal algorithm.
	//y2 and u are used for temporary storage of the decomposed factors.
	for(i = 1; i < n-1; i++){
		sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p = sig*y2[i-1]+2.0;
		y2[i] = (sig-1.0)/p;
		u[i] = (f[i+1]-f[i])/(x[i+1]-x[i]) - (f[i]-f[i-1])/(x[i]-x[i-1]);
		u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	//The upper boundary condition is set either to be "natural"
	// or else to have a specified first derivative.
	if (ypn > 0.99e30){
		qn = 0.0;
		un = 0.0;
	}
	else{
		qn = 0.5;
		un = (3.0/(x[n-1]-x[n-2]))*(ypn-(f[n-1]-f[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	//This is the backsubstitution loop of the tridiagonal algorithm.
	for(k = n-2; k >= 0; k--)
		y2[k] = y2[k]*y2[k+1]+u[k];
	return y2;
}

/**
 * @brief Calculates a cubic spline interpolated value from a function, and its second derivative.
 * @param x Vector of function args.
 * @param f Vector of function values given a vector of function arguments.
 * @param f2 Second derivative of function.
 * @param x2 Interpolation location.
 * @return Cubic spline interpolated value.
 */
double NumericalRecipes::splint(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& f2, double x2){
	int klo, khi, k;
	double h,b,a;
	int n = x.size();
	//We will find the right place in the table by means of bisection. This is
	//optimal if sequential calls to this routine are at random values of x.
	//If sequential calls are in order, and closely spaced, one would do better
	//to store previous values of klo and khi and test if they remain appropriate
	//on the next call.
	/*
	cerr << "x = {" << x[0];
	for(unsigned int i = 1; i < x.size(); i++)
		cerr << ", " << x[i];
	cerr << "}" << endl;
	cerr << "x2 = " << x2 << endl;
	*/
	klo = 0;
	khi = n-1;
	while(khi-klo > 1){
		k = (khi+klo) >> 1;
		if(x[k] > x2)
			khi = k;
		else
			klo = k;
	}
	//cerr << "klo = " << klo << '\t' << "khi = " << khi << endl;
	//klo and khi now bracket the input value of x.
	h = x[khi]-x[klo];
	if(h == 0.0){
		std::cerr << "ERROR: Bad x input to routine splint()" << std::endl; //The xa's must be distinct.
		for(int i = 0; i < n; i++)
			std::cerr << "x[" << i << "] = " << x[i];
		exit(EXIT_FAILURE);
	}
	a = (x[khi]-x2)/h;
	b = (x2-x[klo])/h; //Cubic spline polynomial is now evaluated.

	return a*f[klo]+b*f[khi]+((a*a*a-a)*f2[klo]+(b*b*b-b)*f2[khi])*(h*h)/6.0;
}

/**
 * @brief Calculates cubic spline interpolated values from a function, and its second derivative.
 * @param x Vector of function args.
 * @param f Vector of function values given a vector of function arguments.
 * @param f2 Second derivative of function.
 * @param x2 Interpolation locations.
 * @return Vector of cubic spline interpolated values.
 */
std::vector<double> NumericalRecipes::vsplint(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& f2, const std::vector<double>& x2){
	std::vector<double> result(x2.size(), 0.0);
	for(unsigned int i = 0; i < x2.size(); i++)
		result[i] = splint(x, f, f2, x2[i]);
	return result;
}
