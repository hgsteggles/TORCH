#include <iostream> //std::cout
#include <vector>


namespace NumericalRecipes {

std::vector<double> spline(const std::vector<double>& x, const std::vector<double>& f, double yp1, double ypn);
double splint(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& f2, double x2);
std::vector<double> vsplint(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& f2, const std::vector<double>& x2);

}
