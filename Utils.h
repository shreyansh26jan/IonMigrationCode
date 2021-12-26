#ifndef Utils_h
#define Utils_h

#include <cmath>
#include <vector>

namespace Utils
{

std::vector<double> LinearSpacedArray(double a, double b, std::size_t N)
{
  const double h = (b - a) / static_cast<double>(N-1);
  std::vector<double> xs(N);
  std::vector<double>::iterator x;
  double val;
  for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
    *x = val;
  }
  return xs;
}

// grid clustering increased rms error in n, J

double cluster1sided(double x, double L, double b)
{
  const double r = std::pow((b+1.)/(b-1.), 1.-x/L);
  return L * ((b+1.) - (b-1.)*r) / (r+1.);
}

double cluster2sided(double x, double L, double b)
{
  const double r = std::pow((b+1.)/(b-1.), 2.*x/L-1.);
  return L * (-(b-1.) + (b+1.)*r) / (2.*(r+1.));
}

}

#endif