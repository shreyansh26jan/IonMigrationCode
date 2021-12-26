#ifndef PostProcessing_h
#define PostProcessing_h

#include "Charge.h"

#include <eigen3/Eigen/Sparse>
#include <vector>

class PostProcessing
{
public:
  PostProcessing(std::vector<double> const &xi, Eigen::VectorXd const &si,
                 Eigen::VectorXd const &n, Eigen::VectorXd const &p, Eigen::VectorXd const &P);

  static double J0_1st_order(double x0, double x1, double x2, double J1, double J2);
  static double J0_2nd_order(double x0, double x1, double x2, double x3, double J1, double J2, double J3);

  std::vector<double> Ji(Charge::ChargeSign chargeSign) const;
  std::vector<double> xi_J() const;

  void plotCharges(double V_a) const;
  void plotCurrent(double V_a) const;
  void plotR(double V_a) const;
  static void plotJV(std::vector<double> const &J, std::vector<double> const &V,std::vector<double> const &justSi);

private:
  double integralfx(std::vector<double> const &f) const;

  std::vector<double> const &_xi;
  Eigen::VectorXd const &_si;
  Eigen::VectorXd const &_n;
  Eigen::VectorXd const &_p;
  Eigen::VectorXd const &_P;
};

#endif