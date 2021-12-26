#ifndef CurrentContinuityEquation_h
#define CurrentContinuityEquation_h

#include "Charge.h"

#include <eigen3/Eigen/Sparse>
#include <memory>
#include <vector>

class CurrentContinuityEquation
{
public:
  CurrentContinuityEquation(std::vector<double> const &xi, Eigen::VectorXd const &si,
                            Eigen::VectorXd const &n, Eigen::VectorXd const &p, Eigen::VectorXd const &P);

  void initialize(Charge::ChargeSign chargeSign, Eigen::VectorXd &sol, std::vector<double> &_si0, std::vector<double> &_n0, std::vector<double> &_p0, std::vector<double> &_P0);
  double solve(Charge::ChargeSign chargeSign, Eigen::VectorXd &sol, std::vector<double> &_si0, std::vector<double> &_n0, std::vector<double> &_p0, std::vector<double> &_P0);

private:
  std::vector<double> const &_xi;
  Eigen::VectorXd const &_si;
  Eigen::VectorXd const &_n;
  Eigen::VectorXd const &_p;
  Eigen::VectorXd const &_P;

  std::vector<Eigen::Triplet<double> > _tripletList;
  Eigen::SparseMatrix<double> _mat;
  Eigen::VectorXd _vec;
  Eigen::VectorXd _res;
  Eigen::VectorXd _dsol;

  std::unique_ptr<Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > > _solver;
};

#endif