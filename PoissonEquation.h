#ifndef PoissonEquation_h
#define PoissonEquation_h

#include <eigen3/Eigen/Sparse>
#include <memory>
#include <vector>

class PoissonEquation
{
public:
  PoissonEquation(std::vector<double> const &xi, Eigen::VectorXd const &n, Eigen::VectorXd const &p, Eigen::VectorXd const &P);

  void initialize(double V_a, Eigen::VectorXd &sol);
  double solve(double V_a, Eigen::VectorXd &sol);

private:
  std::vector<double> const &_xi;
  Eigen::VectorXd const &_n;
  Eigen::VectorXd const &_p;
  Eigen::VectorXd const &_P;

  std::vector<Eigen::Triplet<double> > _tripletList;
  Eigen::VectorXd _vec;
  Eigen::VectorXd _dsol;

  Eigen::SparseMatrix<double> _mat;
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > _solver;
};

#endif