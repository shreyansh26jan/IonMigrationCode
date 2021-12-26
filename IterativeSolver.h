#ifndef IterativeSolver_h
#define IterativeSolver_h

#include "CurrentContinuityEquation.h"
#include "PoissonEquation.h"

#include <eigen3/Eigen/Sparse>
#include <memory>
#include <vector>

struct Solution
{
  double J_0;
  double J_1;
  double dsol_l2norm_sum;
  
};

class IterativeSolver
{
public:
  explicit IterativeSolver(std::vector<double> const &xi);
  Solution solve(double V_a, std::vector<double> &_si0, std::vector<double> &_n0, std::vector<double> &_p0, std::vector<double> &_P0  );//std::vector<double> &_si0, std::vector<double> &_n0,std::vector<double> &_p0, std::vector<double> &_P0

private:
  std::vector<double> const &_xi;
  Eigen::VectorXd _n;
  Eigen::VectorXd _p;
  Eigen::VectorXd _P;
  Eigen::VectorXd _si;

  std::unique_ptr<CurrentContinuityEquation> _ccEqn;
  std::unique_ptr<CurrentContinuityEquation> _ccEqn1;
  std::unique_ptr<PoissonEquation> _pEqn;
};

#endif