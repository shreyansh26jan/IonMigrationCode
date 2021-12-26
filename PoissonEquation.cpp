#include "PoissonEquation.h"
#include "BoundaryConditions.h"
#include "DebugMessages.h"

#include <iostream>

PoissonEquation::
PoissonEquation(std::vector<double> const &xi, Eigen::VectorXd const &n, Eigen::VectorXd const &p, Eigen::VectorXd const &P)
: _xi(xi)
, _n(n)
, _p(p)
, _P(P)
{
  const int Ni = _xi.size();
  _tripletList.resize((Ni-2)*3+2);
  _vec.resize(Ni);
  _dsol.resize(Ni);
  _mat.resize(Ni, Ni);
}


void
PoissonEquation::
initialize(double V_a, Eigen::VectorXd &sol)
{
  const int Ni = _xi.size();
  sol.setLinSpaced(Ni, BoundaryConditions::si1(V_a), BoundaryConditions::si2(V_a));
}


double
PoissonEquation::
solve(double V_a, Eigen::VectorXd &sol)
{
  const double si1 = BoundaryConditions::si1(V_a);
  const double si2 = BoundaryConditions::si2(V_a);

  // Solution using SG method
  const int Ni = _xi.size();

  // assembly
  typedef Eigen::Triplet<double> T;

  // boundary conditions
  int tripletListIndex = 0;
  _tripletList[tripletListIndex++] = T(0, 0, 1);
  _tripletList[tripletListIndex++] = T(Ni-1, Ni-1, 1);
  _vec[0] = 0;
  _vec[Ni-1] = 0;

  // fill A in Ax=b
  for(int i = 1; i < Ni-1; i++)
  {
    const double x0 = Normalization::denormalized_x(_xi[i-1]);
    const double x1 = Normalization::denormalized_x(_xi[i]);
    const double x2 = Normalization::denormalized_x(_xi[i+1]);
    const double dxi = Parameters::lambda2 * 2. / (_xi[i+1] - _xi[i-1]);

    const double eps_plus  = 0.5 * (Parameters::eps(x1) + Parameters::eps(x2)) / Parameters::eps_max;
    const double eps_minus = 0.5 * (Parameters::eps(x0) + Parameters::eps(x1)) / Parameters::eps_max;

    // i - 1
    const double c0 = dxi * eps_minus / (_xi[i] - _xi[i-1]);
    _tripletList[tripletListIndex++] = T( i, i-1, c0 );

    // i
    const double c1 = (-dxi * eps_plus  / (_xi[i+1] - _xi[i])
                       -dxi * eps_minus / (_xi[i] - _xi[i-1]))
                       -_n[i] - _p[i] -_P[i];                      //check for P
    _tripletList[tripletListIndex++] = T( i, i, c1 ); 

    // i + 1
    const double c2 = dxi * eps_plus / (_xi[i+1] - _xi[i]);
    _tripletList[tripletListIndex++] = T( i, i+1, c2 );
  }

#ifdef printInnerConvergence
  std::cout << "potential [si]:" << std::endl;
#endif

  _mat.setFromTriplets(_tripletList.begin(), _tripletList.end());

  // Compute the ordering permutation vector from the structural pattern of A
  _solver.analyzePattern(_mat);
  // Compute the numerical factorization
  _solver.factorize(_mat);

  constexpr int iter_max = 0;
  constexpr double dsol_abs_tol = 1e-20;
  constexpr double URF = 1;
  double dsol_l2norm_sum = 0;
  {
    // solve \delta x in A\delta x = b - Ax_n
    int iter = 0,iteration=0;
    do
    {
      //std::cout << "iteration = " << iteration  << std::endl;
      // fill b in Ax=b
      for(int i = 1; i < Ni-1; i++)
      {
        const double x0 = Normalization::denormalized_x(_xi[i-1]);
        const double x1 = Normalization::denormalized_x(_xi[i]);
        const double x2 = Normalization::denormalized_x(_xi[i+1]);
        const double dxi = Parameters::lambda2 * 2. / (_xi[i+1] - _xi[i-1]);

        const double eps_plus  = 0.5 * (Parameters::eps(x1) + Parameters::eps(x2)) / Parameters::eps_max;
        const double eps_minus = 0.5 * (Parameters::eps(x0) + Parameters::eps(x1)) / Parameters::eps_max;

        // i - 1
        const double c0 = dxi * eps_minus / (_xi[i] - _xi[i-1]);
        // i
        const double c1 = (-dxi * eps_plus  / (_xi[i+1] - _xi[i])
                           -dxi * eps_minus / (_xi[i] - _xi[i-1]));
        // i + 1
        const double c2 = dxi * eps_plus / (_xi[i+1] - _xi[i]);

        const double D = Normalization::normalized_n(Parameters::N_D_etl * Parameters::is_etl(x1)) -
                         Normalization::normalized_n(Parameters::N_A_htl * Parameters::is_htl(x1)) -
                         Normalization::normalized_n(Parameters::N0 * Parameters::is_perovskite(x1)) //N0 value
                         ;

        //std::cout << " D = " << D << std::endl;
        double normalised_carriers=0.0;                 
         if(Parameters::is_perovskite(x1)){
          normalised_carriers= -_P[i];
         }                
         else if(Parameters::is_htl(x1)){
          normalised_carriers = -_p[i];

         }
         else if(Parameters::is_etl(x1)){
          normalised_carriers=_n[i];
         }
        _vec[i] = -(c0*sol[i-1] + c1*sol[i] + c2*sol[i+1]) + normalised_carriers - D;
      }

      // Use the factors to solve the linear system
      _dsol = _solver.solve(_vec);
      // apply correction
      sol += URF*_dsol;
      // correct boundary conditions
      sol[0] = si1;
      sol[Ni-1] = si2;
      // convergence criteria
      const double dsol_l2norm = _dsol.squaredNorm();
      dsol_l2norm_sum += dsol_l2norm;
#ifdef printInnerConvergence
      std::cout << "iter = " << iter << ", dsol_l2norm = " << dsol_l2norm << std::endl;
#endif
      if (dsol_l2norm < dsol_abs_tol)
        break;
      iteration++;
      //std::cout << "iteration = " << iteration << ", dsol_l2norm = " << dsol_l2norm << std::endl;
    } while (iter++ < iter_max);
  }

  return dsol_l2norm_sum;
}