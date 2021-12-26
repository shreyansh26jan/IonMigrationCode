#include "IterativeSolver.h"
#include "BoundaryConditions.h"
#include "DebugMessages.h"
#include "PostProcessing.h"

#include <iostream>

IterativeSolver::
IterativeSolver(std::vector<double> const &xi)
: _xi(xi)
{
  const int Ni = _xi.size();
  
  _n.resize(Ni);
  
  _p.resize(Ni);
  _P.resize(Ni);
  /*std::cout << "The vector elements are : ";
  std::cout << Ni  <<","<<_n <<","<<_p <<","<<_si <<","<< std::endl;
  */
  _si.resize(Ni);
  
  _ccEqn.reset(new CurrentContinuityEquation(_xi, _si, _n, _p,_P));
  //_ccEqn1.reset(new CurrentContinuityEquation(_xi, _si, _n, _p,_P));
  _pEqn.reset(new PoissonEquation(_xi, _n, _p,_P));
  

}


Solution
IterativeSolver::
solve(double V_a, std::vector<double> &_si0, std::vector<double> &_n0, std::vector<double> &_p0, std::vector<double> &_P0) //receive vectors _si0,_n0,_p0,_P0
{
  _ccEqn->initialize(Charge::P, _P,_si0, _n0, _p0, _P0); //pass kth time paramters as well
  _ccEqn->initialize(Charge::n, _n,_si0, _n0, _p0, _P0);
  _ccEqn->initialize(Charge::p, _p,_si0, _n0, _p0, _P0);
  


  _pEqn->initialize(V_a, _si);
  
  const int Ni = _xi.size();
  const int outer_iter_max = 2000;
  int outer_iter = 0;
  const double dsol_l2norm_sum_tol = 1e-15;
  double dsol_l2norm_sum = 0;

  do
  {
    dsol_l2norm_sum = 0;
    
    
    
       
    dsol_l2norm_sum += _ccEqn->solve(Charge::p, _p,_si0, _n0, _p0, _P0);
    dsol_l2norm_sum += _ccEqn->solve(Charge::P, _P,_si0, _n0, _p0, _P0);
    dsol_l2norm_sum += _ccEqn->solve(Charge::n, _n,_si0, _n0, _p0, _P0);
    dsol_l2norm_sum += _pEqn->solve(V_a, _si);
    //std::cout << "do while loop " << dsol_l2norm_sum << std::endl;
    //std::cout <<dsol_l2norm_sum <<" and "<<outer_iter  << std::endl;
    
#ifdef printOuterConvergence
    std::cout << "[" << outer_iter <<  "] dsol_l2norm_sum = " << dsol_l2norm_sum << std::endl;
#endif

    if (dsol_l2norm_sum < dsol_l2norm_sum_tol)
      break;

  } while (outer_iter++ < outer_iter_max);

  PostProcessing pp(_xi, _si, _n, _p, _P); //check for _P here

  const std::vector<double> Jni = pp.Ji(Charge::n);
  const std::vector<double> Jpi = pp.Ji(Charge::p);
  const double Jn0_2 = Normalization::denormalized_J(Jni[0]);
  const double Jp0_2 = Normalization::denormalized_J(Jpi[0]);
  const double Jn1_2 = Normalization::denormalized_J(Jni[Ni]);
  const double Jp1_2 = Normalization::denormalized_J(Jpi[Ni]);

#ifdef printOuterConvergence
  std::cout << "Jn0 = " << Jn0_2 << ", Jp0 = " << Jp0_2 << ", Jn1 = " << Jn1_2 << ", Jp1 = " << Jp1_2 << std::endl;
  std::cout << "J_0 = " << Jn0_2 + Jp0_2 << ", J_1 = " << Jn1_2 + Jp1_2 << ", J_0 - J_1 = " << Jn0_2 + Jp0_2 - Jn1_2 - Jp1_2 << std::endl;
#endif

#ifdef plot_n_p_V
  //pp.plotCharges(V_a);
  //pp.plotCurrent(V_a);
  //pp.plotR(V_a);
#endif
  for(int i=0;i<Ni;i++){
    _si0[i]=_si[i];             //storing the value of n,p,P,si for current time for future use 
    _n0[i]=_n[i];
    _p0[i]=_p[i];
    _P0[i]=_P[i];
  }
  return {Jn0_2 + Jp0_2, Jn1_2 + Jp1_2, dsol_l2norm_sum };
}
