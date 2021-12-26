#include "CurrentContinuityEquation.h"
#include "BoundaryConditions.h"
#include "DebugMessages.h"
#include "GenerationRecombination.h"
#include "PostProcessing.h"

#include <cmath>
#include <iomanip>
#include <iostream>

CurrentContinuityEquation::
CurrentContinuityEquation(std::vector<double> const &xi, Eigen::VectorXd const &si,
                          Eigen::VectorXd const &n, Eigen::VectorXd const &p, Eigen::VectorXd const &P)
: _xi(xi)
, _si(si)
, _n(n)
, _p(p)
,_P(P)
{
  const int Ni = _xi.size();
  _tripletList.resize((Ni-2)*3+2);
  _mat.resize(Ni, Ni);
  _vec.resize(Ni);
  _res.resize(Ni);
  _dsol.resize(Ni);
}


void
CurrentContinuityEquation::
initialize(Charge::ChargeSign chargeSign, Eigen::VectorXd &sol, std::vector<double> &_si0, std::vector<double> &_n0, std::vector<double> &_p0, std::vector<double> &_P0)
{
  const int Ni = _xi.size();
  //const double b1 = (chargeSign == Charge::n) ? BoundaryConditions::n1() : BoundaryConditions::p1();
  //const double b2 = (chargeSign == Charge::n) ? BoundaryConditions::n2() : BoundaryConditions::p2();
  const double b1 = (chargeSign == Charge::n) ? BoundaryConditions::n1() : ((chargeSign == Charge::p) ? BoundaryConditions::p1():BoundaryConditions::P1());
  const double b2 = (chargeSign == Charge::n) ? BoundaryConditions::n2() : ((chargeSign == Charge::p) ? BoundaryConditions::p2():BoundaryConditions::P2());



  const std::vector<double> val= (chargeSign == Charge::n) ? _n0 : ((chargeSign == Charge::p) ? _p0 :_P0);
  bool check=false; //false means 0s are present and true means non 0s are present
  for (int i=0;i<Ni;i++)
  {
    if(val[i]!=0)
    {
      check = true;
      break;
    }
  } 

  if(check==false){
  // linear
   if(chargeSign == Charge::P){
      sol.setLinSpaced(Ni, b1+0.5, b2+0.5);
   }
   else{
    sol.setLinSpaced(Ni, b1, b2);
   }
 }
 else
 {
  for(int i=0;i<Ni;i++){
    sol[i]=val[i];

  }
  sol[0]=b1;
  sol[Ni-1]=b2;
 }
   
  //check boundary condition of anion vacancies and arrage sol for loop accordingly 
  // exponential
  /*const double c2b = std::log(b2/b1)/(_xi[Ni-1]-_xi[0]);
  for (int i=0; i<Ni; i++)
  {
    sol[i] = b1 * std::exp(c2b * _xi[i]);
  }*/
}


namespace
{
  // Bernoulli generating function
  // B(z) = z / (exp(z) - 1), B(0) = 1
  double B(double z)
  {
    constexpr double tol = 1e-16;

    if (fabs(z) < tol)
    {
      return 1;
    }

    return z / expm1(z);
  }
}


double
CurrentContinuityEquation::
solve(Charge::ChargeSign chargeSign, Eigen::VectorXd &sol, std::vector<double> &_si0, std::vector<double> &_n0, std::vector<double> &_p0, std::vector<double> &_P0)
{
  const double b1 = (chargeSign == Charge::n) ? BoundaryConditions::n1() : ((chargeSign == Charge::p) ? BoundaryConditions::p1():BoundaryConditions::P1());
  const double b2 = (chargeSign == Charge::n) ? BoundaryConditions::n2() : ((chargeSign == Charge::p) ? BoundaryConditions::p2():BoundaryConditions::P2());
  const double lambda_t=Parameters::l_c * Parameters::l_c / Parameters::tou_ion / Parameters::mu_max / Parameters::V_T; //value needs to be added
  double delta_t;
  const std::vector<double> val= (chargeSign == Charge::n) ? _n0 : ((chargeSign == Charge::p) ? _p0 :_P0);
  bool check=false; //false means 0s are present and true means non 0s are present
  for (int i=0;i<val.size();i++)
  {
    if(val[i]!=0)
    {
      check = true;
      break;
    }
  } 

  if(check==false){
  // linear
   delta_t=1e15;
 }
 else
 {
  delta_t=0.1;
  }

  const int Ni = _xi.size();
  //const std::vector<double> val= (chargeSign == Charge::n) ? _n0 : ((chargeSign == Charge::p) ? _p0 :_P0);
  // Solution using SG method
  // assembly
  typedef Eigen::Triplet<double> T;

  // boundary conditions
  int tripletListIndex = 0;
  _tripletList[tripletListIndex++] = T(0, 0, 1);
  _tripletList[tripletListIndex++] = T(Ni-1, Ni-1, 1);
  /*std::cout << "The vector elements are : ";
  std::cout << T(0,0,1) <<","<< std::endl;
  */
  _vec[0] = b1;
  _vec[Ni-1] = b2;
  
  constexpr double RFactorNormalization = Parameters::l_c * Parameters::l_c * Parameters::N_max / Parameters::D_max;
  int bernouliifunction; // used it in B in co,c1,c2 calculation
  // fill A,b in Ax=b
  for(int i = 1; i < Ni-1; i++)
  {
    const double x0 = Normalization::denormalized_x(_xi[i-1]);
    const double x1 = Normalization::denormalized_x(_xi[i]);
    const double x2 = Normalization::denormalized_x(_xi[i+1]);

    const double mu0 = (chargeSign == Charge::n) ? Parameters::mu_n(x0) : ((chargeSign == Charge::p) ? Parameters::mu_p(x0) : Parameters::mu_P(x0));
    const double mu1 = (chargeSign == Charge::n) ? Parameters::mu_n(x1) : ((chargeSign == Charge::p) ? Parameters::mu_p(x1) : Parameters::mu_P(x1));
    const double mu2 = (chargeSign == Charge::n) ? Parameters::mu_n(x2) : ((chargeSign == Charge::p) ? Parameters::mu_p(x2) : Parameters::mu_P(x2));

    const double dxi_minus = 0.5 * (mu0 + mu1) / (_xi[i] - _xi[i-1]) / Parameters::mu_max;
    const double dxi_plus  = 0.5 * (mu1 + mu2) / (_xi[i+1] - _xi[i]) / Parameters::mu_max;
    const double dxi       = 0.5 * (_xi[i+1] - _xi[i-1]);

    
    // bimolecular
    const double k_b = (chargeSign == Charge::P) ? 0 : Parameters::k_b;

    // srh
    const double k_srh0 = (chargeSign == Charge::P) ? 0 : GenerationRecombination::k_srh(_n[i-1], _p[i-1], x1); //check srh and other charges for positive ion vacancy
    const double k_srh1 = (chargeSign == Charge::P) ? 0 : GenerationRecombination::k_srh(_n[i],   _p[i],   x1);
    const double k_srh2 = (chargeSign == Charge::P) ? 0 : GenerationRecombination::k_srh(_n[i+1], _p[i+1], x1);

    const double other_charge0 = (chargeSign == Charge::n) ? _p[i-1] : _n[i-1];
    const double other_charge1 = (chargeSign == Charge::n) ? _p[i]   : _n[i];
    const double other_charge2 = (chargeSign == Charge::n) ? _p[i+1] : _n[i+1];

    const double RFactor0 = 0.25 * dxi * (k_b + k_srh0) * other_charge0 * RFactorNormalization;
    const double RFactor1 = 0.5  * dxi * (k_b + k_srh1) * other_charge1 * RFactorNormalization;
    const double RFactor2 = 0.25 * dxi * (k_b + k_srh2) * other_charge2 * RFactorNormalization;
    if(chargeSign == Charge::P){
      bernouliifunction = 1 ;  // changed the value back to 1 as in charge headr file I set it as 2 to distinguish between P and p

    }
    else{
      bernouliifunction = chargeSign;
    }
    // i - 1
    const double c0 = dxi_minus * B(bernouliifunction * (_si[i] - _si[i-1])) - RFactor0 - lambda_t*dxi/(4*delta_t);
    _tripletList[tripletListIndex++] = T( i, i-1, c0 );

    // i
    const double c1 = (-dxi_plus  * B(bernouliifunction * (_si[i+1] - _si[i]))
                       -dxi_minus * B(bernouliifunction * (_si[i-1] - _si[i]))) - RFactor1 - lambda_t*dxi/(2*delta_t);
    _tripletList[tripletListIndex++] = T( i, i, c1 );

    // i + 1
    const double c2 = dxi_plus * B(bernouliifunction*(_si[i] - _si[i+1])) - RFactor2 - lambda_t*dxi/(4*delta_t);
    _tripletList[tripletListIndex++] = T( i, i+1, c2 );

    const double n_int_squared = Parameters::n_int()*Parameters::n_int();
    const double integral_k_n_int_squared_denormalized = dxi * (0.25*(k_b + k_srh0) + 0.5*(k_b + k_srh1) + 0.25*(k_b + k_srh2)) * n_int_squared;
    const double integral_k_n_int_squared = GenerationRecombination::G_normalized(integral_k_n_int_squared_denormalized);

    const double G0 = (chargeSign == Charge::P) ? 0 : GenerationRecombination::G_normalized(x0);
    const double G1 = (chargeSign == Charge::P) ? 0 : GenerationRecombination::G_normalized(x1);
    const double G2 = (chargeSign == Charge::P) ? 0 : GenerationRecombination::G_normalized(x2);
    const double integral_G = dxi * (0.25*G0 + 0.5*G1 + 0.25*G2);

    _vec[i] = - integral_G - integral_k_n_int_squared-(lambda_t*dxi*(val[i]+(val[i-1]+val[i+1])/2)/(2*delta_t));
    //if(chargeSign == Charge::p && Parameters::is_perovskite(x0) )
      //std::cout << _vec[i]<<" total "<<- integral_G - integral_k_n_int_squared<< "left" << -(lambda_t*dxi*(val[i]+(val[i-1]+val[i+1])/2)/(2*delta_t))<<"right"<< std::endl;
    if(i==0  )
    { 
      std::cout << chargeSign <<  std::endl; 
      std::cout << "Rfactor 0 1 2 " << RFactor0 << "  " <<RFactor1 <<"  " <<RFactor2 << std::endl;
      std::cout << "dxi - 0 + " << dxi_minus <<"  " << dxi <<"  " << dxi_plus << std::endl;
      std::cout << "c0 c1 c2 " << c0 <<"  " << c1 <<"  " << c2 << std::endl;
      std::cout << "vec " << _vec[i] << std::endl;
    }
  }
  
#ifdef printInnerConvergence
  std::cout << Charge::ChargeName(chargeSign)
  << " [" << Charge::ChargeSymbol(chargeSign) << "]:" << std::endl;
#endif

  // convergence criteria
  constexpr int iter_max = 0;
  constexpr double dsol_abs_tol = 1e-20;
  constexpr double URF = 1;

  _mat.setFromTriplets(_tripletList.begin(), _tripletList.end());
  /*std::cout << "The vector elements are : ";
  std::cout << _mat <<","<< std::endl;
  */
  if (!_solver)
  {
    _solver.reset(new Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >());
    // Compute the ordering permutation vector from the structural pattern of A
    _solver->analyzePattern(_mat);
  }
  
  // Compute the numerical factorization
  _solver->factorize(_mat);
  

  double dsol_l2norm_sum = 0;
  {
    // solve \delta x in A\delta x = b - Ax_n
    int iter = 0;
    do
    {
      // b-Ax_n in A\delta x = b-Ax_n
      _res = _vec - _mat*sol;
      //std::cout << "inside " << _res << std::endl;
      //Use the factors to solve the linear system
      _dsol = _solver->solve(_res);
      //std::cout << "inside1" <<  std::endl;
      // apply correction
      sol += URF*_dsol;
      // correct boundary conditions
      sol[0] = b1;
      sol[Ni-1] = b2;
      // convergence criteria
      const double dsol_l2norm = _dsol.squaredNorm();
      dsol_l2norm_sum += dsol_l2norm;
#ifdef printInnerConvergence
      std::cout << "iter = " << iter << ", dsol_l2norm = " << dsol_l2norm << std::endl;
#endif
      if (dsol_l2norm < dsol_abs_tol)
        break;
    } while (iter++ < iter_max);
  }
  //std::cout << "after do" << std::endl;
  if (_solver->info() != Eigen::Success)
  {
    std::cout << "Failed to solve !" << std::endl;
  }
  //std::cout << "end" << std::endl;
  return dsol_l2norm_sum;
}

