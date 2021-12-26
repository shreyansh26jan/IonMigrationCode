#include "Utils.h"
#include "IterativeSolver.h"
#include "Parameters.h"
#include "PostProcessing.h"

#include <iostream>
#include <iomanip>
#include<math.h>
#include "GenerationRecombination.h"

int main(int argc, char **argv)
{
  // grid
  const int Ni = std::floor(Parameters::l_c * 1e9);
  const double L = 1;
  std::vector<double> xi = Utils::LinearSpacedArray(0, L, Ni);
  /*std::cout << Ni<<"The vector elements are : ";
  for(int i=0; i < xi.size(); i++)
  std::cout << xi[i] << std::endl;*/
  IterativeSolver solver(xi);
  //IterativeSolver unsteadysolver(xi);
  std::vector<double> V;
  std::vector<double> J;

  std::vector<double> _si0;     //initially these contain 0's for first iteration
  std::vector<double> _n0;
  std::vector<double> _p0;
  std::vector<double> _P0;
  
  std::vector<double> _si1;     //initially these contain 0's for first iteration of reverse scan
  std::vector<double> _n1;
  std::vector<double> _p1;
  std::vector<double> _P1;

  std::vector<double> justSi;

  //initialising the vectors of t=0
  for(int i=1;i<=Ni;i++){
    _si0.push_back(0);
    _n0.push_back(0);
    _p0.push_back(0);
    _P0.push_back(0);
    _si1.push_back(0);
    _n1.push_back(0);
    _p1.push_back(0);
    _P1.push_back(0);

  }

  std::cout << std::setprecision(3) << std::scientific;
  std::cout << "V_a,        J0,        J1,         dsol_l2norm_sum" << std::endl;
  //unsteady state
  double scan_rate=0.08; // (V/s) 1) mv try different scan rate - potential vs 
  scan_rate=scan_rate*0.1;
  for (double V_a = 0; V_a < 0.97; V_a +=scan_rate)//0.94
  {
    /*if(V_a>0.95){
    std::cout << "###" << std::endl;
    std::cout<< V_a << ","<<_si0[0]<< ","<<_n0[0]<< ","<<_p0[0]<< ","<<_P0[0]<<std::endl;
    std::cout << "###" << std::endl;
    std::cout<< V_a << ","<<_si0[1]<< ","<<_n0[1]<< ","<<_p0[1]<< ","<<_P0[1]<<std::endl;
    std::cout << "##end#" << std::endl;
    }*/
    Solution sol = solver.solve(V_a,_si0,_n0,_p0,_P0);     //add the parameter _si0,_n0,_p0,_P0
    std::cout << V_a << ", " << sol.J_0 << ", " << sol.J_1 << ", " << sol.dsol_l2norm_sum << std::endl;
    
    //std::cout << "###" <<V_a<<" , "<<(int)(V_a/0.08)<<"##"<< std::endl;

    /*if( std::fmod( V_a,0.08)<1.0e-10){
      std::cout << "###" << std::endl;
      for(int i = 1; i < Ni-1; i++){
        justSi.push_back(_si0[i]);
        std::cout << _si0[i] << std::endl;

      }
      std::cout << "###" << std::endl;
    }*/
    V.push_back(V_a);
    
    J.push_back(sol.J_0/10); // A/m^2 to mA/cm^2
    
    /*_si0.clear();
    _n0.clear();
    _p0.clear();
    _P0.clear();


    _si0.push_back();               //we will get a vector here of current time, we will pass this along with V_a for next time step
    _n0.push_back();
    _p0.push_back();                 // enter paramter here
    _P0.push_back();

    //break;*/
  }

  for (double V_a = 0.96; V_a >= 0.0; V_a -=scan_rate)
  {
    /*if(V_a>0.95){
    std::cout << "###" << std::endl;
    std::cout<< V_a << ","<<_si0[0]<< ","<<_n0[0]<< ","<<_p0[0]<< ","<<_P0[0]<<std::endl;
    std::cout << "###" << std::endl;
    std::cout<< V_a << ","<<_si0[1]<< ","<<_n0[1]<< ","<<_p0[1]<< ","<<_P0[1]<<std::endl;
    std::cout << "##end#" << std::endl;
    }*/
    Solution sol = solver.solve(V_a,_si1,_n1,_p1,_P1);     //add the parameter _si0,_n0,_p0,_P0
    std::cout << V_a << ", " << sol.J_0 << ", " << sol.J_1 << ", " << sol.dsol_l2norm_sum << std::endl;
    
    V.push_back(V_a);
    
    J.push_back(sol.J_0/10); // A/m^2 to mA/cm^2

  }  
  // unsteadystate
  
  
  PostProcessing::plotJV(J, V,justSi);

  return 0;
}
