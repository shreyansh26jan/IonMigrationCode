#ifndef BoundaryConditions_h
#define BoundaryConditions_h

#include "Normalization.h"
#include <cmath>

namespace BoundaryConditions
{
  constexpr double n_0() { return Parameters::N_C * exp(-Parameters::phi_n_etl/Parameters::V_T); }
  constexpr double n_L() { return Parameters::N_C * exp(-(Parameters::E_gap - Parameters::phi_p_htl)/Parameters::V_T); }

  constexpr double p_0() { return Parameters::N_V * exp(-(Parameters::E_gap - Parameters::phi_n_etl)/Parameters::V_T); }
  constexpr double p_L() { return Parameters::N_V * exp(-Parameters::phi_p_htl/Parameters::V_T); }

  constexpr double P_0() { return 0; }
  constexpr double P_L() { return 0; }


  constexpr double si_0(double V_a) { return 0.5 * (Parameters::V_bi - V_a); }
  constexpr double si_L(double V_a) { return -si_0(V_a); }

  // normalized values for convenience
  constexpr double n1() { return Normalization::normalized_n(BoundaryConditions::n_0()); }
  constexpr double n2() { return Normalization::normalized_n(BoundaryConditions::n_L()); }

  constexpr double p1() { return Normalization::normalized_p(BoundaryConditions::p_0()); }
  constexpr double p2() { return Normalization::normalized_p(BoundaryConditions::p_L()); }

  constexpr double P1() { return Normalization::normalized_P(BoundaryConditions::P_0()); }
  constexpr double P2() { return Normalization::normalized_P(BoundaryConditions::P_L()); }


  constexpr double si1(double V_a) { return Normalization::normalized_si(BoundaryConditions::si_0(V_a)); }
  constexpr double si2(double V_a) { return Normalization::normalized_si(BoundaryConditions::si_L(V_a)); }
}

#endif