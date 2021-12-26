#ifndef Normalization_h
#define Normalization_h

#include "Parameters.h"

class Normalization
{
public:
  static constexpr double normalized_n(double n) { return n / Parameters::N_max; }
  static constexpr double normalized_p(double p) { return p / Parameters::N_max; }
  static constexpr double normalized_P(double P) { return P / Parameters::N_max; }

  static constexpr double normalized_si(double si) { return si / Parameters::V_T; }
  static constexpr double normalized_J(double J) { return J * Parameters::l_c / (Parameters::D_max * Parameters::q * Parameters::N_max); }
  static constexpr double normalized_G(double G) { return G * Parameters::l_c * Parameters::l_c / (Parameters::D_max * Parameters::N_max); }

  static constexpr double denormalized_n(double n_normalized) { return n_normalized * Parameters::N_max; }
  static constexpr double denormalized_p(double p_normalized) { return p_normalized * Parameters::N_max; }
  static constexpr double denormalized_P(double P_normalized) { return P_normalized * Parameters::N_max; }

  static constexpr double denormalized_si(double si_normalized) { return si_normalized * Parameters::V_T; }
  static constexpr double denormalized_x(double x_normalized) { return x_normalized * Parameters::l_c; }

  // NOTE: denormalized value is supplied as parameter x
  static constexpr double denormalized_J(double J_normalized) { return J_normalized / Parameters::l_c * (Parameters::D_max * Parameters::q * Parameters::N_max); }
  // NOTE: denormalized value is supplied as parameter x
  static constexpr double denormalized_G(double G_normalized) { return G_normalized / Parameters::l_c / Parameters::l_c * (Parameters::D_max * Parameters::N_max); }
};

#endif