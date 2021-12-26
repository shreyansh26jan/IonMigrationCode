#ifndef Parameters_h
#define Parameters_h

#include <iostream>
#include <cmath>

class Parameters
{
public:
  static constexpr double q = 1.602176634e-19; // elementary charge, Coulomb = Ampere Second
  static constexpr double k_B = 1.38064852e-23; // Boltzmann constant, m^2 kg s^{-2} K^{-1}
  static constexpr double eps0 = 8.85418782e-12; // vacuum permittivity, m^{-3} kg^{-1} s^4 A^2, F/m
  static constexpr double T = 300; // room temperature, K
  static constexpr double V_T = k_B * T / q; // thermal voltage, Volt

  // perovskite absorber layer
  static constexpr double L_perovskite = 400e-9; //m
  static constexpr double eps_rel = 24.1; // perovskite relative permittivity
  static constexpr double E_gap = 1.7; // perovskite band gap, eV
  static constexpr double N_C = 8.1e24; // density of states, m^{-3}
  static constexpr double N_V = 5.8e24; // density of states, m^{-3}
  static constexpr double E_C = -5.4; //-5.43 perovskite conduction band gap, eV
  static constexpr double E_V = -3.7; //-3.88 perovskite valence band gap, eV
  static constexpr double C_n = 5e-17; // electron capture coefficient, m^3 s^{-1}
  static constexpr double C_p = 4e-17; // hole capture coefficient, m^3 s^{-1}
  static constexpr double N0 = 1.6e21; // Density of ion vacancies, m^{-3} // Richardson paper

// HP
  // static constexpr double mu_n_perovskite = 3.23e-6; // electron mobility, m^2 V^{-1} s^{-1} // 3.23e-6 HP, 4.07e-6 MW
  // static constexpr double mu_p_perovskite = 3e-5; // hole mobility, m^2 V^{-1} s^{-1}
  // static constexpr double N_trap_perovskite = 1e21; // perovskite trap density, m^{-3} // 2.58e22 HP, 1.42e22 MW
  // static constexpr double N_trap_etl = 3.65e24; // perovskite/etl interface trap density, m^{-3} // 3.65e24 HP, 1.4e24 MW
  // static constexpr double P_trap_perovskite = 1e21; // perovskite trap density, m^{-3} fit
  // static constexpr double mu_n_etl = 1e-6; // m^2 Vs^{-1}, fit by koster
  // static constexpr double P_trap_htl = 1.0e22; // perovskite/htl interface trap density, m^{-3}

// MW
  static constexpr double DiffusionCoefficient_electron=1.7e-4;
  static constexpr double DiffusionCoefficient_hole=1.7e-4;

  static constexpr double mu_n_perovskite = DiffusionCoefficient_electron/V_T; //actual-1e-6 6e-6 1e-5 electron mobility, m^2 V^{-1} s^{-1} // 3.23e-6 HP, 4.07e-6 MW
  static constexpr double N_trap_perovskite = 9e20; //actual 9e20 1.42perovskite trap density, m^{-3} // 2.58e22 HP, 1.42e22 MW
  static constexpr double mu_p_perovskite = DiffusionCoefficient_hole/V_T; //actual 1e-5 4.7 hole mobility, m^2 V^{-1} s^{-1}
  static constexpr double N_trap_etl = 1.4e23; //1.4e23 perovskite/etl interface trap density, m^{-3} // 3.65e24 HP, 1.4e24 MW
  static constexpr double P_trap_perovskite = 1e21; //actual 1e21 perovskite trap density, m^{-3} fit
  static constexpr double mu_n_etl = 1e-6; //important 1e-6 m^2 Vs^{-1}, fit by koster
  static constexpr double P_trap_htl = 1.0e22; //e22 perovskite/htl interface trap density, m^{-3}
  
  static constexpr double DiffusionCoefficient_ionvacancy = 1e-17; //m^2 s^{-1}
  static constexpr double mu_P_perovskite = DiffusionCoefficient_ionvacancy/V_T;//value needed
  static constexpr double G_max =  5.38661687656549e+27; 

  static constexpr double tou_ion=1e-5; //characteristic timescale for ion motion




; //4.473,4.649,4.74original(1.65) maximum charge generation rate, m^{-3}s^{-1}
  static constexpr double k_b = 1e-17; // bimolecular recombination constant, m^3 s^{-1}

  // ETL (electron transporting layer) c60 TiO2
  static constexpr double L_etl = 35e-9; // m
  static constexpr double eps_n_rel_etl = 10; // TiO2 relative permittivity=20
  static constexpr double L_interface_etl = 5e-9; // m
  static constexpr double E_lumo_etl = -4.0; // PCBM LUMO level, eV
  static constexpr double mu_p_etl = 1e-5/V_T; // m^2 Vs^{-1}
  static constexpr double N_D_etl = 5e20; // Ionized doping in PCBM, m^{-3} // fit by Koster
  static constexpr double phi_n_etl = 0.05; // injection barrier height for electrons, eV
  static constexpr double mu_P_etl = 1e-15;


  // HTL (hole transporting layer) pedot pss (spiro-OMeTAD)
  static constexpr double L_htl = 100e-9; // m
  static constexpr double eps_p_rel_htl = 3; // polyTPD relative permittivity
  static constexpr double L_interface_htl = 5e-9; // m
  static constexpr double E_homo_htl = -5.1; // polyTPD HOMO level, eV
  static constexpr double mu_P_htl = 1e-15;
  static constexpr double mu_p_htl = 4e-7; //7e-5 m^2 Vs^{-1}, fit by koster
  static constexpr double mu_n_htl = 1e-9; // m^2 Vs^{-1}
  static constexpr double N_A_htl = 5e20; // Ionized doping in polyTPD, m^{-3} // fit by Koster
  static constexpr double phi_p_htl = 0.05; // injection barrier height for holes, eV

  static constexpr double n1() { return N_C * exp(-0.5*E_gap / V_T); }
  static constexpr double p1() { return N_C * exp(-0.5*E_gap / V_T); }

  // normalization parameters
  static constexpr double N_max = 1e20;//std::max(N_D_etl,std::max(N_A_htl, N0) ); // maximum dopant density, m^{-3}
  static constexpr double l_c = L_perovskite + L_etl + L_htl; // characteristic length, m
  static constexpr double mu_max = std::max(mu_n_perovskite, std::max(mu_p_perovskite, std::max(mu_n_etl, std::max(mu_p_htl, mu_P_perovskite))));
  static constexpr double D_max = mu_max * V_T;
  static constexpr double eps_max = eps0 * std::max(eps_p_rel_htl, std::max(eps_rel, eps_n_rel_etl));
  static constexpr double lambda2 = eps_max * V_T / (l_c*l_c * N_max * q);

  static constexpr double n_int() { return N_C * exp(-0.5 * E_gap / V_T); }
  static constexpr double V_bi = 1.1; // 1.127 built-in voltage at short circuit, V

  // parameter fields
  // ASSUMPTION: etl - perovskite layer - htl
  // NOTE: denormalized value is supplied as parameter x
  static constexpr double is_perovskite(double x) { return (x >= L_etl && x <= L_etl+L_perovskite) ? 1 : 0; }
  static constexpr double is_etl(double x) { return (x < L_etl) ? 1 : 0; }
  static constexpr double is_htl(double x) { return (x > L_htl+L_perovskite) ? 1 : 0; }

  static constexpr double eps(double x) { return eps0 * (is_htl(x) ? eps_p_rel_htl : (is_etl(x) ? eps_n_rel_etl : eps_rel)); }
  static constexpr double mu_n(double x) { return is_etl(x) ? mu_n_etl : (is_htl(x) ? mu_n_htl : mu_n_perovskite); }
  static constexpr double mu_p(double x) { return is_htl(x) ? mu_p_htl : (is_etl(x) ? mu_p_etl : mu_p_perovskite); }
  static constexpr double mu_P(double x) { return is_htl(x) ? mu_P_htl : (is_etl(x) ? mu_P_etl : mu_P_perovskite); }

  static constexpr double is_interface_etl(double x)
  { return (x >= L_etl - 0.5*L_interface_etl && x <= L_etl + 0.5*L_interface_etl) ? 1 : 0; }

  static constexpr double is_interface_htl(double x)
  { return (x >= L_etl+L_perovskite - 0.5*L_interface_htl && x <= L_etl+L_perovskite + 0.5*L_interface_htl) ? 1 : 0; }

  static constexpr double p_trap_htl(double x) { return P_trap_htl * is_interface_htl(x); }
  static constexpr double n_trap_etl(double x) { return N_trap_etl * is_interface_etl(x); }

  static constexpr double p_trap_perovskite(double x) { return P_trap_perovskite * is_perovskite(x); }
  static constexpr double n_trap_perovskite(double x) { return N_trap_perovskite * is_perovskite(x); }
};

#endif
