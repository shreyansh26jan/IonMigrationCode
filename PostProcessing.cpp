#include "GenerationRecombination.h"
#include "PostProcessing.h"

//#include "matplotlibcpp.h"
#include <assert.h>
#include <fstream>

PostProcessing::
PostProcessing(std::vector<double> const &xi, Eigen::VectorXd const &si,
               Eigen::VectorXd const &n, Eigen::VectorXd const &p, Eigen::VectorXd const &P)
: _xi(xi)
, _si(si)
, _n(n)
, _p(p)
, _P(P)
{
  const int Ni = _xi.size();
  assert(_si.size() == Ni && _n.size() == Ni && _p.size() == Ni && _P.size() == Ni);
}


double
PostProcessing::
J0_1st_order(double x0, double x1, double x2, double J1, double J2)
{
  return J1 + (J1-J2)*(x0-x1)/(x0-x2);
}


double
PostProcessing::
J0_2nd_order(double x0, double x1, double x2, double x3, double J1, double J2, double J3)
{
  double x12 = x1 + x2 - 2*x0;
  double x23 = x2 + x3 - 2*x0;
  double x01 = x0 + x1 - 2*x0;
  double x2_0 = x2 - x0;
  double x3_1 = x3 - x1;
  double x23_01 = x2 + x3 - x0 - x1;
  return J1*( x12*x23/x2_0/x23_01) +
          J2*(-x01*x23/x2_0/x3_1) +
          J3*( x01*x12/x23_01/x3_1);
}


std::vector<double>
PostProcessing::
Ji(Charge::ChargeSign chargeSign) const
{
  const int Ni = _xi.size();
  std::vector<double> Ji(Ni+1);
  Eigen::VectorXd const &sol = (chargeSign == Charge::n) ? _n : _p;
  for (int i = 1; i < Ni; i++)
  {
    const double x0 = Normalization::denormalized_x(_xi[i-1]);
    const double x1 = Normalization::denormalized_x(_xi[i]);

    const double mu0 = (chargeSign == Charge::n) ? Parameters::mu_n(x0) : Parameters::mu_p(x0);
    const double mu1 = (chargeSign == Charge::n) ? Parameters::mu_n(x1) : Parameters::mu_p(x1);

    const double mu_minus = 0.5 * (mu0 + mu1) / Parameters::mu_max;

    Ji[i] = - mu_minus * (_si[i] - _si[i-1]) / (_xi[i] - _xi[i-1]) *
            (sol[i] * exp(chargeSign*(_si[i] - _si[i-1])) - sol[i-1])  /
            expm1(chargeSign*(_si[i] - _si[i-1]));
    Ji[i] *= -1;
  }

  Ji[0] = PostProcessing::J0_2nd_order(_xi[0],_xi[1],_xi[2],_xi[3],Ji[1],Ji[2],Ji[3]);
  Ji[Ni] = PostProcessing::J0_2nd_order(_xi[Ni-1],_xi[Ni-2],_xi[Ni-3],_xi[Ni-4],Ji[Ni-1],Ji[Ni-2],Ji[Ni-3]);

  return Ji;
}


std::vector<double>
PostProcessing::
xi_J() const
{
  const int Ni = _xi.size();
  std::vector<double> xi_J(Ni+1);
  xi_J[0] = _xi[0];
  xi_J[Ni] = _xi[Ni-1];
  for (int i = 1; i < Ni; i++)
  {
    xi_J[i] = 0.5 * (_xi[i] + _xi[i-1]);
  }

  return xi_J;
}


namespace
{
  class CSVWriter
  {
  public:
    CSVWriter(std::string const filename)
    : _filename(filename)
    {
      _file.open(filename);
    }

    ~CSVWriter()
    {
      _file.close();
    }

    void addColumn(std::string const name, std::vector<double> const &values)
    {
      _header.push_back(name);
      _data.push_back(values);
    }

    void write()
    {
      int c=0;
      for (auto &val : _header)
      {
        if (c++ > 0)
          _file << ",";

        _file << val;
      }
      _file << "\n";

      // ASSUMPTION: all columns have same number of values
      if (_data.size() > 0)
      {
        for (int rowIdx=0; rowIdx<_data[0].size(); rowIdx++)
        {
          c=0;
          for (auto &col : _data)
          {
            if (c++ > 0)
              _file << ",";

            _file << col[rowIdx];
          }
          _file << "\n";
        }
      }
    }

  private:
    std::vector<std::vector<double> > _data;
    std::vector<std::string> _header;
    std::string const _filename;
    std::ofstream _file;
  };
}


void
PostProcessing::
plotCharges(double V_a) const
{
  std::cout << "Plotting..." << std::endl;
  const int Ni = _xi.size();
  std::vector<double> n(Ni);
  std::vector<double> p(Ni);
  std::vector<double> si(Ni);
  for (int i = 0; i < Ni; i++)
  {
    n[i] = Normalization::denormalized_n(_n[i]);
    p[i] = Normalization::denormalized_p(_p[i]);
    si[i] = Normalization::denormalized_si(_si[i]);
  }

  {
    //matplotlibcpp::named_semilogy("electron density [m^-3]",_xi, n);
    //matplotlibcpp::named_semilogy("hole density [m^-3]",_xi, p);
    //matplotlibcpp::ylim(1e18,1e22);
    //matplotlibcpp::legend();
    //matplotlibcpp::title("Applied voltage Va = " + std::to_string(V_a) + " volts");
    //matplotlibcpp::grid(true);
    //matplotlibcpp::show();
    //matplotlibcpp::save("chargeDensitiesVa="+ std::to_string(V_a) + ".png");
  }

  {
    //matplotlibcpp::named_plot("potential [V]",_xi, si);
    //matplotlibcpp::legend();
    //matplotlibcpp::title("Applied voltage Va = " + std::to_string(V_a) + " volts");
    //matplotlibcpp::grid(true);
    //matplotlibcpp::show();
    //matplotlibcpp::save("potentialVa="+ std::to_string(V_a) + ".png");
  }

  CSVWriter writer("npV_Va=" + std::to_string(V_a) + ".csv");
  writer.addColumn("normalized x", _xi);
  writer.addColumn("electron density [m^-3]", n);
  writer.addColumn("hole density [m^-3]", p);
  writer.addColumn("potential [V]", si);
  writer.write();
}


void
PostProcessing::
plotCurrent(double V_a) const
{
  const std::vector<double> x_J = xi_J();

  std::vector<double> Jn = Ji(Charge::n);
  for (int i = 0; i < x_J.size(); i++)
  {
    Jn[i] = Normalization::denormalized_J(Jn[i]);
  }

  std::vector<double> Jp = Ji(Charge::p);
  for (int i = 0; i < x_J.size(); i++)
  {
    Jp[i] = Normalization::denormalized_J(Jp[i]);
  }

  //matplotlibcpp::named_plot(Charge::ChargeName(Charge::n) + " current density [A/m^2]", x_J, Jn);
  //matplotlibcpp::named_plot(Charge::ChargeName(Charge::p) + " current density [A/m^2]", x_J, Jp);
  //matplotlibcpp::legend();
  //matplotlibcpp::title("Applied voltage Va = " + std::to_string(V_a) + " volts");
  //matplotlibcpp::grid(true);
  //matplotlibcpp::show();

  CSVWriter writer("J_Va=" + std::to_string(V_a) + ".csv");
  writer.addColumn("normalized x", _xi);
  writer.addColumn(Charge::ChargeName(Charge::n) + " current density [A/m^2]", Jn);
  writer.addColumn(Charge::ChargeName(Charge::p) + " current density [A/m^2]", Jp);
  writer.write();
}


void
PostProcessing::
plotR(double V_a) const
{
  const int Ni = _xi.size();
  std::vector<double> Rb(Ni);
  std::vector<double> Rsrh(Ni);
  std::vector<double> Rsrh_interface(Ni);
  std::vector<double> G(Ni);
  for (int i = 0; i < Ni; i++)
  {
    const double x0 = Normalization::denormalized_x(_xi[i]);
    Rb[i] = Normalization::denormalized_G(GenerationRecombination::Rb_normalized(_n[i], _p[i], x0));
    Rsrh[i] = Normalization::denormalized_G(GenerationRecombination::Rsrh_normalized(_n[i], _p[i], x0));
    Rsrh_interface[i] = Rsrh[i] * (Parameters::is_interface_etl(x0) + Parameters::is_interface_htl(x0));
    G[i] = Normalization::denormalized_G(GenerationRecombination::G_normalized(x0));
  }

  double integral_Rb = integralfx(Rb);
  double integral_Rsrh = integralfx(Rsrh);
  double integral_Rsrh_interface = integralfx(Rsrh_interface);
  double integral_G = integralfx(G);

  double total_recombination = 100 * (integral_Rb + integral_Rsrh) / integral_G;
  double bimolecular_recombination = 100 * integral_Rb / (integral_Rb + integral_Rsrh);
  double srh_recombination = 100 * integral_Rsrh / (integral_Rb + integral_Rsrh);
  double srh_interface_recombination = 100 * integral_Rsrh_interface / (integral_Rb + integral_Rsrh);
  double srh_activeLayer_recombination = 100 * (integral_Rsrh - integral_Rsrh_interface) / (integral_Rb + integral_Rsrh);

  std::cout << "V_a = " << V_a << ", total_R = " << total_recombination
  << "%\nbimolecular_R = " << bimolecular_recombination
  << "%, srh_R = " << srh_recombination
  << "%\nsrh_activeLayer_R = " << srh_activeLayer_recombination
  << "%, srh_interface_R = " << srh_interface_recombination << "%" << std::endl;

  {
    //matplotlibcpp::named_semilogy("R_bimolecular",_xi, Rb);
    //matplotlibcpp::named_semilogy("R_srh",_xi, Rsrh);
    //matplotlibcpp::named_semilogy("G",_xi, G);
    //matplotlibcpp::ylabel("Generation/Recombination rate [m^-3 s^-1]");
    //matplotlibcpp::legend();
    //matplotlibcpp::suptitle("Applied voltage Va = " + std::to_string(V_a) + " volts");

    const std::string total_R_str = "Total R = " + std::to_string(total_recombination);
    const std::string Rb_str = "Rb = " + std::to_string(bimolecular_recombination);
    const std::string Rsrh_str = "Rsrh = " + std::to_string(srh_recombination);

    //matplotlibcpp::title(total_R_str + "%, " + Rb_str + "%, " + Rsrh_str + "%");
    //matplotlibcpp::grid(true);
    //matplotlibcpp::show();
  }

  CSVWriter writer("R_Va=" + std::to_string(V_a) + ".csv");
  writer.addColumn("normalized x", _xi);
  writer.addColumn("R_bimolecular", Rb);
  writer.addColumn("R_srh", Rsrh);
  writer.addColumn("G", G);
  writer.write();
}

double
PostProcessing::
integralfx(std::vector<double> const &f) const
{
  std::vector<double> _xi_denormalized(_xi);

  for (int i=0; i<_xi_denormalized.size(); i++)
  {
    _xi_denormalized[i] = Normalization::denormalized_x(_xi_denormalized[i]);
  }

  if (f.size() != _xi.size())
    std::cout << "ERROR: PostProcessing::integralfx, f.size != _xi.size" << std::endl;

  double integral=0;
  for (int i=1; i<f.size(); i++)
  {
    integral += 0.5 * (f[i] + f[i-1]) * (_xi_denormalized[i] - _xi_denormalized[i-1]);
  }

  return integral;
}


void
PostProcessing::
plotJV(std::vector<double> const &J, std::vector<double> const &V,std::vector<double> const &justSi)
{
  if (V.size() < 2)
    return;

  double max_JV = 0;
  double max_JV_J = 0;
  double max_JV_V = 0;
  double Voc = 1e-3;

  if (J.size() != V.size())
    std::cout << "ERROR: PostProcessing::printJV, J.size != V.size" << std::endl;

  for (int i=0; i< J.size(); i++)
  {
    double const JV = J[i] * V[i];
    if (fabs(max_JV) < fabs(JV) && V[i] < 1.)
    {
      max_JV = JV;
      max_JV_J = J[i];
      max_JV_V = V[i];
    }

    if (i > 0)
    {
      if ((J[i] > 0) && (J[i-1] < 0) ||
          (J[i] < 0) && (J[i-1] > 0))
      {
        Voc = V[i-1] - J[i-1] * (V[i] - V[i-1]) / (J[i] - J[i-1]);
      }
    }
  }

  const std::string ff_str = "Jsc = " + std::to_string(J[0]) + ", Voc = " + std::to_string(Voc) + ", FF = " + std::to_string(max_JV / (J[0] * Voc));
  std::cout << ff_str << std::endl;
  std::cout << "J_max = " << max_JV_J << ", V_max = " << max_JV_V << std::endl;

  //matplotlibcpp::plot(V, J);
  //matplotlibcpp::ylabel("Current density [mA/cm^2]");
  //matplotlibcpp::xlabel("Voltage (V)");
  //matplotlibcpp::title(ff_str);
  //matplotlibcpp::grid(true);
  //matplotlibcpp::ylim(-25,10);
  //matplotlibcpp::show();

  CSVWriter writer("JV.csv");
  //writer.addColumn("Si (V)", justSi);
  
  writer.addColumn("Voltage (V)", V);
  writer.addColumn("Current density [mA/cm^2]", J);

  writer.write();
}

