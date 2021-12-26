#ifndef Charge_h
#define Charge_h

#include <string>

class Charge
{
public:
  enum ChargeSign
  {
    n = -1,
    p = 1,
    P = 2 //number
  };

  static std::string ChargeName(ChargeSign sign)
  {
    return (sign==n) ? "electron" : ((sign==p) ? "hole" : "ion");
  }

  static std::string ChargeSymbol(ChargeSign sign)
  {
    return (sign==n) ? "n" : ((sign==p) ? "p" : "P");
  }
};

#endif