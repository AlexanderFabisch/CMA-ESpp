#include "utils.h"
#include <cmath>

double myhypot(double a, double b)
{
  const register double fabsa = fabs(a), fabsb = fabs(b);
  if(fabsa > fabsb)
  {
    const register double r = b / a;
    return fabsa*sqrt(1.0+r*r);
  }
  else if(b != 0.0)
  {
    const register double r = a / b;
    return fabsb*sqrt(1.0+r*r);
  }
  else
    return 0.0;
}
