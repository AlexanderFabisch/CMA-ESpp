#include "utils.h"
#include <string.h>
#include <math.h>
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>

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

void FATAL(const std::string& msg)
{
  time_t t = time(NULL);
  ERRORMESSAGE(msg);
  ERRORMESSAGE("*** Exiting CMAES ***");
  std::cout << std::endl << " -- " << asctime(localtime(&t)) << " " << msg
      << std::endl << " *** CMA-ES ABORTED, see errcmaes.err *** " << std::endl;
  exit(1);
}

void ERRORMESSAGE(const std::string& msg)
{
  time_t t = time(NULL);
  std::ofstream file("errcmaes.err", std::ios_base::app);
  if(!file.is_open())
  {
    std::cout << std::endl << "FATAL ERROR: " << msg << std::endl
        << "CMAES could not open file 'errcmaes.err'." << std::endl
        << " *** CMA-ES ABORTED *** ";
    exit(1);
  }
  else
  {
    file << std::endl << " -- " << asctime(localtime(&t)) << " " << msg;
    file.close();
  }
}
