#include "utils.h"
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>

std::vector<std::string> split(const std::string line, char separator)
{
  std::vector<std::string> elements;
  size_t last = 0;
  size_t i;
  for(i = 0; i < line.length(); ++i)
  {
    if(line[i] == separator)
    {
      if(i - last > 0)
        elements.push_back(line.substr(last, i - last));
      last = i + 1;
    }
  }
  if(i - last > 0)
    elements.push_back(line.substr(last, i - last));
  return elements;
}

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
