#include "random.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>

Random::Random(long unsigned inseed)
{
  stored = false;
  if(inseed < 1)
    inseed = (long unsigned) abs(100*time(NULL) + clock());
  start(inseed);
}

void Random::start(long unsigned inseed)
{
  stored = false;
  startseed = inseed;
  if(inseed < 1) inseed = 1;
  aktseed = inseed;
  for(int i = 39; i >= 0; --i)
  {
    long tmp = aktseed / 127773;
    aktseed = 16807* (aktseed - tmp* 127773) - 2836* tmp;
    if(aktseed < 0) aktseed += 2147483647;
    if(i < 32) rgrand[i] = aktseed;
  }
  aktrand = rgrand[0];
}

double Random::gauss()
{
  if(stored)
  {
    stored = false;
    return hold;
  }
  stored = true;
  double x1, x2, rquad;
  do {
    x1 = 2.0*uniform() - 1.0;
    x2 = 2.0*uniform() - 1.0;
    rquad = x1*x1 + x2*x2;
  } while(rquad >= 1 || rquad <= 0);
  const register double fac = sqrt(-2.0*log(rquad)/rquad);
  hold = fac*x1;
  return fac*x2;
}

double Random::uniform()
{
  long tmp = aktseed / 127773;
  aktseed = 16807 * (aktseed - tmp * 127773) - 2836 * tmp;
  if(aktseed < 0)
    aktseed += 2147483647;
  tmp = aktrand / 67108865;
  aktrand = rgrand[tmp];
  rgrand[tmp] = aktseed;
  return (double) aktrand / 2.147483647e9;
}
