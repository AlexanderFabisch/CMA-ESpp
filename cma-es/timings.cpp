#include "timings.h"
#include "utils.h"
#include <cassert>

Timing::Timing()
{
  totaltotaltime = 0;
  start();
}

void Timing::start()
{
  totaltime = 0;
  tictoctime = 0;
  lasttictoctime = 0;
  istic = 0;
  lastclock = clock();
  lasttime = time(NULL);
  lastdiff = 0;
  tictoczwischensumme = 0;
  isstarted = 1;
}

double Timing::update()
{
  double diffc, difft;
  clock_t lc = lastclock; // measure CPU in 1e-6s
  time_t lt = lasttime; // measure time in s

  assert(isstarted == 1 && "timings_started() must be called before using timings... functions");

  lastclock = clock(); // measures at most 2147 seconds, where 1s = 1e6 CLOCKS_PER_SEC
  lasttime = time(NULL);
  diffc = (double) (lastclock - lc) / CLOCKS_PER_SEC; // is presumably in [-21??, 21??]
  difft = difftime(lasttime, lt); // is presumably an integer
  lastdiff = difft; // on the "save" side
  // use diffc clock measurement if appropriate
  if(diffc > 0 && difft < 1000)
    lastdiff = diffc;
  assert(lastdiff >= 0 && "BUG in time measurement");
  totaltime += lastdiff;
  totaltotaltime += lastdiff;
  if(istic)
  {
    tictoczwischensumme += lastdiff;
    tictoctime += lastdiff;
  }
  return lastdiff;
}

void Timing::tic()
{
  assert(!istic && "Timingic called twice without toc");
  update();
  istic = 1;
}

double Timing::toc()
{
  assert(istic && "Timingoc called without tic");
  update();
  lasttictoctime = tictoczwischensumme;
  tictoczwischensumme = 0;
  istic = 0;
  return lasttictoctime;
}
