#pragma once

#include <time.h>

/**
 * @class Timing 
 * A class for time measurements of the eigen decomposition.
 * Timing measures overall time and times between calls of tic and toc. For
 * small time spans (up to 1000 seconds) CPU time via clock() is used. For large
 * time spans the fall-back to elapsed time from time() is used.
 * timings_update() must be called often enough to prevent the fallback.
 */
class Timing
{
  clock_t lastclock;
  time_t lasttime;
  clock_t ticclock;
  time_t tictime;
  short istic;
  short isstarted;

  double lastdiff;
  double tictoczwischensumme;

public:
  double totaltime; //! zeroed by re-calling timings_start
  double totaltotaltime;
  double tictoctime;
  double lasttictoctime;

  Timing();
  void start();
  /**
   * @return time between last call of timings_*() and now
   */
  double update();
  void tic();
  double toc();
};
