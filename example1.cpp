/**
 * @file example1.cpp
 * Very short example source code. The purpose of the example codes is to be
 * edited/extended.
 */

#include <stdlib.h>
#include <iostream>
#include "cmaes.h"

/** the objective (fitness) function to be minized */
double fitfun(double const *x, int N)
{ // function "cigtab"
  double sum = 1e4*x[0]*x[0] + 1e-4*x[1]*x[1];
  for(int i = 2; i < N; ++i)  
    sum += x[i]*x[i]; 
  return sum;
}

/** the optimization loop */
int main(int argn, char **args) {
  CMAES<double> evo;
  double *arFunvals, *const*pop, *xfinal;

  // Initialize everything into the struct evo.
  const int dim = 22;
  double xstart[dim];
  for(int i=0; i<dim; i++) xstart[i] = 0.5;
  double stddev[dim];
  for(int i=0; i<dim; i++) stddev[i] = 0.3;
  Parameters<double> parameters;
  // TODO adjust parameters here
  parameters.init(dim, xstart, stddev);
  arFunvals = evo.init(parameters);

  std::cout << evo.sayHello() << std::endl;

  // Iterate until stop criterion holds
  while(!evo.testForTermination())
  {
    // generate lambda new search points, sample population
    pop = evo.samplePopulation(); /* do not change content of pop */

    /* Here you may resample each solution point pop[i] until it
       becomes feasible, e.g. for box constraints (variable
       boundaries). function is_feasible(...) needs to be
       user-defined.  
       Assumptions: the feasible domain is convex, the optimum is
       not on (or very close to) the domain boundary, initialX is
       feasible and initialStandardDeviations are sufficiently small
       to prevent quasi-infinite looping.
    */
    /* for (i = 0; i < cmaes_Get(&evo, "popsize"); ++i) 
         while (!is_feasible(pop[i])) 
           cmaes_ReSampleSingle(&evo, i); 
    */

    // evaluate the new search points using fitfun from above
    for (int i = 0; i < evo.get(CMAES<double>::Lambda); ++i)
      arFunvals[i] = fitfun(pop[i], (int) evo.get(CMAES<double>::Dimension));

    // update the search distribution used for sampleDistribution()
    evo.updateDistribution(arFunvals);
  }
  std::cout << "Stop:" << std::endl << evo.getStopMessage() << std::endl;
  evo.writeToFile(CMAES<double>::WKResume, "resumeevo1.dat"); // write resumable

  // get best estimator for the optimum, xmean
  xfinal = evo.getNew(CMAES<double>::XMean); // "xbestever" might be used as well

  // do something with final solution and finally release memory
  delete[] xfinal;

  return 0;
}

