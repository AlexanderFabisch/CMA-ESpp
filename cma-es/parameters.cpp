#include "parameters.h"
#include "utils.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <limits>

Parameters::Parameters()
    : N(-1),
      xstart(0),
      typicalX(0),
      typicalXcase(false),
      rgInitialStds(0),
      rgDiffMinChange(0),
      stopMaxFunEvals(-1.0),
      facmaxeval(1.0),
      stopMaxIter(-1.0),
      stopTolFun(1e-12),
      stopTolFunHist(1e-13),
      stopTolX(0), // 1e-11*insigma would also be reasonable
      stopTolUpXFactor(1e3),
      lambda(-1),
      mu(-1),
      mucov(-1),
      mueff(-1),
      weights(0),
      damps(-1),
      cs(-1),
      ccumcov(-1),
      ccov(-1),
      facupdateCmode(1),
      weightMode(UNINITIALIZED_WEIGHTS),
      resumefile("")
{
  stStopFitness.flg = false;
  stStopFitness.val = -std::numeric_limits<double>::max();
  updateCmode.modulo = -1;
  updateCmode.maxtime = -1;
}

Parameters::Parameters(const Parameters& parameters)
{
  assign(parameters);
}

Parameters& Parameters::operator=(const Parameters& parameters)
{
  assign(parameters);
  return *this;
}

void Parameters::init(int dimension, const double* inxstart, const double* inrgsigma)
{
  if(dimension <= 0 && N <= 0)
    FATAL("Problem dimension N undefined.");
  else if(dimension > 0)
    N = dimension;

  if(weightMode == UNINITIALIZED_WEIGHTS)
    weightMode = LOG_WEIGHTS;

  diagonalCov = 0; // default is 0, but this might change in future

  if(xstart == NULL && inxstart == NULL && typicalX == NULL)
  {
    ERRORMESSAGE("Warning: initialX undefined. typicalX = 0.5...0.5 used.");
    printf("\nWarning: initialX undefined. typicalX = 0.5...0.5 used.\n");
  }
  if(rgInitialStds == NULL && inrgsigma == NULL)
  {
    ERRORMESSAGE("Warning: initialStandardDeviations undefined. 0.3...0.3 used.");
    printf("\nWarning: initialStandardDeviations. 0.3...0.3 used.\n");
  }

  if(!xstart)
  {
    xstart = new double[N];
    if(inxstart)
    {
      for(int i = 0; i < N; ++i)
        xstart[i] = inxstart[i];
    }
    else if(typicalX)
    {
      typicalXcase = true;
      for(int i = 0; i < N; ++i)
        xstart[i] = typicalX[i];
    }
    else
    {
      typicalXcase = true;
      for(int i = 0; i < N; i++)
        xstart[i] = 0.5;
    }
  }

  if(!rgInitialStds)
  {
    rgInitialStds = new double[N];
    if(inrgsigma)
    {
      for(int i = 0; i < N; ++i)
        rgInitialStds[i] = inrgsigma[i];
    }
    else
    {
      for(int i = 0; i < N; ++i)
        rgInitialStds[i] = 0.3;
    }
  }

  supplementDefaults();
}

Parameters::~Parameters()
{
  if(xstart)
    delete[] xstart;
  if(typicalX)
    delete[] typicalX;
  if(rgInitialStds)
    delete[] rgInitialStds;
  if(rgDiffMinChange)
    delete[] rgDiffMinChange;
  if(weights)
    delete[] weights;
}

void Parameters::assign(const Parameters& p)
{
  N = p.N;

  if(xstart)
    delete[] xstart;
  if(p.xstart)
  {
    xstart = new double[N];
    for(int i = 0; i < N; i++)
      xstart[i] = p.xstart[i];
  }

  if(typicalX)
    delete[] typicalX;
  if(p.typicalX)
  {
    typicalX = new double[N];
    for(int i = 0; i < N; i++)
      typicalX[i] = p.typicalX[i];
  }

  typicalXcase = p.typicalXcase;

  if(rgInitialStds)
    delete[] rgInitialStds;
  if(p.rgInitialStds)
  {
    rgInitialStds = new double[N];
    for(int i = 0; i < N; i++)
      rgInitialStds[i] = p.rgInitialStds[i];
  }

  if(rgDiffMinChange)
    delete[] rgDiffMinChange;
  if(p.rgDiffMinChange)
  {
    rgDiffMinChange = new double[N];
    for(int i = 0; i < N; i++)
      rgDiffMinChange[i] = p.rgDiffMinChange[i];
  }

  stopMaxFunEvals = p.stopMaxFunEvals;
  facmaxeval = p.facmaxeval;
  stopMaxIter = p.stopMaxIter;

  stStopFitness.flg = p.stStopFitness.flg;
  stStopFitness.val = p.stStopFitness.val;

  stopTolFun = p.stopTolFun;
  stopTolFunHist = p.stopTolFunHist;
  stopTolX = p.stopTolX;
  stopTolUpXFactor = p.stopTolUpXFactor;

  lambda = p.lambda;
  mu = p.mu;
  mucov = p.mucov;
  mueff = p.mueff;

  if(weights)
    delete[] weights;
  if(p.weights)
  {
    weights = new double[mu];
    for(int i = 0; i < mu; i++)
      weights[i] = p.weights[i];
  }

  damps = p.damps;
  cs = p.cs;
  ccumcov = p.ccumcov;
  ccov = p.ccov;
  diagonalCov = p.diagonalCov;

  updateCmode.modulo = p.updateCmode.modulo;
  updateCmode.maxtime = p.updateCmode.maxtime;

  facupdateCmode = p.facupdateCmode;

  weightMode = p.weightMode;

  resumefile = p.resumefile;
}

void Parameters::supplementDefaults()
{
  if(lambda < 2)
    lambda = 4 + (int) (3.0*log((double) N));
  if(mu <= 0)
    mu = lambda / 2;
  if(!weights)
    setWeights(weightMode);

  if(cs > 0)
    cs *= (mueff + 2.) / (N + mueff + 3.);
  if(cs <= 0 || cs >= 1)
    cs = (mueff + 2.) / (N + mueff + 3.);

  if(ccumcov <= 0 || ccumcov > 1)
    ccumcov = 4. / (N + 4);

  if(mucov < 1)
    mucov = mueff;
  double t1 = 2. / ((N + 1.4142)*(N + 1.4142));
  double t2 = (2.* mueff - 1.) / ((N + 2.)*(N + 2.) + mueff);
  t2 = (t2 > 1) ? 1 : t2;
  t2 = (1. / mucov)* t1 + (1. - 1. / mucov)* t2;
  if(ccov >= 0)
    ccov *= t2;
  if(ccov < 0 || ccov > 1)
    ccov = t2;

  if(diagonalCov < 0)
    diagonalCov = 2 + 100. * N / sqrt((double) lambda);

  if(stopMaxFunEvals <= 0)
    stopMaxFunEvals = facmaxeval * 900 * (N + 3)*(N + 3);
  else
    stopMaxFunEvals *= facmaxeval;

  if(stopMaxIter <= 0)
    stopMaxIter = ceil((double) (stopMaxFunEvals / lambda));

  if(damps < 0)
    damps = 1;
  damps = damps
     * (1 + 2* std::max(0., sqrt((mueff - 1.) / (N + 1.)) - 1))
     * std::max(0.3, 1. - // modify for short runs
      (double) N / (1e-6 + std::min(stopMaxIter, stopMaxFunEvals / lambda)))
      + cs;

  if(updateCmode.modulo < 0)
    updateCmode.modulo = 1. / ccov / (double) N / 10.;
  updateCmode.modulo *= facupdateCmode;
  if(updateCmode.maxtime < 0)
    updateCmode.maxtime = 0.20; // maximal 20% of CPU-time
}

void Parameters::setWeights(Weights mode)
{
  if(weights)
    delete[] weights;
  weights = new double[mu];
  switch(mode)
  {
  case LINEAR_WEIGHTS:
    for(int i = 0; i < mu; ++i) weights[i] = mu - i;
    break;
  case EQUAL_WEIGHTS:
    for(int i = 0; i < mu; ++i) weights[i] = 1;
    break;
  case LOG_WEIGHTS:
  default:
    for(int i = 0; i < mu; ++i) weights[i] = log(mu + 1.) - log(i + 1.);
    break;
  }

  // normalize weights vector and set mueff
  double s1 = 0, s2 = 0;
  for(int i = 0; i < mu; ++i)
  {
    s1 += weights[i];
    s2 += weights[i]*weights[i];
  }
  mueff = s1*s1/s2;
  for(int i = 0; i < mu; ++i)
    weights[i] /= s1;

  if(mu < 1 || mu > lambda || (mu == lambda && weights[0] == weights[mu - 1]))
    FATAL("setWeights(): invalid setting of mu or lambda");
}

std::ostream& operator<<(std::ostream& stream, const Parameters& p)
{
  stream << "N= " << p.N;

  if(p.xstart)
  {
    stream << ", xstart= [";
    for(int i = 0; i < p.N; i++)
    {
      stream << p.xstart[i];
      if(i != p.N-1)
        stream << ", ";
    }
    stream << "]";
  }

  if(p.typicalX)
  {
    stream << ", typicalX= [";
    for(int i = 0; i < p.N; i++)
    {
      stream << p.typicalX[i];
      if(i != p.N-1)
        stream << ", ";
    }
    stream << "]";
  }

  stream << ", typicalXcase= " << p.typicalXcase;

  if(p.rgInitialStds)
  {
    stream << ", rgInitialStds= [";
    for(int i = 0; i < p.N; i++)
    {
      stream << p.rgInitialStds[i];
      if(i != p.N-1)
        stream << ", ";
    }
    stream << "]";
  }

  if(p.rgDiffMinChange)
  {
    stream << ", rgDiffMinChange= [";
    for(int i = 0; i < p.N; i++)
    {
      stream << p.rgDiffMinChange[i];
      if(i != p.N-1)
        stream << ", ";
    }
    stream << "]";
  }

  stream << ", stopMaxFunEvals= " << p.stopMaxFunEvals;
  stream << ", facmaxeval= " << p.facmaxeval;
  stream << ", stopMaxIter= " << p.stopMaxIter;

  stream << ", stStopFitness.flg= " << p.stStopFitness.flg;
  stream << ", stStopFitness.val= " << p.stStopFitness.val;

  stream << ", stopTolFun= " << p.stopTolFun;
  stream << ", stopTolFunHist= " << p.stopTolFunHist;
  stream << ", stopTolX= " << p.stopTolX;
  stream << ", stopTolUpXFactor= " << p.stopTolUpXFactor;

  stream << ", lambda= " << p.lambda;
  stream << ", mu= " << p.mu;
  stream << ", mucov= " << p.mucov;
  stream << ", mueff= " << p.mueff;

  if(p.weights)
  {
    stream << ", weights= [";
    for(int i = 0; i < p.mu; i++)
    {
      stream << p.weights[i];
      if(i != p.mu-1)
        stream << ", ";
    }
    stream << "]";
  }

  stream << ", damps= " << p.damps;
  stream << ", cs= " << p.cs;
  stream << ", ccumcov= " << p.ccumcov;
  stream << ", ccov= " << p.ccov;
  stream << ", diagonalCov= " << p.diagonalCov;

  stream << ", p.updateCmode.modulo= " << p.updateCmode.modulo;
  stream << ", p.updateCmode.maxtime= " << p.updateCmode.maxtime;

  stream << ", p.facupdateCmode= " << p.facupdateCmode;

  stream << ", p.weightMode= " << p.weightMode;

  stream << ", p.resumefile= " << p.resumefile;
  return stream;
}
