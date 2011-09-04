/**
 * @file cmaes.cpp
 * @author Nikolaus Hansen, ported to C++ by Alexander Fabisch
 *
 * CMA-ES for non-linear function minimization.
 */

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <limits>
#include "cmaes.h"
#include "utils.h"

CMAES::~CMAES()
{
  delete[] rgpc;
  delete[] rgps;
  delete[] rgdTmp;
  delete[] rgBDz;
  delete[] --xmean;
  delete[] --rgxold;
  delete[] --rgxbestever;
  delete[] --rgout;
  delete[] rgD;
  for(int i = 0; i < sp.N; ++i)
  {
    delete[] C[i];
    delete[] B[i];
  }
  for(int i = 0; i < sp.lambda; ++i)
    delete[] --rgrgx[i];
  delete[] rgrgx;
  delete[] C;
  delete[] B;
  delete[] index;
  delete[] publicFitness;
  delete[] --rgFuncValue;
  delete[] --arFuncValueHist;
}

std::string CMAES::getTimeStr(void)
{
  time_t tm = time(NULL);
  static char s[33];

  strncpy(s, ctime(&tm), 24);
  s[24] = '\0'; // cut the \n
  return std::string(s);
}

double* CMAES::init(const Parameters& parameters)
{
  sp = parameters;

  version = "1.0alpha";
  stopMessage = "";

  double trace = 0.;
  for(int i = 0; i < sp.N; ++i)
    trace += sp.rgInitialStds[i]*sp.rgInitialStds[i];
  sigma = sqrt(trace/sp.N);

  chiN = sqrt((double) sp.N)* (1. - 1./(4.*sp.N) + 1./(21.*sp.N*sp.N));
  eigensysIsUptodate = true;
  doCheckEigen = false;
  genOfEigensysUpdate = 0;
  isResumeDone = false;

  double dtest;
  for(dtest = 1.; dtest && dtest < 1.1*dtest; dtest *= 2.)
    if(dtest == dtest + 1.)
      break;
  dMaxSignifKond = dtest / 1000.; // not sure whether this is really save, 100 does not work well enough

  gen = 0;
  countevals = 0;
  state = INITIALIZED;
  dLastMinEWgroesserNull = 1.0;
  printtime = writetime = firstwritetime = firstprinttime = 0;

  rgpc = new double[sp.N];
  rgps = new double[sp.N];
  rgdTmp = new double[sp.N+1];
  rgBDz = new double[sp.N];
  xmean = new double[sp.N+2];
  xmean[0] = sp.N;
  ++xmean;
  rgxold = new double[sp.N+2];
  rgxold[0] = sp.N;
  ++rgxold;
  rgxbestever = new double[sp.N+3];
  rgxbestever[0] = sp.N;
  ++rgxbestever;
  rgxbestever[sp.N] = std::numeric_limits<double>::max();
  rgout = new double[sp.N+2];
  rgout[0] = sp.N;
  ++rgout;
  rgD = new double[sp.N];
  C = new double*[sp.N];
  B = new double*[sp.N];
  publicFitness = new double[sp.lambda];
  rgFuncValue = new double[sp.lambda+1];
  rgFuncValue[0] = sp.lambda;
  ++rgFuncValue;
  const int historySize = 10 + (int) ceil(3.*10.*sp.N/sp.lambda);
  arFuncValueHist = new double[historySize + 1];
  arFuncValueHist[0] = (double) historySize;
  arFuncValueHist++;

  for(int i = 0; i < sp.N; ++i)
  {
    C[i] = new double[i+1];
    B[i] = new double[sp.N];
  }
  index = new int[sp.lambda];
  for(int i = 0; i < sp.lambda; ++i)
      index[i] = i;
  rgrgx = new double*[sp.lambda];
  for(int i = 0; i < sp.lambda; ++i)
  {
    rgrgx[i] = new double[sp.N+2];
    rgrgx[i][0] = sp.N;
    rgrgx[i]++;
    for(int j = 0; j < sp.N; j++)
      rgrgx[i][j] = 0.0;
  }

  // initialize newed space
  for(int i = 0; i < sp.lambda; i++)
  {
    rgFuncValue[i] = std::numeric_limits<double>::max();
  }
  for(int i = 0; i < historySize; i++)
  {
    arFuncValueHist[i] = std::numeric_limits<double>::max();
  }
  for(int i = 0; i < sp.N; ++i)
    for(int j = 0; j < i; ++j)
      C[i][j] = B[i][j] = B[j][i] = 0.;

  for(int i = 0; i < sp.N; ++i)
  {
    B[i][i] = 1.;
    C[i][i] = rgD[i] = sp.rgInitialStds[i]*sqrt(sp.N/trace);
    C[i][i] *= C[i][i];
    rgpc[i] = rgps[i] = 0.;
  }
  minEW = minElement(rgD, sp.N);
  minEW = minEW*minEW;
  maxEW = maxElement(rgD, sp.N);
  maxEW = maxEW*maxEW;

  maxdiagC = C[0][0];
  for(int i = 1; i < sp.N; ++i) if(maxdiagC < C[i][i]) maxdiagC = C[i][i];
  mindiagC = C[0][0];
  for(int i = 1; i < sp.N; ++i) if(mindiagC > C[i][i]) mindiagC = C[i][i];

  for(int i = 0; i < sp.N; ++i)
    xmean[i] = rgxold[i] = sp.xstart[i];
  // use in case xstart as typicalX
  if(sp.typicalXcase)
    for(int i = 0; i < sp.N; ++i)
      xmean[i] += sigma*rgD[i]*rand.gauss();

  if(sp.resumefile != "")
    resumeDistribution(sp.resumefile);

  return publicFitness;
}

std::string CMAES::sayHello()
{
  std::stringstream stream;
  stream << "(" << sp.mu << "," << sp.lambda << ")-CMA-ES(mu_eff="
      << std::setprecision(1) << sp.mueff << "), Ver=\"" << version
      << "\", dimension=" << sp.N << ", diagonalIterations="
      << (long) sp.diagonalCov << " (" << getTimeStr() << ")";
  return stream.str();
}

void CMAES::resumeDistribution(const std::string& filename)
{
  int i, j, res, n;
  double d;
  FILE *fp = fopen(filename.c_str(), "r");
  if(fp == NULL)
  {
    ERRORMESSAGE("resumeDistribution(): could not open '" + filename + "'");
    return;
  }
  // count number of "resume" entries
  i = 0;
  res = 0;
  while(1)
  {
    if((res = fscanf(fp, " resume %lg", &d)) == EOF)
      break;
    else if(res == 0)
      fscanf(fp, " %*s");
    else if(res > 0)
      i += 1;
  }

  // go to last "resume" entry
  n = i;
  i = 0;
  res = 0;
  rewind(fp);
  while(i < n)
  {
    if((res = fscanf(fp, " resume %lg", &d)) == EOF)
      FATAL("resumeDistribution(): Unexpected error, bug");
    else if(res == 0)
      fscanf(fp, " %*s");
    else if(res > 0)
      ++i;
  }
  if(d != sp.N)
    FATAL("resumeDistribution(): Dimension numbers do not match");

  // find next "xmean" entry
  while(1)
  {
    if((res = fscanf(fp, " xmean %lg", &d)) == EOF)
      FATAL("resumeDistribution(): 'xmean' not found");
    else if(res == 0)
      fscanf(fp, " %*s");
    else if(res > 0)
      break;
  }

  // read xmean
  xmean[0] = d;
  res = 1;
  for(i = 1; i < sp.N; ++i)
    res += fscanf(fp, " %lg", &xmean[i]);
  if(res != sp.N)
    FATAL("resumeDistribution(): xmean: dimensions differ");

  // find next "path for sigma" entry
  while(1)
  {
    if((res = fscanf(fp, " path for sigma %lg", &d)) == EOF)
      FATAL("resumeDistribution(): 'path for sigma' not found");
    else if(res == 0)
      fscanf(fp, " %*s");
    else if(res > 0)
      break;
  }

  // read ps
  rgps[0] = d;
  res = 1;
  for(i = 1; i < sp.N; ++i)
    res += fscanf(fp, " %lg", &rgps[i]);
  if(res != sp.N)
    FATAL("resumeDistribution(): ps: dimensions differ");

  // find next "path for C" entry
  while(1)
  {
    if((res = fscanf(fp, " path for C %lg", &d)) == EOF)
      FATAL("resumeDistribution(): 'path for C' not found");
    else if(res == 0)
      fscanf(fp, " %*s");
    else if(res > 0)
      break;
  }
  // read pc
  rgpc[0] = d;
  res = 1;
  for(i = 1; i < sp.N; ++i)
    res += fscanf(fp, " %lg", &rgpc[i]);
  if(res != sp.N)
    FATAL("resumeDistribution(): pc: dimensions differ");

  // find next "sigma" entry
  while(1)
  {
    if((res = fscanf(fp, " sigma %lg", &d)) == EOF)
      FATAL("resumeDistribution(): 'sigma' not found");
    else if(res == 0)
      fscanf(fp, " %*s");
    else if(res > 0)
      break;
  }
  sigma = d;

  // find next entry "covariance matrix"
  while(1)
  {
    if((res = fscanf(fp, " covariance matrix %lg", &d)) == EOF)
      FATAL("resumeDistribution(): 'covariance matrix' not found");
    else if(res == 0)
      fscanf(fp, " %*s");
    else if(res > 0)
      break;
  }
  // read C
  C[0][0] = d;
  res = 1;
  for(i = 1; i < sp.N; ++i)
    for(j = 0; j <= i; ++j)
      res += fscanf(fp, " %lg", &C[i][j]);
  if(res != (sp.N* sp.N + sp.N) / 2)
    FATAL("resumeDistribution(): C: dimensions differ");

  eigensysIsUptodate = false;
  isResumeDone = true;
  updateEigensystem(true);
}

double const* CMAES::setMean(const double *newxmean)
{
  if(state == SAMPLED)
    FATAL("setMean: mean cannot be set inbetween the calls of samplePopulation"
        " and updateDistribution");

  if(newxmean != NULL && newxmean != xmean)
    for(int i = 0; i < sp.N; ++i)
      xmean[i] = newxmean[i];
  else
    newxmean = xmean;

  return newxmean;
}

void CMAES::addMutation(double* x, double eps)
{
  for(int i = 0; i < sp.N; ++i)
    rgdTmp[i] = rgD[i]*rand.gauss();
  for(int i = 0; i < sp.N; ++i)
  {
    double sum = 0.0;
    for(int j = 0; j < sp.N; ++j)
      sum += B[i][j]*rgdTmp[j];
    x[i] = xmean[i] + eps*sigma*sum;
  }
}

double* const* CMAES::samplePopulation()
{
  bool diag = sp.diagonalCov == 1 || sp.diagonalCov >= gen;

  // calculate eigensystem
  if(!eigensysIsUptodate)
  {
    if(!diag)
      updateEigensystem(false);
    else
    {
      for(int i = 0; i < sp.N; ++i)
        rgD[i] = sqrt(C[i][i]);
      minEW = square(minElement(rgD, sp.N));
      maxEW = square(maxElement(rgD, sp.N));
      eigensysIsUptodate = true;
      eigenTimings.start();
    }
  }

  testMinStdDevs();

  for(int iNk = 0; iNk < sp.lambda; ++iNk)
  { // generate scaled random vector D*z
    double* rgrgxink = rgrgx[iNk];
    for(int i = 0; i < sp.N; ++i)
      if(diag)
        rgrgxink[i] = xmean[i] + sigma*rgD[i]*rand.gauss();
      else
        rgdTmp[i] = rgD[i]*rand.gauss();
    if(!diag)
      for(int i = 0; i < sp.N; ++i) // add mutation sigma*B*(D*z)
      {
        double sum = 0.0;
        for(int j = 0; j < sp.N; ++j)
          sum += B[i][j]*rgdTmp[j];
        rgrgxink[i] = xmean[i] + sigma*sum;
      }
  }

  if(state == UPDATED || gen == 0)
    ++gen;
  state = SAMPLED;

  return rgrgx;
}

double const* CMAES::reSampleSingleOld(double *x)
{
  if(x == NULL)
    FATAL("reSampleSingle(): Missing input double *x");
  addMutation(x);
  return x;
}

double* const* CMAES::reSampleSingle(int iindex)
{
  double *x;
  static char s[99];
  if(iindex < 0 || iindex >= sp.lambda)
  {
    sprintf(s, "index==%d must be between 0 and %d", iindex, sp.lambda);
    FATAL("reSampleSingle(): Population member ");
  }
  x = rgrgx[iindex];
  addMutation(x);
  return rgrgx;
}

double* CMAES::sampleSingleInto(double *x)
{
  if(x == NULL)
    x = new double[sp.N];
  addMutation(x);
  return x;
}

double* CMAES::perturbSolutionInto(double *x, double const *pxmean, double eps)
{
  if(x == NULL)
    x = new double[sp.N];
  if(pxmean == NULL)
    FATAL("perturbSolutionInto(): xmean was not given");
  addMutation(x, eps);
  return x;
}

double* CMAES::updateDistribution(const double* rgFunVal)
{
  const int N = sp.N;
  bool diag = sp.diagonalCov == 1 || sp.diagonalCov >= gen;

  if(state == UPDATED)
    FATAL("updateDistribution(): You need to call \n"
        "samplePopulation() before update can take place.");
  if(rgFunVal == NULL)
    FATAL("updateDistribution(): No fitness function value array input.");

  if(state == SAMPLED) // function values are delivered here
    countevals += sp.lambda;
  else
    ERRORMESSAGE("updateDistribution(): unexpected state");

  // assign function values
  for(int i = 0; i < sp.lambda; ++i)
    rgrgx[i][N] = rgFuncValue[i] = rgFunVal[i];

  // Generate index
  sortIndex(rgFunVal, index, sp.lambda);

  // Test if function values are identical, escape flat fitness
  if(rgFuncValue[index[0]] == rgFuncValue[index[(int) sp.lambda / 2]])
  {
    sigma *= exp(0.2 + sp.cs / sp.damps);
    ERRORMESSAGE("Warning: sigma increased due to equal function values\n"
        "   Reconsider the formulation of the objective function");
  }

  // update function value history
  for(int i = (int) *(arFuncValueHist - 1) - 1; i > 0; --i)
    arFuncValueHist[i] = arFuncValueHist[i - 1];
  arFuncValueHist[0] = rgFunVal[index[0]];

  // update xbestever
  if(rgxbestever[N] > rgrgx[index[0]][N] || gen == 1)
    for(int i = 0; i <= N; ++i)
    {
      rgxbestever[i] = rgrgx[index[0]][i];
      rgxbestever[N+1] = countevals;
    }

  const double sqrtmueffdivsigma = sqrt(sp.mueff) / sigma;
  // calculate xmean and rgBDz~N(0,C)
  for(int i = 0; i < N; ++i)
  {
    rgxold[i] = xmean[i];
    xmean[i] = 0.;
    for(int iNk = 0; iNk < sp.mu; ++iNk)
      xmean[i] += sp.weights[iNk]*rgrgx[index[iNk]][i];
    rgBDz[i] = sqrtmueffdivsigma*(xmean[i]-rgxold[i]);
  }

  // calculate z := D^(-1)* B^(-1)* rgBDz into rgdTmp
  for(int i = 0; i < N; ++i)
  {
    double sum;
    if(diag)
      sum = rgBDz[i];
    else
    {
      sum = 0.;
      for(int j = 0; j < N; ++j)
        sum += B[j][i]*rgBDz[j];
    }
    rgdTmp[i] = sum/rgD[i];
  }

  // cumulation for sigma (ps) using B*z
  const double sqrtFactor = sqrt(sp.cs*(2.-sp.cs));
  const double invps = 1.-sp.cs;
  for(int i = 0; i < N; ++i)
  {
    double sum;
    if(diag)
      sum = rgdTmp[i];
    else
    {
      sum = 0.;
      double* Bi = B[i];
      for(int j = 0; j < N; ++j)
        sum += Bi[j]*rgdTmp[j];
    }
    rgps[i] = invps*rgps[i] + sqrtFactor*sum;
  }

  // calculate norm(ps)^2
  double psxps = 0.;
  for(int i = 0; i < N; ++i)
  {
    const double& rgpsi = rgps[i];
    psxps += rgpsi*rgpsi;
  }

  // cumulation for covariance matrix (pc) using B*D*z~N(0,C)
  int hsig = sqrt(psxps) / sqrt(1. - pow(1. - sp.cs, 2.* gen)) / chiN < 1.4 + 2. / (N + 1);
  const double ccumcovinv = 1.-sp.ccumcov;
  const double hsigFactor = hsig*sqrt(sp.ccumcov*(2.-sp.ccumcov));
  for(int i = 0; i < N; ++i)
    rgpc[i] = ccumcovinv* rgpc[i] + hsigFactor*rgBDz[i];

  // update of C
  adaptC2(hsig);

  // update of sigma
  sigma *= exp(((sqrt(psxps) / chiN) - 1.)* sp.cs / sp.damps);

  state = UPDATED;
  return xmean;
}

void CMAES::adaptC2(const int hsig) {
  const int N = sp.N;
  bool diag = sp.diagonalCov == 1 || sp.diagonalCov >= gen;

  if(sp.ccov != 0.)
  {
    // definitions for speeding up inner-most loop
    const double mucovinv = 1./sp.mucov;
    const double commonFactor = sp.ccov * (diag ? (N + 1.5) / 3. : 1.);
    const double ccov1 = std::min(commonFactor*mucovinv, 1.);
    const double ccovmu = std::min(commonFactor*(1.-mucovinv), 1.-ccov1);
    const double sigmasquare = sigma*sigma;
    const double onemccov1ccovmu = 1.-ccov1-ccovmu;
    const double longFactor = (1.-hsig)*sp.ccumcov*(2.-sp.ccumcov);

    eigensysIsUptodate = false;

    // update covariance matrix
    for(int i = 0; i < N; ++i)
      for(int j = diag ? i : 0; j <= i; ++j)
      {
        double& Cij = C[i][j];
        Cij = onemccov1ccovmu*Cij + ccov1 * (rgpc[i]*rgpc[j] + longFactor*Cij);
        for(int k = 0; k < sp.mu; ++k)
        { // additional rank mu update
          const double* rgrgxindexk = rgrgx[index[k]];
          Cij += ccovmu*sp.weights[k] * (rgrgxindexk[i] - rgxold[i])
              * (rgrgxindexk[j] - rgxold[j]) / sigmasquare;
        }
      }
    // update maximal and minimal diagonal value
    maxdiagC = mindiagC = C[0][0];
    for(int i = 1; i < N; ++i)
    {
      const double& Cii = C[i][i];
      if(maxdiagC < Cii)
        maxdiagC = Cii;
      else if(mindiagC > Cii)
        mindiagC = Cii;
    }
  }
}

void CMAES::testMinStdDevs()
{
  int i, N = this->sp.N;
  if(this->sp.rgDiffMinChange == NULL)
    return;

  for(i = 0; i < N; ++i)
    while(this->sigma* sqrt(this->C[i][i]) < this->sp.rgDiffMinChange[i])
      this->sigma *= exp(0.05 + this->sp.cs / this->sp.damps);
}

void CMAES::writeToFile(int key, const std::string& filename)
{
  std::ofstream file;
  file.open(filename.c_str(), std::ios_base::app);

  if(file.is_open())
  {
    if(gen > 0 || filename.substr(0, 11) != "outcmaesfit")
      writeToStream(key, file); /* do not write fitness for gen==0 */
    file.close();
  }
  else
  {
    ERRORMESSAGE("writeToFile(): could not open '" + filename + "'");
    return;
  }
}

void CMAES::writeToStream(int key, std::ostream& file)
{
  if(key & WKResume)
  {
    file << std::endl << "# resume " << sp.N << std::endl;
    file << "xmean" << std::endl;
    writeToStream(WKXMean, file);
    file << "path for sigma" << std::endl;
    for(int i = 0; i < sp.N; ++i)
      file << rgps[i] << (i == sp.N-1 ? "\n" : "\t");
    file << "path for C" << std::endl;
    for(int i = 0; i < sp.N; ++i)
      file << rgpc[i] << (i == sp.N-1 ? "\n" : "\t");
    file << "sigma " << sigma << std::endl;
    // note than B and D might not be up-to-date
    file << "covariance matrix" << std::endl;
    writeToStream(WKC, file);
  }
  if(key & WKXMean)
  {
    for(int i = 0; i < sp.N; ++i)
      file << (i == 0 ? "" : "\t") << xmean[i];
    file << std::endl;
  }
  if(key & WKC)
  {
    for(int i = 0; i < sp.N; ++i)
      for(int j = 0; j <= i; ++j)
      {
        file << C[i][j];
        if(j == i)
          file << std::endl;
        else
          file << '\t';
      }
    file << std::endl;
  }
  if(key & WKAll)
  {
    time_t ti = time(NULL);
    file << std::endl << "# --------- " << asctime(localtime(&ti)) << std::endl;
    file << " N " << sp.N << std::endl;
    file << "function evaluations " << (long) countevals << std::endl;
    file << "elapsed (CPU) time [s] " << std::setprecision(2) << eigenTimings.totaltotaltime << std::endl;
    file << "function value f(x)=" << rgrgx[index[0]][sp.N] << std::endl;
    file << "maximal standard deviation " << sigma*sqrt(maxdiagC) << std::endl;
    file << "minimal standard deviation " << sigma*sqrt(mindiagC) << std::endl;
    file << "sigma " << sigma << std::endl;
    file << "axisratio " << (maxElement(rgD, sp.N) / minElement(rgD, sp.N)) << std::endl;
    file << "xbestever found after " << std::setprecision(0) << rgxbestever[sp.N+1]
        << "evaluations, function value " << rgxbestever[sp.N] << std::endl;
    for(int i = 0; i < sp.N; ++i)
      file << " " << std::setw(12) << rgxbestever[i] << (i % 5 == 4 || i == sp.N-1 ? '\n' : ' ');
    file << "xbest (of last generation, function value " << rgrgx[index[0]][sp.N] << ")" << std::endl;
    for(int i = 0; i < sp.N; ++i)
      file << " " << std::setw(12) << rgrgx[index[0]][i] << (i % 5 == 4 || i == sp.N-1 ? '\n' : ' ');
    file << "xmean" << std::endl;
    for(int i = 0; i < sp.N; ++i)
      file << " " << std::setw(12) << xmean[i] << (i % 5 == 4 || i == sp.N-1 ? '\n' : ' ');
    file << "Standard deviation of coordinate axes (sigma*sqrt(diag(C)))" << std::endl;
    for(int i = 0; i < sp.N; ++i)
      file << " " << std::setw(12) << sigma*sqrt(C[i][i]) << (i % 5 == 4 || i == sp.N-1 ? '\n' : ' ');
    file << "Main axis lengths of mutation ellipsoid (sigma*diag(D))" << std::endl;
    for(int i = 0; i < sp.N; ++i)
        rgdTmp[i] = rgD[i];
    std::sort(rgdTmp, rgdTmp + sp.N);
    for(int i = 0; i < sp.N; ++i)
      file << " " << std::setw(12) << sigma*rgdTmp[sp.N-1-i] << (i % 5 == 4 || i == sp.N-1 ? '\n' : ' ');
    file << "Longest axis (b_i where d_ii=max(diag(D))" << std::endl;
    int k = maxIndex(rgD, sp.N);
    for(int i = 0; i < sp.N; ++i)
      file << " " << std::setw(12) << B[i][k] << (i % 5 == 4 || i == sp.N-1 ? '\n' : ' ');
    file << "Shortest axis (b_i where d_ii=max(diag(D))" << std::endl;
    k = minIndex(rgD, sp.N);
    for(int i = 0; i < sp.N; ++i)
      file << " " << std::setw(12) << B[i][k] << (i % 5 == 4 || i == sp.N-1 ? '\n' : ' ');
    file << std::endl;
  }
  if(key & WKFewInfo)
  {
    file << " Iter\tFevals\tFunction Value\tSigma\tMaxCoorDev\tMinCoorDev\t"
        << "AxisRatio\tMinDii\tTime in eig" << std::endl;
    file << std::endl;
  }
  if(key & WKFew)
  {
    file << (int) gen << "\t" << (int) countevals << "\t"
        << rgFuncValue[index[0]] << "\t\t" << sigma << "  "
        << sigma*sqrt(maxdiagC) << "\t" << sigma*sqrt(mindiagC)
        << "\t" << std::scientific << std::setprecision(2)
        << sqrt(maxEW / minEW) << "\t" << sqrt(minEW)
        << "  " << eigenTimings.totaltotaltime;
    file << std::endl;
  }
  if(key & WKEval)
  {
    file << countevals;
    file << std::endl;
  }
  if(key & WKFitness)
  {
    for(int i = 0; i < sp.N; ++i)
      file << (i == 0 ? "" : "\t") << rgFuncValue[index[i]];
    file << std::endl;
  }
  if(key & WKFBestEver)
  {
    file << rgxbestever[sp.N] << std::endl;
  }
  if(key & WKCGeneration)
  {
    file << gen << std::endl;
  }
  if(key & WKSigma)
  {
    file << sigma << std::endl;
  }
  if(key & WKLambda)
  {
    file << sp.lambda << std::endl;
  }
  if(key & WKB)
  {
    int* iindex = new int[sp.N];
    sortIndex(rgD, iindex, sp.N);
    for(int i = 0; i < sp.N; ++i)
      for(int j = 0; j < sp.N; ++j)
      {
        file << B[j][iindex[sp.N-1-i]];
        if(j != sp.N-1)
          file << '\t';
        else
          file << std::endl;
      }
    delete[] iindex;
    iindex = 0;
    file << std::endl;
  }
  if(key & WKXBest)
  {
    for(int i = 0; i < sp.N; ++i)
      file << (i == 0 ? "" : "\t") << rgrgx[index[0]][i];
    file << std::endl;
  }
  if(key & WKClock)
  {
    eigenTimings.update();
    file << eigenTimings.totaltotaltime << " " << eigenTimings.tictoctime
        << std::endl;
  }
  if(key & WKDim)
  {
    file << sp.N;
    file << std::endl;
  }
}

double CMAES::get(GetScalar key)
{
  switch(key)
  {
    case AxisRatio:
    {
      return maxElement(rgD, sp.N) / minElement(rgD, sp.N);
    }
    case Eval:
    {
      return countevals;
    }
    case Fitness:
    {
      return rgFuncValue[index[0]];
    }
    case FBestEver:
    {
      return rgxbestever[sp.N];
    }
    case Generation:
    {
      return gen;
    }
    case MaxEval:
    {
      return sp.stopMaxFunEvals;
    }
    case MaxIter:
    {
      return ceil(sp.stopMaxIter);
    }
    case MaxAxisLength:
    {
      return sigma*sqrt(maxEW);
    }
    case MinAxisLength:
    {
      return sigma*sqrt(minEW);
    }
    case MaxStdDev:
    {
      return sigma*sqrt(maxdiagC);
    }
    case MinStdDev:
    {
      return sigma*sqrt(mindiagC);
    }
    case Dimension:
    {
      return sp.N;
    }
    case SampleSize:
    {
      return sp.lambda;
    }
    case Sigma:
    {
      return sigma;
    }
    default:
      FATAL("get(): No match found for key");
      return 0;
  }
}

double* CMAES::getInto(GetVector key, double *res)
{
  double const* res0 = getPtr(key);
  if(res == NULL)
    res = new double[sp.N];
  for(int i = 0; i < sp.N; ++i)
    res[i] = res0[i];
  return res;
}

double* CMAES::getNew(GetVector key)
{
  return getInto(key, NULL);
}

const double* CMAES::getPtr(GetVector key)
{
  switch(key)
  {
    case DiagC:
    {
      for(int i = 0; i < sp.N; ++i)
        rgout[i] = C[i][i];
      return rgout;
    }
    case DiagD:
    {
      return rgD;
    }
    case StdDev:
    {
      for(int i = 0; i < sp.N; ++i)
        rgout[i] = sigma* sqrt(C[i][i]);
      return rgout;
    }
    case XBestEver:
      return rgxbestever;
    case XBest:
      return rgrgx[index[0]];
    case XMean:
      return xmean;
    default:
      FATAL("getPtr(): No match found for key");
      return 0;
  }
}

bool CMAES::testForTermination()
{
  double range, fac;
  int iAchse, iKoo;
  int diag = sp.diagonalCov == 1 || sp.diagonalCov >= gen;
  int i, cTemp, N = sp.N;
  std::stringstream message;

  if(stopMessage != "")
  {
    message << stopMessage << std::endl;
  }

  // function value reached
  if((gen > 1 || state > SAMPLED) && sp.stStopFitness.flg &&
      rgFuncValue[index[0]] <= sp.stStopFitness.val)
  {
    message << "Fitness: function value " << rgFuncValue[index[0]]
        << " <= stopFitness (" << sp.stStopFitness.val << ")" << std::endl;
  }

  // TolFun
  range = std::max(maxElement(arFuncValueHist, (int) std::min(gen, *(arFuncValueHist - 1))),
      maxElement(rgFuncValue, sp.lambda)) -
      std::min(minElement(arFuncValueHist, (int) std::min(gen, *(arFuncValueHist - 1))),
      minElement(rgFuncValue, sp.lambda));

  if(gen > 0 && range <= sp.stopTolFun)
  {
    message << "TolFun: function value differences " << range
        << " < stopTolFun=" << sp.stopTolFun << std::endl;
  }

  // TolFunHist
  if(gen > *(arFuncValueHist - 1))
  {
    range = maxElement(arFuncValueHist, (int) *(arFuncValueHist - 1))
        - minElement(arFuncValueHist, (int) *(arFuncValueHist - 1));
    if(range <= sp.stopTolFunHist)
      message << "TolFunHist: history of function value changes " << range
          << " stopTolFunHist=" << sp.stopTolFunHist << std::endl;
  }

  // TolX
  for(i = 0, cTemp = 0; i < N; ++i)
  {
    cTemp += (sigma* sqrt(C[i][i]) < sp.stopTolX) ? 1 : 0;
    cTemp += (sigma* rgpc[i] < sp.stopTolX) ? 1 : 0;
  }
  if(cTemp == 2* N)
  {
    message << "TolX: object variable changes below " << sp.stopTolX << std::endl;
  }

  // TolUpX
  for(i = 0; i < N; ++i)
  {
    if(sigma* sqrt(C[i][i]) > sp.stopTolUpXFactor* sp.rgInitialStds[i])
      break;
  }
  if(i < N)
  {
    message << "TolUpX: standard deviation increased by more than "
        << sp.stopTolUpXFactor << ", larger initial standard deviation recommended \n"
        << std::endl;
  }

  // Condition of C greater than dMaxSignifKond
  if(maxEW >= minEW* dMaxSignifKond)
  {
    message << "ConditionNumber: maximal condition number " << dMaxSignifKond
        << " reached. maxEW=" << maxEW <<  ",minEW=" << minEW << ",maxdiagC="
        << maxdiagC << ",mindiagC=" << mindiagC << std::endl;
  }

  // Principal axis i has no effect on xmean, ie. x == x + 0.1* sigma* rgD[i]* B[i]
  if(!diag)
  {
    for(iAchse = 0; iAchse < N; ++iAchse)
    {
      fac = 0.1* sigma* rgD[iAchse];
      for(iKoo = 0; iKoo < N; ++iKoo)
      {
        if(xmean[iKoo] != xmean[iKoo] + fac* B[iKoo][iAchse])
          break;
      }
      if(iKoo == N)
      {
        message << "NoEffectAxis: standard deviation 0.1*" << (fac / 0.1)
            << " in principal axis " << iAchse << " without effect" << std::endl;
        break;
      }
    }
  }
  // Component of xmean is not changed anymore
  for(iKoo = 0; iKoo < N; ++iKoo)
  {
    if(xmean[iKoo] == xmean[iKoo] + 0.2* sigma* sqrt(C[iKoo][iKoo]))
    {
      message << "NoEffectCoordinate: standard deviation 0.2*"
          << (sigma*sqrt(C[iKoo][iKoo])) << " in coordinate " << iKoo
          << " without effect" << std::endl;
      break;
    }
  }

  if(countevals >= sp.stopMaxFunEvals)
  {
    message << "MaxFunEvals: conducted function evaluations " << countevals
        << " >= " << sp.stopMaxFunEvals << std::endl;
  }
  if(gen >= sp.stopMaxIter)
  {
    message << "MaxIter: number of iterations " << gen << " >= "
        << sp.stopMaxIter << std::endl;
  }

  stopMessage = message.str();
  return stopMessage != "";
}

std::string CMAES::getStopMessage()
{
  return stopMessage;
}

int CMAES::checkEigen(double* diag, double** Q)
{
  // compute Q diag Q^T and Q Q^T to check
  int res = 0;
  static char s[324];
  for(int i = 0; i < sp.N; ++i)
    for(int j = 0; j < sp.N; ++j) {
      double cc = 0., dd = 0.;
      for(int k = 0; k < sp.N; ++k)
      {
        cc += diag[k]*Q[i][k]*Q[j][k];
        dd += Q[i][k]*Q[j][k];
      }
      // check here, is the normalization the right one?
      const bool cond1 = fabs(cc - C[i > j ? i : j][i > j ? j : i]) / sqrt(C[i][i]* C[j][j]) > 1e-10;
      const bool cond2 = fabs(cc - C[i > j ? i : j][i > j ? j : i]) > 3e-14;
      if(cond1 && cond2)
      {
        sprintf(s, "%d %d: %.17e %.17e, %e",
                i, j, cc, C[i > j ? i : j][i > j ? j : i], cc - C[i > j ? i : j][i > j ? j : i]);
        ERRORMESSAGE("Eigen(): imprecise result detected " + std::string(s));
        ++res;
      }
      if(fabs(dd - (i == j)) > 1e-10)
      {
        sprintf(s, "%d %d %.17e ", i, j, dd);
        ERRORMESSAGE("Eigen(): imprecise result detected (Q not orthog.)" + std::string(s));
        ++res;
      }
    }
  return res;
}

void CMAES::updateEigensystem(bool force)
{
  eigenTimings.update();

  if(!force)
  {
    if(eigensysIsUptodate)
      return;
    // return on modulo generation number
    if(gen < genOfEigensysUpdate + sp.updateCmode.modulo)
      return;
    // return on time percentage
    if(sp.updateCmode.maxtime < 1.00
        && eigenTimings.tictoctime > sp.updateCmode.maxtime* eigenTimings.totaltime
        && eigenTimings.tictoctime > 0.0002)
      return;
  }

  eigenTimings.tic();
  Eigen(rgD, B, rgdTmp);
  eigenTimings.toc();

  // find largest and smallest eigenvalue, they are supposed to be sorted anyway
  minEW = minElement(rgD, sp.N);
  maxEW = maxElement(rgD, sp.N);

  if(doCheckEigen) // needs O(n^3)! writes, in case, error message in error file
    checkEigen(rgD, B);

  for(int i = 0; i < sp.N; ++i)
    rgD[i] = sqrt(rgD[i]);

  eigensysIsUptodate = true;
  genOfEigensysUpdate = gen;
}

void CMAES::Eigen(double *diag, double **Q, double *rgtmp)
{
  if(rgtmp == NULL)
    FATAL("eigen(): input parameter double *rgtmp must be non-NULL");

  if(C != Q) // copy C to Q
  {
    for(int i = 0; i < sp.N; ++i)
      for(int j = 0; j <= i; ++j)
        Q[i][j] = Q[j][i] = C[i][j];
  }

  householder(Q, diag, rgtmp);
  ql(diag, rgtmp, Q);
}

void CMAES::ql(double *d, double *e, double **V)
{
  const int n = sp.N;
  double f = 0.0;
  double tst1 = 0.0;
  const double eps = 2.22e-16; // 2.0^-52.0 = 2.22e-16

  // shift input e
  double* ep1 = e;
  for(double *ep2 = e+1, *const end = e+n; ep2 != end; ep1++, ep2++)
    *ep1 = *ep2;
  *ep1 = 0.0; // never changed again

  for(int l = 0; l < n; l++)
  {
    // find small subdiagonal element
    double& el = e[l];
    double& dl = d[l];
    const double smallSDElement = fabs(dl) + fabs(el);
    if(tst1 < smallSDElement)
      tst1 = smallSDElement;
    const double epsTst1 = eps*tst1;
    int m = l;
    while(m < n)
    {
      if(fabs(e[m]) <= epsTst1) break;
      m++;
    }

    // if m == l, d[l] is an eigenvalue, otherwise, iterate.
    if(m > l)
    {
      do {
        double h, g = dl;
        double& dl1r = d[l+1];
        double p = (dl1r - g) / (2.0*el);
        double r = myhypot(p, 1.);

        // compute implicit shift
        if(p < 0) r = -r;
        const double pr = p+r;
        dl = el/pr;
        h = g - dl;
        const double dl1 = el*pr;
        dl1r = dl1;
        for(int i = l+2; i < n; i++) d[i] -= h;
        f += h;

        // implicit QL transformation.
        p = d[m];
        double c = 1.0;
        double c2 = 1.0;
        double c3 = 1.0;
        const double el1 = e[l+1];
        double s = 0.0;
        double s2 = 0.0;
        for(int i = m-1; i >= l; i--)
        {
          c3 = c2;
          c2 = c;
          s2 = s;
          const double& ei = e[i];
          g = c*ei;
          h = c*p;
          r = myhypot(p, ei);
          e[i+1] = s*r;
          s = ei/r;
          c = p/r;
          const double& di = d[i];
          p = c*di - s*g;
          d[i+1] = h + s*(c*g + s*di);

          // accumulate transformation.
          for(int k = 0; k < n; k++)
          {
            double& Vki1 = V[k][i+1];
            h = Vki1;
            double& Vki = V[k][i];
            Vki1 = s*Vki + c*h;
            Vki *= c; Vki -= s*h;
          }
        }
        p = -s*s2*c3*el1*el/dl1;
        el = s*p;
        dl = c*p;
      } while(fabs(el) > epsTst1);
    }
    dl += f;
    el = 0.0;
  }
}

void CMAES::householder(double **V, double *d, double *e)
{
  const int n = sp.N;

  for(int j = 0; j < n; j++)
  {
    d[j] = V[n - 1][j];
  }

  // Householder reduction to tridiagonal form

  for(int i = n - 1; i > 0; i--)
  {
    // scale to avoid under/overflow
    double scale = 0.0;
    double h = 0.0;
    for(double *pd = d, *const dend = d+i; pd != dend; pd++)
    {
      scale += fabs(*pd);
    }
    if(scale == 0.0)
    {
      e[i] = d[i-1];
      for(int j = 0; j < i; j++)
      {
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
        V[j][i] = 0.0;
      }
    }
    else
    {
      // generate Householder vector
      for(double *pd = d, *const dend = d+i; pd != dend; pd++)
      {
        *pd /= scale;
        h += *pd * *pd;
      }
      double& dim1 = d[i-1];
      double f = dim1;
      double g = f > 0 ? -sqrt(h) : sqrt(h);
      e[i] = scale*g;
      h = h - f* g;
      dim1 = f - g;
      memset((void *) e, 0, (size_t)i*sizeof(double));

      // apply similarity transformation to remaining columns
      for(int j = 0; j < i; j++)
      {
        f = d[j];
        V[j][i] = f;
        double& ej = e[j];
        g = ej + V[j][j]* f;
        for(int k = j + 1; k <= i - 1; k++)
        {
          double& Vkj = V[k][j];
          g += Vkj*d[k];
          e[k] += Vkj*f;
        }
        ej = g;
      }
      f = 0.0;
      for(int j = 0; j < i; j++)
      {
        double& ej = e[j];
        ej /= h;
        f += ej* d[j];
      }
      double hh = f / (h + h);
      for(int j = 0; j < i; j++)
      {
        e[j] -= hh*d[j];
      }
      for(int j = 0; j < i; j++)
      {
        double& dj = d[j];
        f = dj;
        g = e[j];
        for(int k = j; k <= i - 1; k++)
        {
          V[k][j] -= f*e[k] + g*d[k];
        }
        dj = V[i-1][j];
        V[i][j] = 0.0;
      }
    }
    d[i] = h;
  }

  // accumulate transformations
  const int nm1 = n-1;
  for(int i = 0; i < nm1; i++)
  {
    double h;
    double& Vii = V[i][i];
    V[n-1][i] = Vii;
    Vii = 1.0;
    h = d[i+1];
    if(h != 0.0)
    {
      for(int k = 0; k <= i; k++)
      {
        d[k] = V[k][i+1] / h;
      }
      for(int j = 0; j <= i; j++) {
        double g = 0.0;
        for(int k = 0; k <= i; k++)
        {
          double* Vk = V[k];
          g += Vk[i+1]* Vk[j];
        }
        for(int k = 0; k <= i; k++)
        {
          V[k][j] -= g*d[k];
        }
      }
    }
    for(int k = 0; k <= i; k++)
    {
      V[k][i+1] = 0.0;
    }
  }
  for(int j = 0; j < n; j++)
  {
    double& Vnm1j = V[n-1][j];
    d[j] = Vnm1j;
    Vnm1j = 0.0;
  }
  V[n-1][n-1] = 1.0;
  e[0] = 0.0;
}

void CMAES::sortIndex(const double *rgFunVal, int *iindex, int n)
{
  int i, j;
  for(i = 1, iindex[0] = 0; i < n; ++i)
  {
    for(j = i; j > 0; --j)
    {
      if(rgFunVal[iindex[j - 1]] < rgFunVal[i])
        break;
      iindex[j] = iindex[j - 1]; // shift up
    }
    iindex[j] = i;
  }
}

